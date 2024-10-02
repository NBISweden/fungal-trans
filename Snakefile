import pandas as pd
import glob
import os
import sys
import platform
from snakemake.utils import min_version, validate
from snakemake.exceptions import WorkflowError
from source.utils.list_jgi_genomes import parse_genomes

include: "source/utils/common.py"

def parse_validation_error(e):
    """
    Catches errors thrown when validating the config file and attempts to
    print more meaningful messages.

    :param e: Error
    :return:
    """
    instance = ""
    print("ERROR VALIDATING CONFIG FILE")
    for item in str(e).split("\n"):
        item = item.replace('"','')
        if "ValidationError:" in item:
            print(item)
        if item[0:11] == "On instance":
            instance = item.replace("On instance['", "INCORRECT CONFIG AT: ")
            instance = instance.replace("']:","")
            print(instance)
    return

# Set temporary dir
if not os.getenv("TMPDIR"):
    os.environ["TMPDIR"] = "tmp"
    os.makedirs(os.environ["TMPDIR"],exist_ok=True)

# Set default config and validate against schema
if os.path.exists("config.yml"):
    configfile: "config.yml"
try:
    validate(config, "config/config.schema.yaml")
except WorkflowError as e:
    parse_validation_error(e)
    sys.exit()

workdir: config["workdir"]

# Parse samples and assemblies
samples, map_dict, assemblies = parse_sample_list(config["sample_file_list"], config)
# Parse JGI genomes from fungi info url, also save to file for quick re-use in future runs
genomes = parse_genomes(config["fungi_info"], f="resources/JGI/genomes.tsv")

wildcard_constraints:
    sample_id = f"({'|'.join(list(samples.keys()))})",
    assembler = "megahit|trinity|transabyss",
    filter_source = "unfiltered|filtered",
    portals = f"({'|'.join(list(genomes.index.tolist()))})",
    taxname = f"({'|'.join(list(config['taxmap'].keys()))})",
    aligner = "bowtie2|star",

# Get environment info
pythonpath = sys.executable
envdir = '/'.join(pythonpath.split("/")[0:-2])
system = platform.system()
config["system"] = system


# Include rules
include: "source/rules/db.smk"
include: "source/rules/preprocess.smk"
include: "source/rules/filter.smk"
include: "source/rules/assembly.smk"
include: "source/rules/map.smk"
include: "source/rules/sourmash.smk"
include: "source/rules/annotate_single.smk"
include: "source/rules/annotate_co.smk"
include: "source/rules/kraken.smk"

# Define targets
## Preprocessing
preprocess = expand("results/preprocess/sortmerna/{sample_id}/{sample_id}_R{i}.cut.trim.mRNA.fastq.gz", sample_id = samples.keys(), i = [1,2])
preprocess += ["results/report/preprocess/preprocess_report.html"]
## Host reads
host_reads = expand("results/host/{sample_id}_R{i}.host.fastq.gz",
            sample_id = samples.keys(), i = [1,2])

## filter
filtered_reads = expand("results/{aligner}/{sample_id}/{sample_id}_R{i}.fungi.nohost.fastq.gz",
            sample_id = samples.keys(), i = [1,2], aligner = config["host_aligner"])
filtered_reads += ["results/report/filtering/filter_report.html"]
filtered_reads += host_reads

## Sourmash
sourmash = expand("results/sourmash/{sample_id}/{sample_id}.{source}.sig",
            sample_id = samples.keys(), source = config["read_source"])
sourmash += expand("results/report/sourmash/{source}_sample.dist.labels.txt", source = config["read_source"])
sourmash_assembly = expand("results/assembly/{assembler}/{source}/{sample_id}/final.sig",
            assembler = config["assembler"], source = config["read_source"], sample_id=samples.keys())
sourmash_assembly += expand("results/report/sourmash/{source}_assembly_{assembler}.dist.labels.txt",
            source = config["read_source"], assembler = config["assembler"])

## Single-assemblies
assembly = expand("results/report/assembly/{source}_{assembler}_stats.tsv",
            source = config["read_source"], assembler = config["assembler"])
assembly += expand("results/assembly/{assembler}/{source}/{sample_id}/final.fa",
                assembler = config["assembler"], source = config["read_source"], sample_id = samples.keys())
## Annotations
taxonomy = expand("results/collated/{assembler}/{source}/taxonomy/taxonomy.{fc}.tsv",
            assembler = config["assembler"],
            fc = ["tpm","raw"],
            source = config["read_source"])
dbCAN = expand("results/collated/{assembler}/{source}/dbCAN/dbCAN.{fc}.tsv",
            assembler = config["assembler"],
            fc = ["tpm","raw"],
            source = config["read_source"])
eggnog = expand("results/collated/{assembler}/{source}/eggNOG/{db}.{fc}.tsv",
            db = ["enzymes","pathways","pathways.norm","modules","kos","tc","cazy"],
            assembler = config["assembler"],
            fc = ["tpm","raw"],
            source = config["read_source"])
mapres = expand("results/assembly/{assembler}/{source}/{sample_id}/{sample_id}_R{i}.fastq.gz",
            assembler = config["assembler"], source = config["read_source"], sample_id = samples.keys(), i = [1,2])
mapres += expand("results/report/map/{assembler}_{source}_map_report.html",
            assembler = config["assembler"], source = config["read_source"])
normalize = expand("results/annotation/{assembler}/{source}/{sample_id}/featureCounts/fc.{fc}.tab",
            assembler = config["assembler"], source = config["read_source"], sample_id = samples.keys(), fc=["tpm","raw"])

### Optional blobtools output
blobtools = expand("results/annotation/{assembler}/{source}/{sample_id}/taxonomy/{sample_id}.bestsum.phylum.p5.span.100.exclude_other.blobplot.bam0.png",
            assembler = config["assembler"], sample_id = samples.keys(), source = config["read_source"])
blobtools += expand("results/annotation/{assembler}/{source}/{sample_id}/taxonomy/{sample_id}.bestsum.order.p6.span.100.exclude_other.blobplot.bam0.png",
            assembler = config["assembler"], sample_id = samples.keys(), source = config["read_source"])

### Optional kraken output
kraken_output = expand("results/kraken/{sample_id}.{suffix}", sample_id = samples.keys(), suffix = ["out.gz","kreport"])

## Co-assemblies
co_assembly = expand("results/co-assembly/{assembler}/{assembly}/final.fa", assembly = assemblies.keys(), assembler = config["assembler"])
co_assembly_stats = expand("results/report/co-assembly/{assembler}.{assembly}_assembly_stats.tsv", assembly = assemblies.keys(), assembler = config["assembler"])
co_assembly.append(co_assembly_stats)

## Annotations
eggnog_co = expand("results/collated/co-assembly/{assembler}/{assembly}/eggNOG/{db}.{fc}.tsv",
            db = ["enzymes","pathways","pathways.norm","modules","kos","tc","cazy"],
            fc = ["raw","tpm"],
            assembly = assemblies.keys(), assembler = config["assembler"])
eggnog_co_tax = expand("results/collated/co-assembly/{assembler}/{assembly}/eggNOG_taxonomy/{tax_rank}.{tax_name}.{db}.{fc}.tsv",
            assembly = assemblies.keys(), assembler = config["assembler"],
            fc = ["raw","tpm"],
            db = ["kos","enzymes","modules","pathways","tc","cazy"],
            tax_rank = config["tax_rank"],
            tax_name = config["tax_name"])
abundance_co = expand("results/collated/co-assembly/{assembler}/{assembly}/abundance/{assembly}.{fc}.tsv",
            assembly = assemblies.keys(), assembler = config["assembler"],
            fc = ["raw","tpm"])
abundance_co_tax = expand("results/collated/co-assembly/{assembler}/{assembly}/abundance_taxonomy/{tax_rank}.{tax_name}.{assembly}.{fc}.tsv",
            assembly = assemblies.keys(), assembler = config["assembler"],
            fc = ["raw","tpm"],
            tax_rank = config["tax_rank"],
            tax_name = config["tax_name"])
dbCAN_co = expand("results/collated/co-assembly/{assembler}/{assembly}/dbCAN/dbCAN.{fc}.tsv",
            assembly = assemblies.keys(), assembler = config["assembler"],
            fc = ["raw","tpm"])
dbCAN_co_tax = expand("results/collated/co-assembly/{assembler}/{assembly}/dbCAN_taxonomy/{tax_rank}.{tax_name}.dbCAN.{fc}.tsv",
            assembly = assemblies.keys(), assembler = config["assembler"],
            fc = ["raw","tpm"],
            tax_rank = config["tax_rank"],
            tax_name = config["tax_name"])
taxonomy_co = expand("results/collated/co-assembly/{assembler}/{assembly}/taxonomy/taxonomy.{fc}.tsv",
            assembly = assemblies.keys(), assembler = config["assembler"],
            fc = ["raw","tpm"])
# Map and normalize
map_co = expand("results/report/map/{assembler}.{assembly}_map_report.html", assembly = assemblies.keys(), assembler=config["assembler"])
normalize_co = expand("results/collated/co-assembly/{assembler}/{assembly}/abundance/{assembly}.{fc}.tsv",
                assembly = assemblies.keys(), assembler = config["assembler"], fc = ["tpm","raw"])
### Optional blobtools output
blobtools_co = expand("results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/blobtools/{sample_id}.bestsum.phylum.p5.span.100.exclude_other.blobplot.bam0.png",
                 assembly = assemblies.keys(), sample_id = samples.keys(), assembler = config["assembler"])
blobtools_co += expand("results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/blobtools/{sample_id}.bestsum.phylum.p5.span.100.exclude_other.blobplot.read_cov.bam0.png",
                assembly = assemblies.keys(), sample_id = samples.keys(), assembler = config["assembler"])

inputs = []
inputs += preprocess

if config["co_assembly"]:
    inputs += [eggnog_co, map_co, co_assembly_stats, abundance_co, abundance_co_tax, dbCAN_co, taxonomy_co, dbCAN_co_tax, eggnog_co_tax]
if config["single_assembly"]:
    inputs += [eggnog, dbCAN, taxonomy]

rule all:
    input: inputs

## Preprocess master rule ##
rule preprocess:
    input: preprocess

## Filter master rule ##
rule filter:
    input: filtered_reads

## Assembly master rules ##
rule assemble:
    input: assembly

rule co_assemble:
    input: co_assembly

## Map reads to contigs master rule
rule map:
    input: mapres

rule normalize:
    input: normalize

rule map_co:
    input: map_co

rule normalize_co:
    input: normalize_co

## Annotation master rules
rule annotate:
    input: eggnog, dbCAN, taxonomy

rule annotate_co:
    input: eggnog_co, dbCAN_co, taxonomy_co

rule eggnog:
    input: eggnog

rule eggnog_co:
    input: eggnog_co

rule dbcan:
    input: dbCAN

rule dbcan_co:
    input: dbCAN_co

rule taxonomy:
    input: taxonomy

rule taxonomy_co:
    input: taxonomy_co

rule blobtools:
    input: blobtools

rule blobtools_co:
    input: blobtools_co

## Additional utility rules
rule sourmash:
    input: sourmash

rule sourmash_assembly:
    input: sourmash_assembly

rule kraken:
    input: kraken_output
