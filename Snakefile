import pandas as pd
import glob
import os
import sys
import platform
from source.utils.parse import parse_sample_list
from snakemake.utils import min_version, validate
from snakemake.exceptions import WorkflowError


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

wildcard_constraints:
    sample_id = "[A-Za-z0-9_\-\.]+",
    assembler = "megahit|trinity|transabyss",
    filter_source = "unfiltered|filtered|taxmapper|bowtie2"

# Set default config and validate against schema
if os.path.exists("config.yml"):
    configfile: "config.yml"
try:
    validate(config, "config/config.schema.yaml")
except WorkflowError as e:
    parse_validation_error(e)
    sys.exit()

workdir: config["workdir"]

# Get environment info
pythonpath = sys.executable
envdir = '/'.join(pythonpath.split("/")[0:-2])
system = platform.system()
config["system"] = system

# Parse samples and assemblies
samples, map_dict, assemblies = parse_sample_list(config["sample_file_list"], config)

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
preprocess = expand("results/preprocess/{sample_id}_R{i}.cut.trim.fastq.gz", sample_id = samples.keys(), i = [1,2])
preprocess += ["results/report/preprocess/preprocess_report.html"]

## Host reads
host_reads = expand("results/host/{sample_id}_R{i}.host.fastq.gz",
            sample_id = samples.keys(), i = [1,2])
## Taxmapper filter
taxmapper_filter = expand("results/taxmapper/{sample_id}/{sample_id}_R1.cut.trim.filtered.fastq.gz", sample_id = samples.keys())
taxmapper_filter += expand("results/report/taxmapper/taxa_freq_norm_lvl{i}.svg", i = [1,2])
## Bowtie filter
bowtie_filter = expand("results/bowtie2/{sample_id}/{sample_id}_R{i}.fungi.nohost.fastq.gz",
            sample_id = samples.keys(), i = [1,2])
bowtie_filter += ["results/report/filtering/filter_report.html"]
bowtie_filter += host_reads
## Union filter
union_filter = expand("results/filtered/{sample_id}/{sample_id}_R{i}.filtered.union.fastq.gz",
            sample_id = samples.keys(), i = [1,2])
union_filter += host_reads

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
mapres += expand("results/report/map/{assembler}/{assembler}_{source}_map_report.html",
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
## Fastuniq
fastuniq = expand("results/fastuniq/{assembly}/R{i}.fastuniq.gz", assembly = assemblies.keys(), i=[1,2])
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
map_co = expand("results/report/map/{assembler}/{assembly}_map_report.html", assembly = assemblies.keys(), assembler=config["assembler"])
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

if config["read_source"] == "bowtie2":
    filter_input = bowtie_filter
elif config["read_source"] == "taxmapper":
    filter_input = taxmapper_filter
elif config["read_source"] == "filtered":
    filter_input = union_filter
else:
    filter_input = preprocess

rule all:
    input: inputs

## Preprocess master rule ##
rule preprocess:
    input: preprocess

## Filter master rules ##
rule taxmapper_filter:
    input: taxmapper_filter

rule bowtie_filter:
    input: bowtie_filter

rule union_filter:
    input: union_filter

rule filter:
    input: filter_input

## Assembly master rules ##
rule assemble:
    input: assembly

rule fastuniq:
    input: fastuniq

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
