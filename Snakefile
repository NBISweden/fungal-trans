import pandas as pd
import glob
import os
import sys
import platform
import yaml
from snakemake.utils import min_version, validate
from snakemake.exceptions import WorkflowError
from source.utils.list_jgi_genomes import parse_genomes

include: "source/utils/common.py"

# Set temporary dir
if not os.getenv("TMPDIR"):
    os.environ["TMPDIR"] = "tmp"
    os.makedirs(os.environ["TMPDIR"],exist_ok=True)

# Set default config and validate against schema
if os.path.exists("config.yml"):
    configfile: "config.yml"

validate(config, "config/config.schema.yaml")


workdir: config["workdir"]

# set jgi account info
with open(config["jgi_account_info"], 'r') as fhin:
    config.update(yaml.safe_load(fhin))

# Parse samples and assemblies
samples, map_dict, assemblies = parse_sample_list(config["sample_file_list"], config)
# Parse JGI genomes from fungi info url, also save to file for quick re-use in future runs
genomes = parse_genomes(config["fungi_info"], f=config["fungi_genomes_file"])

strobealign_chunks = {}

def chunks(lst, n):
    """
    Yield successive n-sized chunks from list
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

for i, chunk in enumerate(chunks(genomes.index.tolist(), 100)):
    strobealign_chunks[f"chunk{i}"] = expand("resources/JGI/genomes/{portal}.transcripts.filt.fna.gz", portal=chunk)

# Reads list of additional genomes to include when building custom mmseqs2 database
extra_genomes = {}
config["mmseqs_db_path"] = os.path.join(config["mmseqs_db_dir"], config["mmseqs_db"])
if config["extra_genomes"] != "":
    extra_genomes = parse_extra_genomes(config["extra_genomes"])
    config["refine_taxonomy"] = True
    config["mmseqs_refine_db"] = f"resources/mmseqs2/combined-fungi-{config['mmseqs_db']}_taxonomy"
else:
    config["mmseqs_refine_db"] = ""
    config["refine_taxonomy"] = False
    
wildcard_constraints:
    sample_id = f"({'|'.join(list(samples.keys()))})",
    assembler = "megahit|trinity|transabyss",
    portals = f"({'|'.join(list(genomes.index.tolist()))})",
    taxname = f"({'|'.join(list(config['taxmap'].keys()))})",
    i = "1|2",
    paired_strategy = "one_mapped|both_mapped",
    kraken_db = config["kraken_db"],
    mmseqs_db = config["mmseqs_db"]

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

def all_input(wildcards):
    """
    Defines all inputs for the workflow
    """
    inputs = []
    
    # Preprocessing
    inputs.extend(
        expand(
            "results/preprocess/sortmerna/{sample_id}/{sample_id}_R{i}.mRNA.fastq.gz", 
            sample_id=samples.keys(), 
            i=[1,2]
        )
    )
    inputs.append(
        "results/report/preprocess/preprocess_report.html"
    )
    # Host reads
    if config["host_filter"]:
        inputs.extend(
            expand(
                "results/filtered/{sample_id}/{sample_id}_R{i}.{h}.fastq.gz",
                sample_id = samples.keys(), 
                i=[1,2],
                h=["host","nohost"]
            )
        )
    if config["fungi_filter"]:
        if config["host_filter"]:
            inputs.extend(
                expand(
                    "results/filtered/{sample_id}/{sample_id}_R{i}.fungi.{paired_strategy}.nohost.fastq.gz",
                    sample_id=samples.keys(),
                    i=["1","2"],
                    paired_strategy=config["paired_strategy"],
                )
            )
        # Filtered reads
        inputs.extend(
            expand(
                "results/filtered/{sample_id}/{sample_id}_R{i}.fungi.{paired_strategy}.fastq.gz",
                sample_id=samples.keys(),
                i=[1,2],
                paired_strategy=config["paired_strategy"],
            )
        )
    
    # Single-assemblies
    if config["single_assembly"]:
        # assembly fasta
        inputs.extend(
            expand(
                "results/assembly/{assembler}/{sample_id}/final.fa",
                assembler=config["assembler"],
                sample_id=samples.keys()
            )
        )
        # stats
        inputs.extend(
            expand(
                "results/report/assembly/{assembler}_stats.tsv",
                assembler=config["assembler"]
            )
        )
        # taxonomy
        inputs.extend(
            expand(
                "results/annotation/{assembler}/{sample_id}/taxonomy/{mmseqs_db}/secondpass.parsed.tsv",
                assembler=config["assembler"],
                sample_id=samples.keys(),
                mmseqs_db=config["mmseqs_db"]
            )
        )
        # taxonomic counts
        inputs.extend(
            expand(
                "results/annotation/{assembler}/{sample_id}/taxonomy/taxonomy.raw.tsv",
                assembler=config["assembler"],
                sample_id=samples.keys(),
            )
        )
        # collated taxonomy counts
        inputs.extend(
            expand(
                "results/collated/{assembler}/taxonomy/taxonomy.raw.tsv",
                assembler=config["assembler"],
            )
        )
        # eggnog
        inputs.extend(
            expand(
                "results/collated/{assembler}/eggNOG/{db}.raw.tsv",
                db=["enzymes","pathways","pathways.norm","modules","kos","tc","cazy"],
                assembler=config["assembler"],
            )
        )
        # count tables
        inputs.extend(
            expand(
                "results/annotation/{assembler}/{sample_id}/featureCounts/fc.raw.tsv",
                assembler=config["assembler"],
                sample_id=samples.keys(),
            )
        )
        # map report
        inputs.extend(
            expand(
                "results/report/map/{assembler}_map_report.html",
                assembler=config["assembler"], 
            )
        )
    
    # Kraken
    inputs.extend(
        expand(
            "results/kraken/{kraken_db}/{sample_id}/{sample_id}.{suffix}",
            kraken_db=config["kraken_db"],
            sample_id=samples.keys(),
            suffix=["out.gz","kreport"]
        )
    )
    # kraken taxbins
    inputs.extend(
        expand(
            "results/kraken/{kraken_db}/{sample_id}/taxbins/{taxname}_R{i}.nohost.fastq.gz",
            kraken_db=config["kraken_db"],
            sample_id=samples.keys(),
            taxname=config["taxmap"].keys(),
            i=[1,2]
        )
    )
    if config["kraken_filter"]:
        inputs.extend(
            expand(
                "results/filtered/{sample_id}/{sample_id}_R{i}.fungi.{paired_strategy}.nohost.kraken.fastq.gz",
                sample_id=samples.keys(),
                i=[1,2],
                paired_strategy=config["paired_strategy"],
            )
        )

    # Co-assemblies
    if config["co_assembly"]:
        # assembly fasta
        inputs.extend(
            expand(
                "results/co-assembly/{assembler}/{assembly}/final.fa", 
                assembly=assemblies.keys(), 
                assembler=config["assembler"]
            )
        )
        # stats
        inputs.extend(
            expand(
                "results/report/co-assembly/{assembler}.{assembly}_assembly_stats.tsv", 
                assembly=assemblies.keys(), 
                assembler=config["assembler"]
            )
        )
        # eggnog
        inputs.extend(
            expand(
                "results/collated/co-assembly/{assembler}/{assembly}/eggNOG/{db}.raw.tsv",
                db=["enzymes","pathways","pathways.norm","modules","kos","tc","cazy"],
                assembly=assemblies.keys(), 
                assembler=config["assembler"]
            )
        )
        # eggnog taxonomy
        inputs.extend(
            expand(
                "results/collated/co-assembly/{assembler}/{assembly}/eggNOG_taxonomy/{tax_rank}.{tax_name}.{db}.raw.tsv",
                assembly=assemblies.keys(), 
                assembler=config["assembler"],
                db=["kos","enzymes","modules","pathways","tc","cazy"],
                tax_rank=config["tax_rank"],
                tax_name=config["tax_name"]
            )
        )
        # count tables
        inputs.extend(
            expand(
                "results/collated/co-assembly/{assembler}/{assembly}/abundance/{assembly}.raw.tsv",
                assembly=assemblies.keys(), 
                assembler=config["assembler"],
            )
        )
        inputs.extend(
            expand(
                "results/collated/co-assembly/{assembler}/{assembly}/abundance/{assembly}.{c}.tsv",
                assembly=assemblies.keys(),
                assembler=config["assembler"],
                c=["est_counts", "tpm"]
            )
        )
        # taxonomy counts
        inputs.extend(
            expand(
                "results/collated/co-assembly/{assembler}/{assembly}/taxonomy/taxonomy.raw.tsv",
                assembly=assemblies.keys(), 
                assembler=config["assembler"],
            )
        )
        # mapping report
        inputs.extend(
            expand(
                "results/report/map/{assembler}_{assembly}_map_report.html", 
                assembly=assemblies.keys(),
                assembler=config["assembler"]
            )
        )
    return inputs



rule all:
    input: all_input
