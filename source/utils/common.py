#!/usr/bin/env

import pandas as pd
import os
import sys


def parse_sample_list(f, config):
    """
    Parse the sample list and set up input file names for the assembly
    based on configuration

    Parameters
    ----------
    f: str
        Path to the sample list
    config: dict
        Configuration dictionary

    Returns
    -------
    dict
        Dictionary with sample names as keys and read files as values
    dict
        Dictionary with sample names as keys and read files as values
    dict
        Dictionary with assembly names as keys and read files as values
    """
    samples = {}
    map_dict = {}
    kraken_db = config["kraken_db"]
    paired_strategy = config["paired_strategy"]
    k = config["strobealign_strobe_len"]
    kraken_db = config["kraken_db"]
    sample_string = ""
    if config["fungi_filter"]:
        sample_string = f"fungi.{paired_strategy}."
    if config["host_filter"]:
        sample_string += "nohost."
    if config["kraken_filter"]:
        sample_string += "kraken."

    dir = config["datadir"]
    # Read sample list
    df = pd.read_csv(f, comment="#", header=0, sep="\t", index_col=0, dtype=str)
    df.fillna("", inplace=True)
    # If it's just a one-column file, expand it by assuming that the
    # sra accessions are in the first column
    if df.shape[1] == 0 or "accession" not in df.columns:
        df = df.assign(accession=pd.Series(df.index, index=df.index))
        df.index.name = "Sample"
    assemblies = {}
    if "assembly" in df.columns:
        for assembly in df.assembly.unique():
            assemblies[assembly] = {"R1": [], "R2": []}
    for sample in df.index:
        accession = ""
        if "Read_file" and "Pair_file" in df.columns:
            R1 = df.loc[sample, "Read_file"]
            R2 = df.loc[sample, "Pair_file"]
        elif "accession" in df.columns:
            accession = df.loc[sample, "accession"]
            R1 = f"{dir}/{sample}_1.fastq.gz"
            R2 = f"{dir}/{sample}_2.fastq.gz"
        # Define reads for assembly
        if sample_string == "":
            R1_f = f"results/preprocess/sortmerna/{sample}/{sample}_R1.mRNA.fastq.gz"
            R2_f = f"results/preprocess/sortmerna/{sample}/{sample}_R2.mRNA.fastq.gz"
        elif sample_string == "nohost.":
            R1_f = f"results/filtered/{sample}/{sample}_R1.nohost.fastq.gz"
            R2_f = f"results/filtered/{sample}/{sample}_R2.nohost.fastq.gz"
        elif sample_string == "kraken.":
            R1_f = f"results/kraken/{kraken_db}/{sample}/taxbins/Fungi_R1.fastq.gz"
            R2_f = f"results/kraken/{kraken_db}/{sample}/taxbins/Fungi_R2.fastq.gz"
        else:
            R1_f = f"results/filtered/{sample}/{sample}_R1.{sample_string}fastq.gz"
            R2_f = f"results/filtered/{sample}/{sample}_R2.{sample_string}fastq.gz"
        map_dict[sample] = {"R1": R1_f, "R2": R2_f}
        if len(assemblies.keys()) > 0:
            if df.loc[sample, "assembly"] != "":
                assembly = df.loc[sample, "assembly"]
                assemblies[assembly]["R1"].append(R1_f)
                assemblies[assembly]["R2"].append(R2_f)

        samples[sample] = {"R1": R1, "R2": R2, "accession": accession}
    return samples, map_dict, assemblies


def fungi_input(wildcards):
    """
    Get the input files for the fungi pipeline

    Parameters
    ----------
    wildcards: dict
        Wildcards from snakemake

    Returns
    -------
    list
        List of input files
    """
    paired_strategy = config["paired_strategy"]
    kraken_db = config["kraken_db"]
    k = config["strobealign_strobe_len"]
    d = "unfiltered"
    if config["filter_reads"]:
        R1 = f"results/filtered/{wildcards.sample_id}/{paired_strategy}/k{k}/{kraken_db}/{wildcards.sample_id}_R1.fungi.nohost.kraken.fastq.gz"
        R2 = f"results/filtered/{wildcards.sample_id}/{paired_strategy}/k{k}/{kraken_db}/{wildcards.sample_id}_R2.fungi.nohost.kraken.fastq.gz"
    else:
        R1 = f"results/unfiltered/{wildcards.sample_id}/{wildcards.sample_id}_R1.fastq.gz"
        R2 = f"results/unfiltered/{wildcards.sample_id}/{wildcards.sample_id}_R2.fastq.gz"
    return [R1, R2]


def parse_extra_genomes(f):
    """
    Parse the extra genomes file

    Parameters
    ----------
    f: str
        Path to the extra genomes file

    Returns
    -------
    dict
        Dictionary with portal names as keys and taxids as values
    """
    extra_genomes = {}
    with open(f, "r") as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            line = line.strip()
            portal, taxid = line.rsplit()[0:2]
            try:
                taxid = int(taxid)
            except ValueError:
                continue
            extra_genomes[portal] = taxid
    return extra_genomes


def get_mmseq_taxdb(wildcards):
    if config["refine_taxonomy"]:
        return config["mmseqs_refine_db"]
    return os.path.join(config["mmseqs_db_dir"], config["mmseqs_db"], "_taxonomy")


def write_fungal_proteins(parsed, indir, outdir):
    """
    Write fungal coding sequences and proteins to separate files

    Parameters
    ----------
    parsed: str
        Path to the parsed taxonomy file
    indir: str
        Path to the input directory
    outdir: str
        Path to the output directory
    """
    import pandas as pd
    import os

    df = pd.read_csv(parsed, sep="\t", header=0, index_col=0)
    fungi = df.loc[df["kingdom"] == "Fungi"].index.tolist()
    gff = os.path.join(indir, "final.fa.transdecoder.gff3")
    with open(f"{outdir}/fungal.gff3", "w") as fhout, open(gff, "r") as fhin:
        for line in fhin:
            if line.startswith("#"):
                fhout.write(line)
            else:
                if line.split("\t")[-1].split(";")[0].replace("ID=", "") in fungi:
                    fhout.write(line)
    bed = os.path.join(indir, "final.fa.transdecoder.bed")
    with open(f"{outdir}/fungal.bed", "w") as fhout, open(bed, "r") as fhin:
        for i, line in enumerate(fhin):
            if i == 0:
                fhout.write(line)
            else:
                if line.split("\t")[3].split(";")[0].replace("ID=", "") in fungi:
                    fhout.write(line)
    with open("fungi.ids", "w") as fhout:
        for i in fungi:
            fhout.write(i + "\n")
    shell(
        "seqkit grep -f fungi.ids {indir}/final.fa.transdecoder.cds > {outdir}/fungal.cds"
    )
    shell(
        "seqkit grep -f fungi.ids {indir}/final.fa.transdecoder.pep > {outdir}/fungal.faa"
    )
    df.loc[fungi].to_csv(f"{outdir}/fungal.taxonomy.tsv", sep="\t", index=True)


def slurm_mem_partition(mb):
    """
    Determine the SLURM memory partition based on the memory

    mb: int
        Memory in MB of the input

    Returns
    -------
    str
        SLURM memory partition
    """
    if mb > 900000:
        return config["slurm_partitions"]["large_memory"]
    elif mb > 220000:
        return config["slurm_partitions"]["medium_memory"]
    return config["slurm_partitions"]["small_memory"]
