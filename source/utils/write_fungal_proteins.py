#!/usr/bin/env python

import polars as pl
from argparse import ArgumentParser
from Bio.SeqIO import parse
from tqdm import tqdm
import sys
from collections import defaultdict


def get_transcript_ids_from_gff(lazy_gff, protein_ids, attribute_type="Parent"):
    # get transcripts from gff
    transcript_ids = (
        (
            # select required columns
            lazy_gff.select(["column_1", "column_3", "column_9"])
            # rename columns
            .rename(
                {
                    "column_1": "transcript_id",
                    "column_3": "type",
                    "column_9": "attributes",
                }
            )
            # filter to only CDS features
            .filter(pl.col("type") == "CDS")
            # split attribute types
            .with_columns(pl.col("attributes").str.split(by=";").alias("attributes"))
            # explode to create rows for each attribute type
            .explode("attributes")
            # filter to only rows with attribute
            .filter(pl.col("attributes").str.starts_with(attribute_type))
            # remove attribute key
            .with_columns(
                query=pl.col("attributes").str.replace(f"{attribute_type}=", "")
            )
            # filter to only show transcripts containing fungal proteins
            .filter(pl.col("query").is_in(protein_ids))
            # select unique transcript_ids
            .select("transcript_id").unique()
        )
        .collect()
        .to_series()
        .to_list()
    )
    return transcript_ids


def write_seqs(infile, outfile, seq_ids):
    d = defaultdict.fromkeys(seq_ids)
    written = 0
    with open(outfile, "w") as fhout:
        for record in tqdm(parse(infile, "fasta"), unit=" seqs", desc="Finding seqs"):
            try:
                d[record.id]
            except KeyError:
                continue
            x = fhout.write(f">{record.description}\n{record.seq}\n")
            written += 1
            if written == len(seq_ids):
                return


def main(args):
    infile = args.infile
    gff_file = args.gff
    bed_file = args.bed
    cds_file = args.cds
    pep_file = args.pep
    attribute_type = args.attribute_type
    outdir = args.outdir
    prefix = args.prefix
    lazy_df = pl.scan_csv(infile, separator="\t")
    # get protein ids from taxonomy file
    sys.stderr.write(f"Loading taxonomy from {infile}\n")
    fungi_proteins = (
        (lazy_df.filter(pl.col("kingdom") == "Fungi").select("query"))
        .collect()
        .to_series()
        .to_list()
    )
    sys.stderr.write(f"Found {len(fungi_proteins)} fungal proteins\n")
    # write taxonomy file
    outfile = f"{outdir}/{prefix}.taxonomy.tsv"
    sys.stderr.write(f"Writing fungal taxonomy to {outfile}\n")
    lazy_df.filter(pl.col("kingdom") == "Fungi").sink_csv(
        f"{outdir}/{prefix}.taxonomy.tsv",
        separator="\t",
        include_header=True,
        quote_style="never",
    )
    if gff_file:
        sys.stderr.write(f"Identifying fungal transcripts using {gff_file}\n")
        lazy_gff = pl.scan_csv(gff_file, separator="\t", has_header=False)
        fungi_transcripts = get_transcript_ids_from_gff(
            lazy_gff, fungi_proteins, attribute_type
        )
        outfile = f"{outdir}/{prefix}.gff3"
        sys.stderr.write(f"Writing fungal gff file to {outfile}\n")
        # write fungal gff
        lazy_gff.filter(pl.col("column_1").is_in(fungi_transcripts)).sink_csv(
            outfile,
            separator="\t",
            include_header=False,
            quote_style="never",
        )
    # write bed file
    if bed_file:
        sys.stderr.write(f"Loading bed file from {bed_file}\n")
        lazy_bed = pl.scan_csv(bed_file, separator="\t", skip_lines=1, has_header=False)
        outfile = f"{outdir}/{prefix}.bed"
        sys.stderr.write(f"Writing fungal bed file to {outfile}\n")
        lazy_bed.filter(pl.col("column_1").is_in(fungi_transcripts)).sink_csv(
            outfile,
            separator="\t",
            include_header=False,
            quote_style="never",
        )
    # write cds file
    if cds_file:
        outfile = f"{outdir}/{prefix}.cds"
        sys.stderr.write(
            f"Identifying fungal sequences in {cds_file} and writing to {outfile}\n"
        )
        write_seqs(infile=cds_file, outfile=outfile, seq_ids=fungi_proteins)
    # write pep file
    if pep_file:
        outfile = f"{outdir}/{prefix}.faa"
        sys.stderr.write(
            f"Identifying fungal sequences in {pep_file} and writing to {outfile}\n"
        )
        write_seqs(infile=pep_file, outfile=outfile, seq_ids=fungi_proteins)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i", "--infile", type=str, required=True, help="Taxonomy results for proteins"
    )
    parser.add_argument("-g", "--gff", type=str, help="GFF3 file from genecaller")
    parser.add_argument("-b", "--bed", type=str, help="BED file from genecaller")
    parser.add_argument("-c", "--cds", type=str, help="CDS file from genecaller")
    parser.add_argument("-p", "--pep", type=str, help="PEP file from genecaller")
    parser.add_argument(
        "-o", "--outdir", type=str, required=True, help="Output directory"
    )
    parser.add_argument(
        "-a",
        "--attribute_type",
        type=str,
        default="Parent",
        help="Attribute type in gff file used to link protein ids",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="fungal",
        help="Prefix for files written to output directory",
    )
    args = parser.parse_args()
    main(args)
