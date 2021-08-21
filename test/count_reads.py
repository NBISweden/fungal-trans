#!/usr/bin/env python
import os.path
import sys

from Bio.SeqIO import parse
import gzip as gz
import pandas as pd
from argparse import ArgumentParser


def count_host_fungi(f):
    counts = {'host': 0, 'fungi': 0}
    with gz.open(f, 'rt') as fh:
        for record in parse(fh, "fastq"):
            if (record.id).startswith("host"):
                counts["host"] += 1
            elif (record.id).startswith("fungi"):
                counts["fungi"] += 1
    return counts


def main(args):
    aligner = args.host_aligner
    c = {}
    for sample in ["needle", "root"]:
        for pair in [1, 2]:
            for f in [f"results/intermediate/{sample}_R{pair}.fastq.gz",
                      f"results/preprocess/{sample}_R{pair}.cut.fastq.gz",
                      f"results/preprocess/{sample}_R{pair}.cut.trim.fastq.gz",
                      f"results/bowtie2/{sample}/{sample}_R{pair}.fungi.fastq.gz",
                      f"results/bowtie2/{sample}/{sample}_R{pair}.nonfungi.fastq.gz",
                      f"results/{aligner}/{sample}/{sample}_R{pair}.fungi.nohost.fastq.gz",
                      f"results/{aligner}/{sample}/{sample}_R{pair}.fungi.putative-host.fastq.gz",
                      f"results/host/{sample}_R{pair}.host.fastq.gz"]:
                c[os.path.basename(f)] = count_host_fungi(f)

    c_paired = {}
    for key in [k for k in c.keys() if "_R1" in k]:
        key2 = key.replace("_R1", "_R2")
        host = "{}/{}".format(c[key]["host"], c[key2]["host"])
        fung = "{}/{}".format(c[key]["fungi"], c[key2]["fungi"])
        newkey = key.replace("_R1", "_R(1/2)")
        c_paired[newkey] = {'host': host, 'fungi': fung}
    df = pd.DataFrame(c_paired)
    df.T.to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--host_aligner", default="star", type=str,
                        help="Host aligner directory")
    args = parser.parse_args()
    main(args)
