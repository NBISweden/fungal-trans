#!/usr/bin/env python

from argparse import ArgumentParser, FileType
import sys
from Bio.SeqIO import parse
import pandas as pd

def main():
    parser = ArgumentParser()
    parser.add_argument("fastafile", nargs='?', type=FileType('r'), default=sys.stdin)
    parser.add_argument("taxidmap", type=str)
    parser.add_argument("-g", "--genome_id", type=str)
    parser.add_argument("--accessionlength", type=int, default=1000000, help="Maximum allowed accessionlength")
    args = parser.parse_args()
    # Read the taxidmap
    taxidmap = pd.read_table(args.taxidmap, header=0, index_col=0, sep="\t")
    if args.genome_id:
        taxid = taxidmap.loc[args.genome_id,"taxid"]

    for i, record in enumerate(parse(args.fastafile, "fasta"), start=1):
        if args.genome_id:
            genome_id = args.genome_id[0:args.accessionlength-6]
        else:
            genome_id = (record.id).split("_")[0]
            genome_id = genome_id[0:args.accessionlength-6]
            try:
                taxid = taxidmap.loc[genome_id,"taxid"]
            except KeyError:
                continue
        seqid = "{genome_id}.{i}".format(genome_id=genome_id, i=i, record_id=record.id)
        sys.stdout.write(">{seqid}\n{seq}\n".format(seqid=seqid, seq=str(record.seq)))
        sys.stderr.write("{seqid}\t{record_id}\t{taxid}\n".format(seqid=seqid, record_id=record.id, taxid=taxid))


if __name__ == '__main__':
    main()
