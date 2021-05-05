#!/usr/bin/env python

# Will create a JGI fasta file with headers like this:
# jgi|Aaoar1|403682|fgenesh1_pg.1_#_134
# into a GTF format file for counting reads with featureCounts:
# jgi|Aaoar1|403682|fgenesh1_pg.1_#_134 JGI CDS 1 630 . + 0 genome_id=Aaoar1

from Bio.SeqIO import parse
from argparse import ArgumentParser, FileType
import sys


def read_idfile(f):
    ids = {}
    with open(f, 'r') as fhin:
        for line in fhin:
            ids[line.rstrip()] = ""
    return ids

def main():
    parser = ArgumentParser()
    parser.add_argument("fastafile", nargs='?', type=FileType('r'),
            help="Fasta file downloaded from JGI", default=sys.stdin)
    parser.add_argument("--idfile", help="Limit records to ids found in this file")
    parser.add_argument("--genome_id", help="Specify optional genome ID to use for all records")
    args = parser.parse_args()
    idfile = False
    idcount = 0
    if args.idfile:
        ids = read_idfile(args.idfile)
        idcount = len(ids)
        if idcount > 0:
            idfile = True
    found_ids = 0
    for record in parse(args.fastafile, 'fasta'):
        if idfile:
            if idcount==found_ids:
                break
            try:
                ids[record.id]
            except KeyError:
                continue
        end = len(record.seq)
        if len(record.id)>99:
            seqid = (record.id)[0:99]
        else:
            seqid = record.id
        if not args.genome_id:
            genome = (record.id).split("|")[1]
        else:
            genome = args.genome_id
        line = "{id}\tJGI\tCDS\t1\t{end}\t.\t+\t0\tgenome_id={genome}\n".format(id=seqid, end=end, genome=genome)
        sys.stdout.write(line)
        found_ids+=1



if __name__ == '__main__':
    main()
