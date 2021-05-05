#!/usr/bin/env python


from Bio.SeqIO import parse
from argparse import ArgumentParser


def parse_fasta(args):
    seqs = {}
    gff = []
    for i, record in enumerate(parse(args.infile, "fasta"), start=1):
        contig = record.id
        gene_id = "{sample}.gene{i} {contig}".format(sample=args.sample, i=i, contig=contig)
        strand, start, stop = (record.description).split("|")[-3:]
        seqs[gene_id] = str(record.seq)
        gff.append([contig, "GeneMark", "CDS", start, stop, ".",
            strand, "0", "gene_id={sample}.gene{i}".format(sample=args.sample, i=i)])
    return seqs, gff


def write_fasta(s, outfasta):
    with open(outfasta, 'w') as fhout:
        for seqid, seq in s.items():
            fhout.write(">{seqid}\n{seq}\n".format(seqid=seqid, seq=seq))


def write_gff(g, outgff):
    with open(outgff, 'w') as fhout:
        for line in g:
            fhout.write("{}\n".format("\t".join(line)))


def main():
    parser = ArgumentParser()
    parser.add_argument("infile", type=str,
                        help="Input protein fasta file")
    parser.add_argument("outfasta", type=str,
                        help="Output reformatted fasta file")
    parser.add_argument("outgff", type=str,
                        help="Output GFF format file")
    parser.add_argument("sample", type=str,
                        help="Sample id to add to headers")
    args = parser.parse_args()
    seqs, gff = parse_fasta(args)
    write_gff(gff, args.outgff)
    write_fasta(seqs, args.outfasta)


if __name__ == '__main__':
    main()
