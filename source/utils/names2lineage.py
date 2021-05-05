#!/usr/bin/env python

from ete3 import NCBITaxa
from argparse import ArgumentParser
import os
import pandas as pd
import sys
import re
v_re = re.compile("[0-9].[0-9]")

def main():
    parser = ArgumentParser()
    parser.add_argument("namesmap", help="File with ids->names mappings")
    parser.add_argument("-t", "--taxdb", help="Taxonomy database formatted using ete3")
    args = parser.parse_args()

    if not os.path.exists(args.taxdb):
        with open(args.taxdb, 'w') as fh:
            pass
    ncbi_taxa = NCBITaxa(args.taxdb)

    # Read names
    lineages = {}
    ids = {}
    ranks = ["superkingdom","phylum","class","order","family","genus","species"]   
    with open(args.namesmap, 'r') as fhin:
        for line in fhin:
            line = line.rstrip()
            genome_id, name = line.split("\t")
            ids[genome_id] = name
    for genome_id, name in ids.items():
        lineages[genome_id] = {}
        for rank in ranks:
            lineages[genome_id][rank] = ""
        items = name.split(" ")
        taxid = False
        for i in range(2,len(items)+1):
            query_name = " ".join(items[0:i])
            name2taxid = ncbi_taxa.get_name_translator([query_name])
            try:
                taxid = name2taxid[query_name][0]
                break
            except KeyError:
                continue
        if not taxid:
            sys.stderr.write("No information found for:{}\n".format(genome_id))
            continue
        lineage = ncbi_taxa.get_lineage(taxid)
        rankdict = ncbi_taxa.get_rank(lineage)
        for taxid,rank in rankdict.items():
            if rank in ranks:
                rank_name = ncbi_taxa.get_taxid_translator([taxid])[taxid]
                lineages[genome_id][rank] = rank_name
    df = pd.DataFrame(lineages).T
    df = df.loc[:,ranks]
    df.index.name="genome_id"
    df.to_csv(sys.stdout, sep="\t")


if __name__ == '__main__':
    main()
