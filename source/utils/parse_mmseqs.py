#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser

def split_ranks(df, ranks):
    uniq_lin = df.lineage.unique()
    d = {}
    for lin in uniq_lin:
        d[lin] = dict(zip(ranks,lin.split(";")))
    return pd.DataFrame(d).T

def parse_mmseqs2(tsv, ranks):
    # Read taxonomic assignments and fill NA values
    taxdf = pd.read_csv(tsv, usecols=[0,1,2,3,4], sep="\t", header=None, index_col=0, names=["query","taxid","lca_rank","scientific_name","lineage"])
    taxdf.fillna("unclassified", inplace=True)
    # Extract ranks from unique assignments
    lineage_df = split_ranks(taxdf, ranks)
    # Merge rank assignments to taxonomy dataframe
    dataframe = pd.merge(taxdf, lineage_df, left_on="lineage", right_index=True)
    # Fill NA values
    dataframe.loc[dataframe.lineage != dataframe.lineage, "lineage"] = ";".join(["unclassified"]*len(ranks))
    dataframe.fillna("unclassified", inplace=True)
    return dataframe

def main(args):
    dataframe = parse_mmseqs2(args.input, args.ranks)
    dataframe.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", help="Input MMseqs2 TSV file", required=True)
    parser.add_argument("-o", "--output", help="Output taxonomy TSV file", required=True)
    parser.add_argument("-r", "--ranks", help="Ranks to extract from lineage", nargs="+", default=["superkingdom","phylum","class","order","family","genus","species"])
    args = parser.parse_args()
    main(args)