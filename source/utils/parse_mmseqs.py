#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser

def split_ranks(df, ranks):
    uniq_lin = df.lineage.unique()
    d = {}
    for lin in uniq_lin:
        d[lin] = dict(zip(ranks,lin.split(";")))
    return pd.DataFrame(d).T

def add_last_known(df, ranks):
    d = {}
    df.fillna("unclassified", inplace=True)
    for row in df.iterrows():
        lineage, items = row
        d[lineage] = {}
        last_known = ""
        for i, rank in enumerate(ranks):
            if not items.iloc[i].startswith("uc_") and not items.iloc[i] in ["unclassified", "unknown"]:
                last_known = items.iloc[i]
                d[lineage][rank] = last_known
            else:
                if last_known != "":
                    d[lineage][rank] = f"unclassified.{last_known}"
                else:
                    d[lineage][rank] = "unclassified"
    return pd.DataFrame(d).T.loc[:, ranks]


def parse_mmseqs2(tsv, ranks):
    # Read taxonomic assignments and fill NA values
    taxdf = pd.read_csv(tsv, usecols=[0,1,2,3,4], sep="\t", header=None, index_col=0, names=["query","taxid","lca_rank","scientific_name","lineage"])
    taxdf.fillna("unclassified", inplace=True)
    # Extract ranks from unique assignments
    lineage_df = split_ranks(taxdf, ranks)
    # Set Eukaryota as first rank for Fungi
    lineage_df.loc[lineage_df.kingdom=="Fungi", ranks[0]] = "Eukaryota"
    # Add last known rank to unclassified assignments
    lineage_df = add_last_known(lineage_df, ranks)
    # Merge rank assignments to taxonomy dataframe
    dataframe = pd.merge(taxdf, lineage_df, left_on="lineage", right_index=True)
    return dataframe.loc[:, ranks]

def main(args):
    dataframe = parse_mmseqs2(args.input, args.ranks)
    dataframe.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", help="Input MMseqs2 TSV file", required=True)
    parser.add_argument("-o", "--output", help="Output taxonomy TSV file", required=True)
    parser.add_argument("-r", "--ranks", help="Ranks to extract from lineage", nargs="+", default=["superkingdom","kingdom","phylum","class","order","family","genus","species"])
    args = parser.parse_args()
    main(args)