#!/usr/bin/env python

import pandas as pd
import sys
from argparse import ArgumentParser


def read_df(f):
    df = pd.read_csv(f, header=0, sep="\t", index_col=0)
    return df


def calc_means(dfsum, sample_info, info_col):
    # Create dictionary to rename columns to their corresponding samplegroup
    rename_dict = sample_info.loc[:,info_col].to_dict()
    sample_means = pd.DataFrame()
    for samplegroup in set(rename_dict.values()):
        # Rename columns on the fly then calculate mean across axis 1
        _ = pd.DataFrame(dfsum.rename(columns=rename_dict).loc[:, samplegroup].mean(axis=1))
        _.columns=[samplegroup]
        sample_means = pd.concat([sample_means,_], axis=1)
    return sample_means


def main(args):
    sample_info = read_df(args.sample_info)
    annot = read_df(args.annot)
    annot_col = annot.columns[0]
    tax = read_df(args.tax)
    abund = read_df(args.abund)
    # Merge annotation with taxonomy
    df = pd.merge(annot, tax, left_index=True, right_index=True, how="left")
    # Merge with abundance
    df = pd.merge(df, abund, left_index=True, right_index=True)
    # Sum to annotation and taxonomy
    dfsum = df.groupby([args.rank,annot_col]).sum()
    # Calculate mean across sample groups
    sample_means = calc_means(dfsum, sample_info, args.info_col)
    # Write to file
    sample_means.to_csv(sys.stdout, sep="\t")


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("sample_info", type=str,
                        help="Sample info file")
    parser.add_argument("info_col", type=str,
                        help="Column to groupby in sample info file")
    parser.add_argument("annot", type=str,
                        help="Annotation file (e.g. kos.parsed.tsv")
    parser.add_argument("tax", type=str,
                        help="Taxonomic classification file")
    parser.add_argument("abund", type=str,
                        help="Abundance file")
    parser.add_argument("-r", "--rank", default="kingdom",
                        help="Taxonomic rank to groupby")
    args = parser.parse_args()
    main(args)
