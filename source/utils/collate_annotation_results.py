#!/usr/bin/env python


from argparse import ArgumentParser
import pandas
from glob import glob
import os


def count_pathways_in_sample(files):
    sample_counts = {}

    cols = ["KEGG Terms"]
    for f in files:
        sample = f.split("/")[-3]
        sample_counts[sample] = {}
        _df = pd.read_table(f, header=0, index_col=0)
        _df_filtered = _df.loc[(_df["Tax Scope"] == "Fungi") & (_df["KEGG Terms"] == _df["KEGG Terms"]), cols]
        total_fungi = len(_df.loc[(_df["Tax Scope"] == "Fungi")])
        sample_counts[sample]["Total_fungi"] = total_fungi
        counts = _df_filtered.reset_index().groupby("KEGG Terms").count()
        for item in counts.index:
            c = counts.loc[item]
            item = item.replace("),", "|").split("|")
            for pathway in item.split("|"):
                pathway = pathway.lstrip()
                if pathway[-1] != ")":
                    pathway += ")"
                try:
                    sample_counts[sample][pathway] += int(c)
                except KeyError:
                    sample_counts[sample][pathway] = int(c)
    return sample_counts


def main(args):
    entap_files = glob(os.path.join(args.dir, "*", "final_results", "final_annotations_lvl0.tsv"))
    sample_counts = count_pathways_in_sample(entap_files)
    sample_counts_df = pd.DataFrame(sample_counts)
    mapids = ["map" + x.split("(")[-1].rstrip(")") for x in counts_df.index]
    sample_counts_df = sample_counts_df.assign(pathway_id = pd.Series(mapids, index=sample_counts_df.index))
    sample_counts_df.to_csv(args.outfile, sep="\t")
    

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-d", "--dir", type=str, default="results/annotation",
                        help="Results dir base path")
    parser.add_argument("-o", "--outfile", type=str, default="results/collated/pathway_counts.tsv",
                        help="Write collated results to this directory")
    args = parser.parse_args()
    main(args)
