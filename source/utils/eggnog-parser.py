#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import polars as pl
import os
import logging

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


# Download functions
####################

def get_kegg_module_hierarchy(s):
    hier = {}
    # First level is 'ko00002'
    for d1 in s['children']:
        c1 = d1['name']
        for d2 in d1['children']:
            c2 = d2['name']
            for d3 in d2['children']:
                c3 = d3['name']
                for module in d3['children']:
                    module_name = module['name']
                    #kos = module['children']
                    hier[module_name] = {"Module_category1": c1, "Module_category2": c2,
                                         "Module_category3": c3}
                    #for ko in kos:
                    #    ko_name = ko['name']
                    #    ko = ko_name.split(" ")[0]
                    #    hier[module_name]["KOs"].append(ko)
    return hier


def get_kegg_ortholog_hierarchy(s):
    hier = {}
    # First level is 'ko00001'
    for d1 in s['children']:
        c1 = d1['name']
        for d2 in d1['children']:
            c2 = d2['name']
            for d3 in d2['children']:
                c3 = d3['name']
                if not "children" in d3.keys():
                    continue
                for ko in d3['children']:
                    ko_name = ko['name'].split("\t")[0]
                    ko_id = ko_name.split(" ")[0]
                    if "[EC:" in ko_name:
                        enzymes = ko_name.split("[")[-1].split("]")[0].lstrip("EC:").split(" ")
                    else:
                        enzymes = []
                    d = {"KO_category1": c1, "KO_category2": c2, "pathway": c3, "name": ko_name, "enzymes": enzymes}
                    try:
                        hier[ko_id].append(d)
                    except KeyError:
                        hier[ko_id] = [d]
    return hier


def get_kegg_ortholog_info(outdir):
    import urllib.request
    import json
    outdir = outdir.rstrip("/")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    url = "https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json"
    logging.info("Fetching ko00001.keg from www.kegg.jp")
    # Download file
    tmp_out = f"{outdir}/ko00001.json"
    urllib.request.urlretrieve(url, tmp_out)
    with open(tmp_out) as fh:
        s = json.load(fh)
    hier = get_kegg_ortholog_hierarchy(s)
    pathways = {}
    ec2path = {}
    ko_out = f"{outdir}/kegg_kos.tsv"
    ko2path_out = f"{outdir}/kegg_ko2pathways.tsv"
    ko2ec_out = f"{outdir}/kegg_ko2ec.tsv"
    ec2path_out = f"{outdir}/kegg_ec2pathways.tsv"
    pathways_out = f"{outdir}/kegg_pathways.tsv"
    for f in [ko_out, ko2path_out, ko2ec_out, ec2path_out, pathways_out]:
        logging.info(f"Writing to {f}")
    # Write KEGG Ortholog names, KEGG Ortholog -> Pathway map, and KEGG Ortholog -> Enzyme map
    with open(ko_out, 'w') as fh_kos, open(ko2path_out, 'w') as fh_ko2path, open(ko2ec_out, 'w') as fh_ko2ec:
        fh_kos.write("ko\tKO_name\n")
        for ko_id, l in hier.items():
            for i, d in enumerate(l):
                if i == 0:
                    fh_kos.write(f"{ko_id}\t{d['name']}\n")
                    for enzyme in d["enzymes"]:
                        fh_ko2ec.write(f"{ko_id}\t{enzyme}\n")
                fh_ko2path.write(f"{ko_id}\tmap{d['pathway'].split(' ')[0]}\n")
                pathways[d["pathway"]] = {"Pathway_category1": d["KO_category1"],
                                          "Pathway_category2": d["KO_category2"]}
                for enzyme in d["enzymes"]:
                    try:
                        ec2path[enzyme].append(f"map{d['pathway'].split(' ')[0]}")
                    except KeyError:
                        ec2path[enzyme] = [f"map{d['pathway'].split(' ')[0]}"]
    # Write Pathway information
    with open(pathways_out, 'w') as fh_pathways:
        fh_pathways.write("Pathway_id\tPathway_name\tPathway_category1\tPathway_category2\n")
        for pathway, d in pathways.items():
            pathway_id = pathway.split(" ")[0]
            pathway_name = pathway.replace(f"{pathway_id} ", "")
            fh_pathways.write(f"map{pathway_id}\t{pathway_name}\t{d['Pathway_category1']}\t{d['Pathway_category2']}\n")
    # Write Enzyme -> Pathway map
    with open(ec2path_out, 'w') as fh_ec2path:
        for enzyme, l in ec2path.items():
            for pathway in set(l):
                fh_ec2path.write(f"{enzyme}\t{pathway}\n")


def get_kegg_module_info(outdir):
    import urllib.request
    import json
    outdir = outdir.rstrip("/")
    # Process KEGG Module information
    logging.info("Fetching ko00002.keg from www.kegg.jp")
    url = "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00002.keg&format=json"
    tmp_out = f"{outdir}/ko00002.json"
    urllib.request.urlretrieve(url, tmp_out)
    modules_out = f"{outdir}/kegg_modules.tsv"
    with open(tmp_out) as fh:
        s = json.load(fh)
    hier = get_kegg_module_hierarchy(s)
    logging.info(f"Writing to {modules_out}")
    with open(modules_out, 'w') as fh_modules:
        fh_modules.write("Module_id\tModule_name\tModule_category1\tModule_category2\tModule_category3\n")
        for module_name, d in hier.items():
            module_key = module_name.split(" ")[0]
            module_name = " ".join(module_name.split(" ")[1:]).lstrip()
            fh_modules.write(f"{module_key}\t{module_name}\t{d['Module_category1']}\t{d['Module_category2']}\t{d['Module_category3']}\n")


# Parse functions
#################
def parse_ko_annotations(annotations, info_file, outfile, col):
    logging.info(f"Reading {col} annotations from {annotations}")
    df = pl.scan_csv(annotations, separator="\t", comment_prefix="##").select(
        ["#query",col]
        ).filter(
            pl.col(col)!="-"
        ).collect(engine="streaming")
    df = df.rename({'#query': 'orf'})
    df = df.with_columns(
        explode=pl.col(col).str.split(",")
        ).explode("explode").drop(col).rename({'explode': col})
    if col=="KEGG_ko":
        df = df.with_columns(pl.col(col).str.replace("ko:", ""))
    elif col=="KEGG_Pathway":
        df = df.filter(pl.col(col).str.starts_with("map"))
    if info_file:
        info_df = pl.read_csv(info_file, separator="\t", ignore_errors=True)
        df = df.join(info_df, right_on=info_df.columns[0], left_on=col, how="left")
        df = df.fill_null("Unknown")
    df = df.sort("orf")
    logging.info(f"Writing annotations to {outfile}")
    df.write_csv(outfile, separator="\t")

# Quantify functions
####################

def process_and_sum(q_df, annot_df):
    # Merge annotations and abundance, keep ORFs without annotation as "Unclassified"
    annot_q_df = pd.merge(annot_df, q_df, left_index=True, right_index=True, how="right")
    annot_q_df.fillna("Unclassified", inplace=True)
    feature_cols = annot_df.columns
    annot_q_sum = annot_q_df.groupby(list(feature_cols)).sum().reset_index()
    annot_q_sum.set_index(feature_cols[0], inplace=True)
    return annot_q_sum


def sum_to_features(abundance, parsed):
    parsed_df = pd.read_csv(parsed, index_col=0, sep="\t")
    abundance_df = pd.read_csv(abundance, index_col=0, sep="\t")
    feature_sum = process_and_sum(abundance_df, parsed_df)
    return feature_sum


def normalize(q_df, parsed, normalize_file):
    info_df = pd.read_csv(normalize_file, header=None, sep="\t")
    info_norm_df = info_df.groupby(1).count()
    info_norm_df.columns = ["norm_factor"]
    annot_df = pd.read_csv(parsed, index_col=0, sep="\t")
    annot_cols = list(set(annot_df.columns).intersection(set(q_df.columns)))
    sample_cols = list(set(q_df.columns).difference(annot_cols))
    q_df = pd.merge(q_df, info_norm_df, left_index=True, right_index=True, how="left")
    q_df.fillna(1, inplace=True)
    q_df.loc[:, sample_cols] = q_df[sample_cols].div(q_df["norm_factor"], axis=0)
    q_df.drop("norm_factor", axis=1, inplace=True)
    return q_df


def merge_files(files, sum_abundance=False):
    df = pd.DataFrame()
    for i, f in enumerate(files):
        _df = pd.read_csv(f, header=0, sep="\t")
        if sum_abundance:
            group_cols = list(_df.columns)[1:-1]
            if len(group_cols) == 0:
                group_cols = [list(_df.columns)[0]]
        else:
            group_cols = list(_df.columns)[0:-1]
        # Sum dataframe
        _df = _df.groupby(group_cols).sum().reset_index()
        if i == 0:
            df = _df.copy()
        else:
            df = pd.merge(df, _df, on=group_cols, how="outer")
    df.set_index(group_cols[0], inplace=True)
    df.fillna(0, inplace=True)
    return df


# Parser default functions
##########################


def download(args):
    if not os.path.exists(args.dldir):
        os.makedirs(args.dldir)
    logging.info(f"Downloading KEGG module info to {args.dldir}")
    get_kegg_module_info(args.dldir)
    logging.info(f"Downloading KEGG ortholog info to {args.dldir}")
    get_kegg_ortholog_info(args.dldir)


def parse(args):
    parse_ko_annotations(args.annotations, args.info_file, args.outfile, args.col)


def quant(args):
    feature_sum = sum_to_features(args.abundance, args.parsed)
    if args.normalize:
        if not os.path.exists(args.normalize):
            dldir = os.path.dirname(args.normalize)
            get_kegg_module_info(dldir)
            get_kegg_ortholog_info(dldir)
        feature_sum = normalize(feature_sum, args.parsed, args.normalize)
    feature_sum.to_csv(args.outfile, sep="\t")


def merge(args):
    df = merge_files(files=args.files, sum_abundance=args.sum)
    df.to_csv(args.outfile, sep="\t", index=True, header=True)


def main():
    parser = ArgumentParser()
    subparser = parser.add_subparsers(title="Subcommands")
    # Download parser
    download_parser = subparser.add_parser("download", help="Download KEGG info files")
    download_parser.add_argument("dldir",
                                 help="Write files to this directory. Will be created if missing.")
    download_parser.set_defaults(func=download)
    # Annotation parser
    annot_parser = subparser.add_parser("parse", help="Parse emapper annotation files")
    annot_parser.add_argument("annotations",
                              help="emapper.py annotation files (Typically *.emapper.annotations)")
    annot_parser.add_argument("outfile",
                              help="Output file for parsed results")
    annot_parser.add_argument("--info_file",
                              help="Extra info file to merge annotations with")
    annot_parser.add_argument("--col", type=str, help="Column to parse", default="KEGG_ko")
    annot_parser.set_defaults(func=parse)
    # Sum annotations parser
    quant_parser = subparser.add_parser("quantify", help="Add abundances and summarize KEGG features")
    quant_parser.add_argument("abundance",
                              help="Abundance table")
    quant_parser.add_argument("parsed",
                              help="Parsed file (from parse command)")
    quant_parser.add_argument("outfile",
                              help="Output file with summed abundances for features")
    quant_parser.add_argument("--normalize",
                              help="Specify either kegg_ko2pathways.tsv or kegg_ko2modules.tsv to normalize abundances"
                                   "by size of the pathway/module.")
    quant_parser.set_defaults(func=quant)
    # Merge files parser
    merge_parser = subparser.add_parser("merge", help="Merge two or more files from quantify step")
    merge_parser.add_argument("files", nargs='+',
                              help="Files from quantify step")
    merge_parser.add_argument("outfile", type=str,
                              help="Output file")
    merge_parser.add_argument("--sum", action="store_true",
                              help="When merging, sum abundances to second column. (default = False)")
    merge_parser.set_defaults(func=merge)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
