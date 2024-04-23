#!/usr/bin/env python

import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
import os
import sys
from argparse import ArgumentParser

def parse_genomes(url, f = None):
    """
    Parses available genome portals at JGI (default: https://mycocosm.jgi.doe.gov/fungi/fungi.info.html)

    :param url: str
    :return: dict
    """
    def extract_name_portal(genome_table):
        df = {}
        for _, items in genome_table.iterrows():
            name, portal = items["Name"]
            if name == "Name":
                continue
            asm_length, _ = items["Assembly Length"]
            genes, _ = items["# Genes"]
            if asm_length == "NO DATA":
                asm_length = np.nan
            else:
                asm_length = int(asm_length.replace(",",""))
            if genes == "NO DATA":
                genes = np.nan
            else:
                genes = int(genes.replace(",",""))
            df[portal.lstrip("/")] = {'Name': name, 'bp': asm_length, 'genes': genes}
        df = pd.DataFrame(df).T
        df.index.name = "portal"
        return df
    if f is not None and os.path.exists(f):
        genomes = pd.read_csv(f, sep="\t", header=0, index_col=0)
        return genomes
    genome_table = pd.read_html(url, extract_links="body", header=None)[0]
    genome_table.columns = ["##","Name","Assembly Length","# Genes", "Published"]
    genome_table.set_index("##", inplace=True)
    genomes = extract_name_portal(genome_table)
    if f is not None:
        os.makedirs(os.path.dirname(f), exist_ok=True)
        genomes.to_csv(f, sep="\t")
    else:
        with sys.stdout as f:
            genomes.to_csv(f, sep="\t")
    return genomes


def main(args):
    """
    Main function
    """
    genomes = parse_genomes(args.info, args.file)
    

if __name__ == '__main__':
    parser = ArgumentParser(description='Download a table of all available genomes from JGI Mycocosm')
    parser.add_argument('-i', '--info', help='Info file at JGI Mycocosm (default: https://mycocosm.jgi.doe.gov/fungi/fungi.info.html)', default="https://mycocosm.jgi.doe.gov/fungi/fungi.info.html")
    parser.add_argument('-f', '--file', help='Output file name. Will print to stdout if not provided.')
    args = parser.parse_args()
    main(args)