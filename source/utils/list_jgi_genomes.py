#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
from argparse import ArgumentParser


def parse_genomes(url, f=None):
    """
    Parses available genome portals at JGI (default: https://mycocosm.jgi.doe.gov/fungi/fungi.info.html)

    If f is given and already exists, it is read from as a cache and url is ignored.

    :param url: str
    :return: dict
    """

    def extract_name_portal(genome_table):
        records = {}
        for _, items in genome_table.iterrows():
            name, portal = items["Name"]
            portal = os.path.basename(portal)
            asm_length, _ = items["Assembly Length"]
            genes, _ = items["# Genes"]
            if asm_length == "NO DATA":
                asm_length = np.nan
            else:
                asm_length = int(asm_length.replace(",", ""))
            if genes == "NO DATA":
                genes = np.nan
            else:
                genes = int(genes.replace(",", ""))
            records[portal.lstrip("/")] = {"Name": name, "bp": asm_length, "genes": genes}
        records = pd.DataFrame(records).T
        records.index.name = "portal"
        return records

    if f is not None and os.path.exists(f):
        genomes = pd.read_csv(f, sep="\t", header=0, index_col=0)
        return genomes
    genome_table = pd.read_html(url, extract_links="body", header=None)[0]
    genome_table.columns = ["##", "Name", "Assembly Length", "# Genes", "Published"]
    genome_table.set_index("##", inplace=True)
    genomes = extract_name_portal(genome_table)
    if f is not None:
        dirname = os.path.dirname(f)
        if dirname:
            os.makedirs(dirname, exist_ok=True)
        genomes.to_csv(f, sep="\t")
    return genomes


def main(args):
    """
    Main function
    """
    genomes = parse_genomes(args.info, args.file)
    if args.file is None:
        print(genomes.to_csv(sep="\t"))


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Download a table of all available genomes from JGI Mycocosm"
    )
    parser.add_argument(
        "-i",
        "--info",
        help="Info file at JGI Mycocosm (default: https://mycocosm.jgi.doe.gov/fungi/fungi.info.html)",
        default="https://mycocosm.jgi.doe.gov/fungi/fungi.info.html",
    )
    parser.add_argument(
        "-f", "--file", help="Output file name. Will print to stdout if not provided."
    )
    args = parser.parse_args()
    main(args)
