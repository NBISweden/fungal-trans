#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser
import sys

def main():
    parser = ArgumentParser()
    parser.add_argument("infile")
    args = parser.parse_args()
    df = pd.read_csv(args.infile, sep="\t", usecols=[0,2,10,11], names=["query","pid","evalue","score"])
    df.to_csv(sys.stdout, sep="\t", index=False)

if __name__ == '__main__':
    main()
