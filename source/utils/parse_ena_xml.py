#!/usr/bin/env python

import xml.etree.ElementTree as ET
from argparse import ArgumentParser
import os
import pandas as pd
import sys


def get_root(f):
    tree = ET.parse(f)
    root = tree.getroot()
    return root


def extract_files(f):
    '''
    Extracts XML REFNAME values and associated fastq files
    :param f: ENA XML file (*.Run.xml)
    :return: dictionary with XML refname -> filename
    '''
    root = get_root(f)
    files = {}
    for r in root.findall('RUN'):
        refname = r[0].attrib['refname']
        R1=r.findall('DATA_BLOCK')[0].getchildren()[0].getchildren()[0].attrib['filename']
        R2=r.findall('DATA_BLOCK')[0].getchildren()[0].getchildren()[1].attrib['filename']
        files[refname] = {'Read_file': os.path.basename(R1),
                         'Pair_file': os.path.basename(R2)}
    df = pd.DataFrame(files).T
    return df


def extract_accessions(f):
    '''
    Extracts alias and accession values from XML receipt
    :param f:
    :return:
    '''
    root = get_root(f)
    accessions = {}
    for e in root.findall('EXPERIMENT'):
        accessions[e.attrib['alias']] = e.attrib['accession']
    df = pd.DataFrame(accessions, index=['accession']).T
    return df


def extract_samples(f):
    """
    Extracts ENA experiment REF to ENA accession
    :param f:
    :return:
    """
    r = get_root(f)
    names = {}
    for c in r.getchildren():
        refname = c.attrib['alias']
        for d in c.findall('DESIGN'):
            lib_name = d.findall('LIBRARY_DESCRIPTOR')[0].findall('LIBRARY_NAME')[0].text
            names[refname] = lib_name
    df = pd.DataFrame(names, index=['Sample']).T
    return df


def main(args):
    files = extract_files(args.runxml)
    accessions = extract_accessions(args.receiptxml)
    samples = extract_samples(args.experimentxml)
    df = pd.merge(files, accessions, left_index=True, right_index=True)
    df = pd.merge(df, samples, left_index=True, right_index=True)
    df = df.reset_index().loc[:,['Sample', 'accession', 'index', 
                                 'Read_file', 'Pair_file']]
    df.set_index("Sample", inplace=True)
    df.to_csv(sys.stdout, sep='\t')

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('runxml', help='XML file with RUN information')
    parser.add_argument('receiptxml', help='XML receipt')
    parser.add_argument('experimentxml', help='XML file with EXPERIMENT info')
    args = parser.parse_args()
    main(args)
