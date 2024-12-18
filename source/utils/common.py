#!/usr/bin/env

import pandas as pd


def parse_sample_list(f, config):
    """
    Parse the sample list and set up input file names for the assembly
    based on configuration
    """
    samples = {}
    map_dict = {}
    if config['filter_reads']:
        input_dir = 'filtered'
    else:
        input_dir = 'unfiltered'
    host_aligner = config["host_aligner"]
    dir = config['datadir']
    suffices = {'unfiltered': 'cut.trim.fastq.gz',
                'filtered': 'filtered.union.fastq.gz',
                'taxmapper': 'cut.trim.filtered.fastq.gz',
                'bowtie2': 'fungi.nohost.fastq.gz',
                'star': 'fungi.nohost.fastq.gz'}
    # Read sample list
    df = pd.read_csv(f, comment='#', header=0, sep='\t', index_col=0, dtype=str)
    df.fillna('', inplace=True)
    # If it's just a one-column file, expand it by assuming that the
    # sra accessions are in the first column
    if df.shape[1] == 0 or 'accession' not in df.columns:
        df = df.assign(accession=pd.Series(df.index, index=df.index))
        df.index.name = "Sample"
    assemblies = {}
    if 'assembly' in df.columns:
        for assembly in df.assembly.unique():
            assemblies[assembly] = {'R1': [], 'R2': []}
    for sample in df.index:
        try:
            R1 = df.loc[sample,'Read_file']
            R2 = df.loc[sample,'Pair_file']
        except KeyError:
            R1 = '{}_1.fastq.gz'.format(sample)
            R2 = '{}_2.fastq.gz'.format(sample)
        if 'accession' in df.columns:
            accession = df.loc[sample, 'accession']
        else:
            accession = ''
        if len(assemblies.keys()) > 0:
            if df.loc[sample,'assembly'] != '':
                assembly = df.loc[sample, 'assembly']
                # Define reads for assembly
                R1_f = f"results/{input_dir}/{sample}/{sample}_R1.fastq.gz"
                R2_f = f"results/{input_dir}/{sample}/{sample}_R2.fastq.gz"
                assemblies[assembly]['R1'].append(R1_f)
                assemblies[assembly]['R2'].append(R2_f)
                map_dict[sample] = {'R1': R1_f, 'R2': R2_f}

        samples[sample]={'R1': f'{dir}/{R1}',
                         'R2': f'{dir}/{R2}',
                         'accession': accession}
    return samples, map_dict, assemblies


def fungi_input(wc):
    suffices = {'unfiltered': 'cut.trim.fastq.gz',
                'filtered': 'filtered.union.fastq.gz',
                'taxmapper': 'cut.trim.filtered.fastq.gz',
                'bowtie2': 'fungi.nohost.fastq.gz'}
    if wc.filter_source == "bowtie2":
        aligndir = config["host_aligner"]
    else:
        aligndir = wc.filter_source
    R1 = f'results/{aligndir}/{wc.sample_id}/{wc.sample_id}_R1.{suffices[wc.filter_source]}'
    R2 = f'results/{aligndir}/{wc.sample_id}/{wc.sample_id}_R2.{suffices[wc.filter_source]}'
    return [R1, R2]

def parse_extra_genomes(f):
    """
    Parse the extra genomes file
    """
    extra_genomes = {}
    with open(f, 'r') as f:
        for line in f:
            line = line.strip()
            portal, taxid = line.rsplit()[0:2]
            try:
                taxid = int(taxid)
            except ValueError:
                continue
            extra_genomes[portal] = taxid
    return extra_genomes
