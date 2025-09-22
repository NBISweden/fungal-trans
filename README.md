# Fungal metatranscriptomics workflow

This is a snakemake workflow originally developed for the N_Street_1801
(Functional insights into the impacts of nitrogen fertilisation on forest
belowground fungal metacommunities) project.

## Table of contents

- [Installation](#installation)
- [Command line interface](#command-line-interface)
  - [Running on the Dardel HPC system](#running-on-the-dardel-hpc-system)
- [Configuration](#configuration)
  - [Essential parameters](#essential-parameters)
    - [Sample list](#sample-list)
    - [JGI login info](#jgi-login-info)
    - [Host filtering](#host-filtering)
    - [Generating a custom taxonomy database](#generating-a-custom-taxonomy-database)
  - [Additional parameters](#additional-parameters)
    - [Data paths and sample info](#data-paths-and-sample-info)

## Installation
To install this workflow do the following:

1. Clone the repository

```
git clone https://github.com/NBISweden/fungal-trans.git
```

2. Install [pixi](https://pixi.sh/latest).

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

3. Install and activate the environment with pixi:

```
pixi shell
```

The `pixi shell` command needs to be run first any time you want to run the
workflow.

## Command line interface

The workflow is executed using the `snakemake` command with the following syntax:

```bash
snakemake --configfile <your-configfile> --profile <profile>
```

Here, `--configfile` specifies a configuration file to use which contains parameters that modify the behaviour of the workflow (*e.g.* what input files to use, see below under [Configuration](#configuration)). The `--profile` flag specifies a configuration profile which modifies the behaviour of snakemake itself (*e.g.* number of cores and how to handle software requirements). There are two profiles included in this repository: 

- The `dardel` profile is tailored for running on the [Dardel](https://www.pdc.kth.se/hpc-services/computing-systems/dardel-hpc-system) HPC system.
- The `local` profile is more general and can be used to test the workflow on your local computer.
- The `kebnekaise` profile is tailored for running on the [HPC2N](https://www.hpc2n.umu.se/) system.

Each profile has a subdirectory in the repository root containing a `config.yaml` file in YAML format:

```
./
├── dardel
│   └── config.yaml
├── kebnekaise
│   └── config.yaml
└── local
    └── config.yaml
```

For all these profiles, the first section of the `config.yaml` file specifies several command line flags for Snakemake, *e.g.*:

```yaml
rerun-triggers: "mtime" # use modification time of files as a trigger for rerunning jobs
keep-going: True # keep the workflow running as long as possible even if one job fails
printshellcmds: True # print shell commands to stdout
cores: 1 # cores to use
software-deployment-method: "apptainer" # manage software requirements with apptainer
```

These options can be overridden by flags you specify on the command line, so for
example to run with the `local` profile but give Snakemake more cores you can
run:

```bash
snakemake --configfile <your-config-file> --profile local --cores 4
```

or to use Conda instead of Apptainer to handle software dependencies:

```bash
snakemake --configfile <your-config-file> --profile local --cores 4 --software-deployment-method conda
```

### Running on the Dardel HPC system

1. Specify compute account

The `dardel/config.yaml` file contains rule-specific resource settings tailored
for the Dardel HPC system. The most important thing you will have to modify in
this file is the `slum_account:` setting under the `default-resources:` section.
By default this section looks like this:

```yaml
default-resources: 
  slurm_account: <your-slurm-account>
  slurm_partition: shared
  runtime: 120
```

You will have to change the `<your-slurm-account>` part to match the compute
account you want to use. This is typically in the form of `naissYYY-NN-NNNN`.

2. Make sure `apptainer` is in your path

By default the `dardel` profile is set up to use
[Apptainer](https://apptainer.org/) to run most rules in isolated containers. For this to work you must either have installed `apptainer` somewhere on the Dardel system and added it to your `$PATH`, or you can simply run the following after you activate the pixi shell (see point 3 under [Installation](#installation) above):

```bash
module load PDC apptainer
```

## Configuration

The workflow can be configured using a configuration file in YAML format. The
default configuration file is included in this repo as `config.yml` in the
repository root. We recommend making a copy of this file and applying any
changes needed in the copy, _e.g._:

```bash
cp config.yml myconfig.yml
```

### Essential parameters

#### Sample list

The workflow requires a list of samples to use as input. By default the workflow
looks for a file called `sample_list.tsv` in the working directory but this can
be configured using the `sample_file_list` parameter in the configuration file.

The sample list should be tab-separated and contain a header as the first line.
An example sample list is shown below:

| Sample | R1 | R2 | assembly |
| ------ | -- | -- | -------- |
| sample1 | data/sample1_R1.fastq.gz | data/sample1_R2.fastq.gz | co-assembly |
| sample2 | data/sample2_R1.fastq.gz | data/sample2_R2.fastq.gz | co-assembly |

In this example there are two samples `sample1` and `sample2`. The `Sample`
column (required) specifies the name of the samples which will be used in the
workflow output. The `R1` and `R2` columns specify the path (absolute or
relative) to the R1 and R2 fastq files, respectively. The `assembly` column
(optional) specifies the co-assemblies to generate. If this column is omitted
then no co-assemblies will be generated, even if you've set `co_assembly: True`
in the configuration file (see below). You can exclude the `assembly` column and
set `single_assembly: True` in the configuration file which will generate one
individual assembly per sample.

If you have data in a public sequencing archive, such as SRA, then you can
exclude the `R1` and `R2` columns and instead add a column called `accession`
which lists the archive accession (one per sample). The workflow will then
download the R1 and R2 files and store them in the directory specified by the
`datadir` config parameter.

#### JGI login info

If your config file has `fungi_filter: True` (to filter reads by mapping against
JGI Mycocosm transcripts) or if the `extra_genomes:` parameter is set to a non-empty
string (pointing to a list of genomes to use for generating a custom taxonomy
database, see below under [Generating a custom taxonomy
database](#generating-a-custom-taxonomy-database)), this requires that files are
downloaded from [JGI Mycocosm](https://mycocosm.jgi.doe.gov/mycocosm/home). To
make this work you need to have login credentials which you can obtain by
registering [here](https://contacts.jgi.doe.gov/registration/new). Once you have
your password, add your email adress and password to a YAML file like in the
example below:

```yaml
jgi_user: "your.email@adress.com"
jgi_password: "your-password"
```

then modify your configuration file so that the `jgi_account_info` parameter
points to the file with your credentials.

The `fungi_info` parameter specifies the path to a JGI Genome portal listing in
HTML format. Previously, this file could be parsed directly from the URL
https://mycocosm.jgi.doe.gov/fungi/fungi.info.html but this is not possible
anymore because of a CAPTCHA introduced by JGI. Instead you can download this
file in your browser and supply the HTML file with the `fungi_info` parameter.

The `fungi_genomes_file` specifies where to save information for genomes parsed
from the `fungi_info` HTML file. This file will be available after a dry-run of the
workflow so you may inspect it before running the workflow. The file should look
like this:

| portal | Name | bp | genes |
|--------|------|----|------|
| Aaoar1 | Aaosphaeria arxii CBS 175.79 v1.0 | 38901049 | 14203 |
| Abobi1 | Abortiporus biennis CCBS 521 v1.0 | 45165060 | 11987 |
| Abobie1 | Abortiporus biennis ​CIRM-BRFM 1778 v1.0 | 33118568 | 11767 |

The `portal` column contains the JGI portal name, `Name` is the name of the genome, `bp` is the size of the genome in base pairs and `genes` is the number of genes in the genome. You may filter this file to your liking to only download transcripts for the genomes you are interested in. Remember to save it as `resources/JGI/genomes.tsv` before running the workflow.

#### Host filtering

In order to filter host reads from your input a host genome fasta file and GFF
file needs to be specified in your configuration file with the `host_fna` and
`host_gff` parameters, respectively. Both these files should be gzipped.

#### Generating a custom taxonomy database

Assembled transcripts are annotated taxonomically with MMSeqs2 using one of the
official databases with taxonomy support (see configuration parameters
`mmseqs_db_dir` and `mmseqs_db`). However, if you want to supplement the
taxonomy database with additional reference proteins from JGI Mycocosm you need
to supply a list of these genomes in JGI Mycocosm. The file should be
tab-separated and have a header as the first row. The first column should
contain the portal name in JGI Mycocosm and the second column should have the
NCBI taxonomy id of the genome. An example file is included in the workflow as
`extra_JGI_genomes.tsv` in the root of the repository. The file of extra genomes
to use can be specified with the `extra_genomes` config parameter.

### Additional parameters

Below are additional parameters which do not *require* modification to be able
to run the workflow on your data.

#### Data paths and sample info

```yaml
datadir: "data/raw"
read_length: 150
```

When the `sample_file_list` contains accession ids for downloading directly from
a public read archive, the `datadir` parameter sets the directory where fastq
files are stored.

The `read_length` parameter is used by the STAR aligner when generating the host
genome index. If your samples have different read lengths, set this to the
maximum among the samples.

#### Preprocessing

```yaml
fastp_adapter_sequence: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
fastp_adapter_sequence_r2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
fastp_length_required: 50
```

The parameters `fastp_adapter_sequence` and `fastp_adapter_sequence_r2` specify
the adapter sequences to trim from the R1 and R2 files, respectively. The
default values shown above are for Illumina TruSeq adapters.

The `fastp_length_required` parameter sets the minimum read length required to
keep reads after quality/adapter trimming.

#### Filtering

Filtering of reads can be done in various ways and is applied after
adapter/quality trimming and removal of rRNA reads.

```yaml
host_filter: False
host_fna: ""
host_gff: ""
star_limitGenomeGenerateRAM: 128
star_extra_build_params: "--genomeChrBinNbits 15 --genomeSAsparseD 3 --genomeSAindexNbases 14 --sjdbGTFtagExonParentTranscript Parent"
star_extra_params: "--outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --seedSearchStartLmax 50"

fungi_filter: False
paired_strategy: "one_mapped"
fungi_info: "https://mycocosm.jgi.doe.gov/fungi/fungi.info.html"
fungi_genomes_file: "resources/JGI/genomes.tsv"
jgi_account_info: "config/jgi_account_info.yml"
strobealign_strobe_len: 17

kraken_filter: False
kraken_db: "20240229-015730_nt"
kraken_db_dir: "/sw/data/Kraken2_data"
kraken_confidence: 0.2
```

The `host_filter` parameter indicates whether reads that align to a host genome
should be removed. If `host_filter` is set to `True` then the parameters
`host_fna` and `host_gff` **must** point to a gzipped fasta and GFF annotation
file, respectively. Reads are aligned using the
[STAR](https://github.com/alexdobin/STAR) aligner. 

The `star_limitGenomeGenerateRAM` parameter sets the maximum allowed RAM usage
(in GB) when generating the host genome index. The `star_extra_build_params` and
`star_extra_params` parameters act as catch-all strings for any settings you
want to pass to the STAR build and STAR align rules.

Setting `fungi_filter` parameter to `True` activates filtering of preprocessed
reads by mapping reads to a collection of transcriptomes downloaded from the JGI
Genome Portal. The `paired_strategy` parameter sets how read fungal filtering is
applied. With `one_mapped` paired end reads where *at least* one of the reads
map to a fungal transcriptome are kept as fungal. With `both_mapped` only reads
where both reads in a pair map to fungi are kept. The `strobealign_strobe_len`
parameter sets the seed length and is used by
[strobealign](https://github.com/ksahlin/StrobeAlign) when filtering fungal
reads. For additional information see the [JGI Login Info](#jgi-login-info)
section.

With `kraken_filter` set to `True` reads are also assigned a taxonomy using kraken2. The kraken database to use is specified with the `kraken_db` and `kraken_db_dir` parameters. For example, if you have the full `NT` kraken database stored in `/data/Kraken2_data/k2_core_nt_20250609/hash.` in the form:

```
/data/Kraken2_data
└── k2_core_nt_20250609
    ├── hash.k2d
    ├── opts.k2d
    └── taxo.k2d
```

then `kraken_db_dir` should be set to `data/Kraken2_data` and `kraken_db` to `k2_core_nt_20250609`.

#### Assembly

The workflow assembles preprocessed and filtered reads with either Trinity,
TransAbyss or Megahit. We recommend Trinity which is also the default.

```yaml
co_assembly: True
single_assembly: False
assembler: trinity
transabyss_kmers: [21, 25, 29, 31, 65, 75, 85, 95]
min_contig_len: 200
insilico_norm_max_cov: 30
insilico_norm_mem: 100
```

With `co_assembly` set to `True`, co-assemblies are generated based on the
`assembly` column in the sample information file. If `single_assembly` is set to
`True` then each sample is assembled separately.

The `assembler` parameter sets what assembly software to use. Choose from
`trinity` (recommended), `transabyss` or `megahit`.

The `transabyss_kmers` parameter is a list of kmer lengths to use in an
iterative assembly with TransAbyss.

The `min_contig_len` parameter specifies the minimum contig length to keep from
assemblies.

When generating co-assemblies, preprocessed reads from all input samples are
passed through an *in silico* normalization step prior to assembly. This is done
for all three supported assembly tools. The `insilico_norm_max_cov` parameter
sets the targeted maximum coverage for reads. The `insilico_norm_mem` parameter
sets the allowed memory usage (in GB) for the normalization step.

#### Read mapping/quantification

Assembled transcripts and predicted proteins are quantified by mapping
preprocessed and filtered reads against the assemblies. The workflow runs both
an alignment-based mapping step using the RSEM program as well as a
pseudo-alignment step using Kallisto. In addition, the BAM file generated as
part of running RSEM is used with featureCounts to assign reads directly to
features in the assembly.

```yaml
fc_params: "-p --countReadPairs -t exon -g Parent -Q 10 -s 0 "
kallisto_params: "-b 100"
```

The `fc_params` parameter passes settings to featureCounts. The default is to:

- run in paired end mode (-p)
- count fragments (--countReadPairs)
- assign reads to exon features (-t exon)
- group by Parent attributes (-g Parent)
- require a minimum mapping quality of 10 (-Q 10), and
- set strandedness to 0 (-s 0) (standard Illumina)

> [!IMPORTANT]  
> If your input data is from IlluminaTruSeq Stranded libraries, add `-s 2` to
> the `fc_params` string and `--rf-stranded` to the `kallisto_params`.