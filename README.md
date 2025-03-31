# Fungal metatranscriptomics workflow

This is a snakemake workflow originally developed for the N_Street_1801
(Functional insights into the impacts of nitrogen fertilisation on forest
belowground fungal metacommunities) project.

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

Each profile has a subdirectory in the repository root containing a `config.yaml` file in YAML format:

```
./
├── dardel
│   └── config.yaml
└── local
    └── config.yaml
```

For both these profiles, the first section of the `config.yaml` file specifies several command line flags for Snakemake, *e.g.*:

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
repository root. The best practice is to make a copy of this file and apply any
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
JGI Mycocosm transcripts) or if the `extra_genomes:` is set to a non-empty
string (pointing to a list of genomes to use for generating a custom taxonomy
database, see below under [Generating a custom taxonomy
database](#generating-a-custom-taxonomy-database)) this requires that files are
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

Information about genomes found in the JGI Mycocosm database will be stored in `resources/JGI/genomes.tsv`. This file will be available after a dry-run of the workflow so you may inspect it before running the workflow. The file should look like this:

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

Assembled transcripts are annotated taxonomically using MMSeqs2 using one of the
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

When the `sample_file_list` contains accession ids for downloading directly from a public read archive, the `datadir` parameter sets the directory where fastq files are stored.

The `read_length` parameter is used by the STAR aligner when generating the host genome index. If your samples have different read lengths, set this to the maximum among the samples.