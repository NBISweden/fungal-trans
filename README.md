# Fungal metatranscriptomics workflow

This is a snakemake workflow originally developed for the N_Street_1801 (Functional 
insights into the impacts of nitrogen fertilisation on forest belowground
 fungal metacommunities) project.

## Installation
To install this workflow:

1. Clone the repository

```
git clone https://github.com/NBISweden/fungal-trans.git
```

2. Install [pixi](https://pixi.sh/latest).

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

3. Install base environment with pixi:

```
pixi shell
```

## Configuration

The workflow can be configured using a configuration file in YAML format. The
default configuration file is included in this repo as `config.yml` in the
repository root. The best practice is to make a copy of this file and apply any
changes needed in the copy, _e.g._:

```bash
cp config.yml myconfig.yml
```

### Sample data
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

If you have data in a public sequencing archive, such as SRA, then you can exclude the `R1` and `R2` columns and instead add a column called `accession` which lists the archive accession (one per sample). The workflow will then download the R1 and R2 files and store them under 

### Option1: Raw data on disk
Let's say you have your raw data downloaded and stored in a directory
called `/home/data/raw` and your study contains 3 samples called 
sample1, sample2 and sample3. Let's also say you want to make a co-assembly
of samples 2 and 3 and assemble sample 1 separately. Then your **tab-separated**
 sample file list should look like this:  

| sample | Read_file | Pair_file | assembly |
|--------|-----------|-----------|----------|
|sample1|sample1_1.fastq.gz|sample1_2.fastq.gz|assembly1|
|sample2|sample2_1.fastq.gz|sample2_2.fastq.gz|assembly2|
|sample3|sample3_1.fastq.gz|sample3_2.fastq.gz|assembly2|

In this example, we'll name sample file list `sample_list.tsv`. You can
then run the workflow for these samples with:
```
snakemake -p --use-conda -n
```

See below for more information on how to **configure** and **run** the
workflow.

### Option2: Data in a sequence read archive

If you have data in a public data repository such as the [SRA](https://www.ncbi.nlm.nih.gov/sra)
then you can add the accession numbers directly to the sample list:

| Sample  |  accession  | assembly  |
| ------- | ----------- | --------- |
| needle | ERX3761257  | needle_asm |
| root | ERX3761470  | root_asm |

This will make the workflow download the fastq data using `sra-tools` prior
to starting.

## Running the workflow

### Configuration
By default the workflow uses the configuration file `config.yml`.
You can either make changes in that config file or copy the file and make
your changes in the copy (_e.g._ `cp config.yml myconfig.yml`). To run the 
workflow with another config file, lets say `myconfig.yml`, specify 
`--configfile myconfig.yml` on the snakemake command line.

### Running the full workflow
To get an overview of the jobs that will be executed do:
```
snakemake --use-conda -j 4 -np
```
`-n` executes a 'dry-run' that doesn't actually do anything and instead
the jobs to be run are printed to stdout (`-p` also prints the bash commands). 
The `-j 4` parameter tells snakemake to run with 4 cores.

### Running parts of the workflow
The workflow is divided into steps that you can run separately. Below
are descriptions of these targets, their output and how to run them.

#### Preprocessing

The preprocessing part will download fastq files for your samples, unless you
start with files locally on disk as described above. Adapter trimming is 
performed using cutadapt followed by quality trimming with trimmomatic. FastQC
is then run on the trimmed sequences and summarized into a report using MultiQC.

Here's an example of running the preprocessing part of the workflow: 
```
snakemake --use-conda -j 10 -p preprocess
```

#### Filter

Preprocessed reads are filtered to separate fungal and host reads. This is done
by first mapping reads to fungal transcripts, thus generating files with 
putative fungi and non-fungi reads. The putative fungal reads are then mapped
to a host database of choice, _e.g._ spruce, either using `bowtie2` or `STAR`. 
Reads mapping to the host are combined with the non-fungal reads from the first
filtering step and put into a 'host' sequence file under 
`results/host/{sample_id}_R{1,2}.fastq.gz`. Reads that do not map to the host db
at this stage are used for downstream assembly and annotation.

The type of filtering performed is determined by the `read_source` config parameter.

- With `read_source: bowtie2` filtering is done as described as above.
- With `read_source: taxmapper` fungal reads are identified using [Taxmapper](https://bitbucket.org/dbeisser/taxmapper/src/master/).
- With `read_source: filtered` both bowtie2 and taxmapper are used to identify 
fungal reads and the union of reads identified by these methods are used for 
downstream analyses.
- With `read_source: unfiltered` no filtering is performed and preprocessed reads 
  are directly used for downstream analyses.

Here's an example of running up to and including the filtering part of the workflow: 
```
snakemake --use-conda -j 10 -p filter
```

##### Configuring login to JGI server

The workflow automatically downloads transcript data from the [JGI Mycocosm](https://mycocosm.jgi.doe.gov/fungi/fungi.info.html) database. For this to work you must have [registered an account](https://contacts.jgi.doe.gov/registration/new). After registering, create a file in YAML format containing your login credentials:

```yaml
jgi_password: <your_password>
jgi_user: <your_username/email>
```

Save this file as `jgi_login.yml` and pass it to the workflow with the `--configfile` flag:

```bash
snakemake --configfile jgi_login.yml <other flags>
```

Information about genomes found in the JGI Mycocosm database will be stored in `resources/JGI/genomes.tsv`. This file will be available after a dry-run of the workflow so you may inspect it before running the workflow. The file should look like this:

| portal | Name | bp | genes |
|--------|------|----|------|
| Aaoar1 | Aaosphaeria arxii CBS 175.79 v1.0 | 38901049 | 14203 |
| Abobi1 | Abortiporus biennis CCBS 521 v1.0 | 45165060 | 11987 |
| Abobie1 | Abortiporus biennis ​CIRM-BRFM 1778 v1.0 | 33118568 | 11767 |

The `portal` column contains the JGI portal name, `Name` is the name of the genome, `bp` is the size of the genome in base pairs and `genes` is the number of genes in the genome. You may filter this file to your liking to only download transcripts for the genomes you are interested in. Remember to save it as `resources/JGI/genomes.tsv` before running the workflow.


#### Assembly

Assemblies can be generated for single samples or by combining multiple samples
into co-assemblies. You can choose to assemble with either [Trinity](https://github.com/trinityrnaseq/trinityrnaseq),
[Trans-ABySS](https://github.com/bcgsc/transabyss) or [Megahit](https://github.com/voutcn/megahit)
assemblers, configurable via the `assembler` config parameter. 

When `single-assembly` is set to `True` in the config, preprocessed and filtered 
reads from each sample are individually assembled using the configured assembler.

When `co-assembly` is set to `True`, the `sample_file_list` must contain an 
'assembly' column that groups samples into co-assembly groups, _e.g._:

| sample | Read_file | Pair_file | assembly |
|--------|-----------|-----------|----------|
|sample1|sample1_1.fastq.gz|sample1_2.fastq.gz|assembly1|
|sample2|sample2_1.fastq.gz|sample2_2.fastq.gz|assembly2|
|sample3|sample3_1.fastq.gz|sample3_2.fastq.gz|assembly2|

In this example, preprocessed and filtered reads from sample2 and sample3 will be
combined into a co-assembly named `assembly2`. Note that prior to generating 
co-assemblies, reads are deduplicated using `fastuniq`.
 
Here's an example of running up to and including the assembly part of the workflow: 
```
snakemake --use-conda -j 10 -p assemble
```

To generate co-assemblies:

```
snakemake --use-conda -j 10 -p co_assemble
```

#### Annotation

Gene calling is done on assemblies using [GeneMarkS-T](http://topaz.gatech.edu/GeneMark/license_download.cgi).
Translated amino acid sequences are then annotated functionally using [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper) and 
the [dbCAN database](https://bcb.unl.edu/dbCAN/), and taxonomically using 
[contigtax](https://github.com/NBISweden/contigtax).

Here's an example of running up to and including the annotation part of the workflow: 
```
snakemake --use-conda -j 10 -p annotate
```

And to annotate co-assemblies:

```
snakemake --use-conda -j 10 -p annotate_co
```

### Running the workflow on a compute cluster

The workflow comes with support for job-handling on compute clusters with the
SLURM workload manager. Simply set your SLURM account id in the `config/cluster.yaml`
file:

````yaml
__default__:
  account: staff # <- REPLACE 'staff' with your account id
````

Then you can run the workflow as:

```
snakemake --use-conda --profile slurm -j 10 -p
```

Another option is to simply submit the entire workflow (in whole or in part) as
one big SLURM job. The benefit of using `--profile slurm` is that resources are
handled in a more fine grained way and that failed jobs can be resubmitted 
automatically with longer run times.