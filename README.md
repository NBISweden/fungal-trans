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

2. Install the necessary dependencies with [conda](https://docs.conda.io/en/latest/miniconda.html).

```
conda env create -f environment.yml
```

**Remember to activate the environment**

3. Install GeneMarkS-T gene caller

The [GeneMarkS-T](http://exon.gatech.edu/GeneMark/) gene caller software used in this workflow can be 
downloaded for Linux from [this download page](http://topaz.gatech.edu/GeneMark/license_download.cgi).
Activate the workflow environment created in the step above, then you can install 
GeneMarkS-T by extracting the archive directly into the environment path:

```bash
tar -C $CONDA_PREFIX/bin/ -xvf gmst_linux_64.tar.gz
```

**IMPORTANT**: Because the workflow relies on conda to handle dependencies for
most steps it is important to always add the `--use-conda` flag to snakemake. 
Even better is actually to install `mamba` (a more efficient replacement for
conda) and have snakemake use that instead by adding `--conda-frontend mamba`
to the command line.

## Sample data
The workflow requires a list of samples to use as input. By default the 
workflow looks for a file called `sample_list.tsv` in the working directory
 but you can name it whatever you like and specify the path with 
 `snakemake --config sample_file_list=<path-to-your-samplefile>`. Data for these
samples can either be present on disk or in a remote sequence read archive.

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
snakemake -j 10 -p preprocess
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
snakemake -j 10 -p filter
```

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
 
Here's an example of running the preprocessing part of the workflow: 
```
snakemake -j 10 -p assemble
```

To generate co-assemblies:

```
snakemake -j 10 -p co_assemble
```

#### Annotation

Gene calling is done on assemblies using [GeneMarkS-T](http://topaz.gatech.edu/GeneMark/license_download.cgi).
Translated amino acid sequences are then annotated using [eggnog-mapper] 