# Test data

To evaluate the filtering performance known fungal/host reads were generated
using `randomreads.sh` from the `bbmap` package:

```bash
conda create -n bbmap -c bioconda bbmap
conda activate bbmap
```

For fungi a subsampled fasta file of fungal transcripts downloaded from JGI was
used:

```bash
randomreads.sh ref=data/fungi_transcripts.100k.fasta \
  out1=data/fungi_reads_R1.fastq \
  out2=data/fungi_reads_R2.fastq \
  simplenames=t length=101 prefix=fungi \
  reads=100 paired=t illuminanames=t minq=30
```

The same was done for host reads using a subsampled fasta file of the 
_Picea abies_ genome (subsampled to 100 random sequences):

```bash
randomreads.sh ref=data/host.sampled100.fna \
  out1=data/host_reads_R1.fastq \
  out2=data/host_reads_R2.fastq \
  simplenames=t length=101 prefix=host \
  reads=100 paired=t illuminanames=t minq=30
```

A `host` prefix was added to these reads:

```bash
seqtk rename temp_host_reads_R1.fastq host > host_reads_R1.fastq
seqtk rename temp_host_reads_R2.fastq host > host_reads_R2.fastq
```