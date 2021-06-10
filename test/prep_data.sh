#!/bin/bash

for sample in needle root;
do
  for read in 1 2;
  do
    cat test/data/fungi_reads_R${read}.fastq test/data/host_reads_R${read}.fastq test/data/${sample}_${read}.fastq | gzip -c > test/data/${sample}_${read}.fastq.gz
  done
done