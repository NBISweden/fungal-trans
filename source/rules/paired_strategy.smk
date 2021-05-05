if config["paired_strategy"] == "concordant":
    localrules:
        collate_fungi,
        collate_spruce

    rule collate_fungi:
        input:
            "results/bowtie2/{sample_id}/{sample_id}_{pair}.fungi.conc.fastq.gz"
        output:
            "results/bowtie2/{sample_id}/{sample_id}_{pair}.fungi.fastq.gz"
        shell:
            """
            mv {input} {output} 
            """
    rule collate_spruce:
        input:
            "results/bowtie2/{sample_id}/{sample_id}_{pair}.fungi.nospruce.conc.fastq.gz"
        output:
            "results/bowtie2/{sample_id}/{sample_id}_{pair}.fungi.nospruce.fastq.gz"
        shell:
            """
            mv {input} {output}
            """

elif config["paired_strategy"] == "both-mapped":
    rule collate_fungi:
        """Extracts paired reads if both ends are mapped"""
        input:
            "results/bowtie2/{sample_id}/{sample_id}.fungi.bam"
        output:
            R1 = "results/bowtie2/{sample_id}/{sample_id}_R1.fungi.fastq.gz",
            R2 = "results/bowtie2/{sample_id}/{sample_id}_R2.fungi.fastq.gz"
        params:
            R1 = "$TMPDIR/{sample_id}_R1.fungi.fastq",
            R2 = "$TMPDIR/{sample_id}_R2.fungi.fastq"
        resources:
            runtime = lambda wildcards, attempt: attempt**2*30
        shell:
            """
            samtools fastq -F 12 -1 {params.R1} -2 {params.R2} -s /dev/null -0 /dev/null {input}
            gzip {params.R1} {params.R2}
            mv {params.R1}.gz {output.R1}
            mv {params.R2}.gz {output.R2}
            """
    rule collate_spruce:
        """Extracts paired reads where neither mate is mapped"""
        input:
            "results/bowtie2/{sample_id}/{sample_id}.spruce.bam"
        output:
            R1="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.nospruce.fastq.gz",
            R2="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.nospruce.fastq.gz"
        params:
            R1="$TMPDIR/{sample_id}_R1.fungi.nospruce.fastq",
            R2="$TMPDIR/{sample_id}_R2.fungi.nospruce.fastq",
            only_this_end="$TMPDIR/{sample_id}.fungi.nospruce.only_this_end.bam",
            only_that_end="$TMPDIR/{sample_id}.fungi.nospruce.only_that_end.bam",
            both_unmapped="$TMPDIR/{sample_id}.fungi.nospruce.both_unmapped.bam",
            merged="$TMPDIR/{sample_id}.fungi.nospruce.merged.bam"
        resources:
            runtime = lambda wildcards, attempt: attempt**2*30
        shell:
            """
            # Get this end unmapped, other end mapped
            samtools view -b -f 4 -F 8 {input} > {params.only_that_end}
            # Get this end mapped, other end unmapped
            samtools view -b -F 4 -f 8 {input} > {params.only_this_end}
            # Get both reads unmapped
            samtools view -b -f 12 {input} > {params.both_unmapped}
            # Merge bam files
            samtools merge {params.merged} {params.only_that_end} {params.only_this_end} {params.both_unmapped}
            # Output fastq
            samtools fastq -1 {params.R1} -2 {params.R2} -0 /dev/null {params.merged}
            gzip {params.R1} {params.R2}
            mv {params.R1}.gz {output.R1}
            mv {params.R2}.gz {output.R2}
            rm {params.merged} {params.both_unmapped} {params.only_this_end} {params.only_that_end}
            """
else:
    rule collate_fungi:
        """Extracts paired reads from bam file if either end is mapped"""
        input:
            "results/bowtie2/{sample_id}/{sample_id}.fungi.bam"
        output:
            R1="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.fastq.gz",
            R2="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.fastq.gz"
        params:
            R1 = "$TMPDIR/{sample_id}_R1.fungi.fastq",
            R2 = "$TMPDIR/{sample_id}_R2.fungi.fastq",
            only_this_end = "$TMPDIR/{sample_id}_onlythisend.bam",
            only_that_end = "$TMPDIR/{sample_id}_onlythatend.bam",
            bothends = "$TMPDIR/{sample_id}_bothends.bam",
            merged = "$TMPDIR/{sample_id}_merged.bam"
        resources:
            runtime = lambda wildcards, attempt: attempt**2*60
        shell:
            """
            # Get this read with mate unmapped
            samtools view -b -F 4 -f 8 {input} > {params.only_this_end}
            # Get unmapped reads with mate mapped
            samtools view -b -f 4 -F 8 {input} > {params.only_that_end}
            # Get both reads mapped
            samtools view -b -F 12 {input} > {params.bothends}
            # Merge bam file
            samtools merge {params.merged} {params.only_this_end} {params.only_that_end} {params.bothends}
            # Extract fastq from merged file 
            samtools fastq -1 {params.R1} -2 {params.R2} -0 /dev/null {params.merged}
            gzip {params.R1} {params.R2}
            mv {params.R1}.gz {output.R1}
            mv {params.R2}.gz {output.R2}
            rm {params.merged} {params.bothends} {params.only_that_end} {params.only_this_end}
            """
    rule collate_spruce:
        """Extracts only unmapped reads with no mapped mate"""
        input:
             "results/bowtie2/{sample_id}/{sample_id}.spruce.bam"
        output:
             R1="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.nospruce.fastq.gz",
             R2="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.nospruce.fastq.gz"
        resources:
            runtime = lambda wildcards, attempt: attempt**2*30
        params:
            R1="$TMPDIR/{sample_id}_R1.fungi.nospruce.fastq",
            R2="$TMPDIR/{sample_id}_R2.fungi.nospruce.fastq"
        shell:
            """
            # Extract only reads where neither read in a pair is mapped
            samtools fastq -f 12 -1 {params.R1} -2 {params.R2} -s /dev/null -0 /dev/null {input}
            gzip {params.R1} {params.R2}
            mv {params.R1}.gz {output.R1}
            mv {params.R2}.gz {output.R2}
            """