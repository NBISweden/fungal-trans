localrules:
    strobealign_collect_chunks,

#############
## MAPPING ##
#############

rule strobealign_map_fungi:
    """
    This rule maps preprocessed reads against fungal genome transcript files.
    This is done in chunks of 100 genomes to avoid memory issues.
    The output is a PAF file with the mapping results. Reads listed twice in the PAF
    file have both ends mapped (properly or not) and these are output to a separate
    file for each chunk.
    """
    output:
        paf="results/strobealign/{sample_id}/{chunk}/{sample_id}.fungi.paf.gz",
        both="results/strobealign/{sample_id}/{chunk}/{sample_id}.fungi.both_mapped.gz",
    input:    
        R1=rules.sortmerna.output.R1,
        R2=rules.sortmerna.output.R2,
        fna=lambda wildcards: strobealign_chunks[wildcards.chunk]
    log:
        "results/strobealign/{sample_id}/{chunk}/{sample_id}.strobealign.log"
    params:
        tmpdir = "$TMPDIR/{sample_id}.{chunk}",
        k = config["strobealign_strobe_len"]
    container: "docker://quay.io/biocontainers/strobealign:0.15.0--h5ca1c30_1"
    threads: 6
    shell:
        """
        mkdir -p {params.tmpdir}
        rm -f {params.tmpdir}/ref.fna
        for f in {input.fna} ; 
        do
            if [ -s $f ]; then
                gunzip -c $f >> {params.tmpdir}/ref.fna
            fi
        done
        strobealign -v --mcs -k {params.k} -x -t {threads} {params.tmpdir}/ref.fna {input.R1} {input.R2} 2>{log} | awk '$11>{params.k}' | igzip > {output.paf}
        igzip -c -d {output.paf} | cut -f1 | sort | uniq -d | igzip > {output.both}
        rm -r {params.tmpdir}
        """

rule strobealign_collect_chunks:
    """
    This rule collects the PAF output from the different chunks and concatenates them.
    """
    output:
        paf="results/strobealign/{sample_id}/{sample_id}.fungi.paf.gz",
        both="results/strobealign/{sample_id}/{sample_id}.fungi.both_mapped.gz",
    input:
        paf=expand("results/strobealign/{{sample_id}}/{chunk}/{{sample_id}}.fungi.paf.gz", chunk=strobealign_chunks.keys()),
        both=expand("results/strobealign/{{sample_id}}/{chunk}/{{sample_id}}.fungi.both_mapped.gz", chunk=strobealign_chunks.keys())
    shell:
        """
        cat {input.paf} > {output.paf}
        cat {input.both} > {output.both}
        """

rule process_fungal_mappings:
    """
    This rule extracts reads from the preprocessed fastq files that map to fungi.
    Reads with both ends and one or both ends mapped are extracted separately.
    """
    output:
        R1one="results/filtered/{sample_id}/{sample_id}_R1.fungi.one_mapped.fastq.gz",
        R2one="results/filtered/{sample_id}/{sample_id}_R2.fungi.one_mapped.fastq.gz",
        R1both="results/filtered/{sample_id}/{sample_id}_R1.fungi.both_mapped.fastq.gz",
        R2both="results/filtered/{sample_id}/{sample_id}_R2.fungi.both_mapped.fastq.gz",
    input:
        R1=rules.sortmerna.output.R1,
        R2=rules.sortmerna.output.R2,
        paf="results/strobealign/{sample_id}/{sample_id}.fungi.paf.gz",
        both="results/strobealign/{sample_id}/{sample_id}.fungi.both_mapped.gz",
    shell:
        """
        # Output reads that map with one or two ends to fungi
        seqkit grep -f <(gunzip -c {input.paf} | cut -f1) {input.R1} -o {output.R1one}
        seqkit grep -f <(gunzip -c {input.paf} | cut -f1) {input.R2} -o {output.R2one}
        # Output reads that map with both ends to fungi
        seqkit grep -f <(gunzip -c {input.both}) {input.R1} -o {output.R1both}
        seqkit grep -f <(gunzip -c {input.both}) {input.R2} -o {output.R2both}
        """

rule star_build_host:
    """
    Builds the host genome index for STAR mapping.
    """
    output:
        expand("resources/host/{f}",
            f = ["Genome", "SA", "SAindex", "chrLength.txt", "chrName.txt",
                 "chrNameLength.txt", "chrStart.txt", "genomeParameters.txt"])
    input:
        fna=config["host_fna"],
        gff=config["host_gff"]
    log:
        "resources/host/star_build_host.log"
    params:
        tmpdir="$TMPDIR/host",
        limitGenomeGenerateRAM = config["star_limitGenomeGenerateRAM"]*1000000000,
        sjdbOverhang=config["read_length"] - 1,
        extra_params = config["star_extra_build_params"],
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    threads: 20
    resources:
        runtime = 60 * 24,
        mem_mb = config["star_limitGenomeGenerateRAM"]*1000, # converts ram in GB to MB
    conda: "../../envs/star.yaml"
    container: "docker://quay.io/biocontainers/star:2.7.11b--h5ca1c30_5"
    shell:
        """
        exec &>{log}
        mkdir -p {params.tmpdir}
        gunzip -c {input.fna} > {params.tmpdir}/host.fna
        gunzip -c {input.gff} > {params.tmpdir}/host.gff
        STAR --runMode genomeGenerate --genomeDir {params.tmpdir} --genomeFastaFiles {params.tmpdir}/host.fna \
            --sjdbGTFfile {params.tmpdir}/host.gff --sjdbOverhang {params.sjdbOverhang} --runThreadN {threads} --limitGenomeGenerateRAM {params.limitGenomeGenerateRAM} \
            {params.extra_params}
        mv {params.tmpdir}/* {params.outdir}/
        rm -rf {params.tmpdir}
        """

rule star_map_host:
    """
    Maps any read that mapped to fungi against the host genome.
    """
    output:
        "results/star/{sample_id}/{sample_id}.host.bam"
    input:
        db=expand("resources/host/{f}",
            f=["Genome", "SA", "SAindex", "chrLength.txt", "chrName.txt",
               "chrNameLength.txt", "chrStart.txt", "genomeParameters.txt"]),
        R1=rules.sortmerna.output.R1,
        R2=rules.sortmerna.output.R2,
    log:
        star="results/star/{sample_id}/{sample_id}.host.log",
        all="results/star/{sample_id}/map.log",
        stat="results/star/{sample_id}/{sample_id}.Log.final.out"
    params:
        tmpdir = "$TMPDIR/{sample_id}",
        prefix = "$TMPDIR/{sample_id}/{sample_id}.",
        genomedir = lambda wildcards, input: os.path.dirname(input.db[0]),
        temp_bam = "$TMPDIR/{sample_id}/{sample_id}.Aligned.out.bam",
        temp_log = "$TMPDIR/{sample_id}/{sample_id}.Log.out",
        temp_stat = "$TMPDIR/{sample_id}/{sample_id}.Log.final.out",
        setting = config["star_extra_params"]
    conda: "../../envs/star.yaml"
    container: "docker://quay.io/biocontainers/star:2.7.11b--h5ca1c30_5"
    threads: 10
    shell:
        """
        mkdir -p {params.tmpdir}
        exec &>{log.all}
        STAR --outFileNamePrefix {params.prefix} --runThreadN {threads} \
            --genomeDir {params.genomedir} --readFilesIn {input.R1} {input.R2} \
            --readFilesCommand 'gunzip -c' --outSAMtype BAM Unsorted \
            --outSAMunmapped Within KeepPairs {params.setting}
        mv {params.temp_bam} {output[0]}
        mv {params.temp_log} {log.star}
        mv {params.temp_stat} {log.stat}
        """

rule process_host_bam:
    """
    This rule extracts reads from the host bam file that did not map to the host genome.
    """
    output:
        R1_nohost="results/filtered/{sample_id}/{sample_id}_R1.nohost.fastq.gz",
        R2_nohost="results/filtered/{sample_id}/{sample_id}_R2.nohost.fastq.gz",
        R1_host="results/filtered/{sample_id}/{sample_id}_R1.host.fastq.gz",
        R2_host="results/filtered/{sample_id}/{sample_id}_R2.host.fastq.gz"
    input:
        bam=rules.star_map_host.output[0],
    log:
        "results/star/{sample_id}/{sample_id}.process_host_bam.log"
    container: "docker://quay.io/biocontainers/samtools:1.21--h96c455f_1"
    shell:
        """
        exec &>{log}
        samtools fastq -f 4 -1 {output.R1_nohost} -2 {output.R2_nohost} -s /dev/null -0 /dev/null {input.bam}
        samtools fastq -F 4 -1 {output.R1_host} -2 {output.R2_host} -s /dev/null -0 /dev/null {input.bam}
        """

rule filter_fungal_reads:
    """
    This rule intersects the reads that mapped to fungi with the reads that did not map to the host.
    """
    output:
        R1="results/filtered/{sample_id}/{sample_id}_R1.fungi.{paired_strategy}.nohost.fastq.gz",
        R2="results/filtered/{sample_id}/{sample_id}_R2.fungi.{paired_strategy}.nohost.fastq.gz"
    input:
        R1_nohost=rules.process_host_bam.output.R1_nohost,
        R2_nohost=rules.process_host_bam.output.R2_nohost,
        R1_fungi="results/filtered/{sample_id}/{sample_id}_R1.fungi.{paired_strategy}.fastq.gz",
        R2_fungi="results/filtered/{sample_id}/{sample_id}_R2.fungi.{paired_strategy}.fastq.gz",
    shell:
        """
        seqkit grep -f <(seqkit seq -n -i {input.R1_nohost}) {input.R1_fungi} -o {output.R1}
        seqkit grep -f <(seqkit seq -n -i {input.R2_nohost}) {input.R2_fungi} -o {output.R2}
        """

rule filter_with_kraken:
    """
    This rule filters reads using Kraken2.
    1. Take reads from 'filter_fungal_reads' rule (mapped to JGI Mycocosm and not mapped to host).
    2. Remove reads classified as 'Prokaryota' by Kraken2.
    3. Add reads classified as 'Fungi' by Kraken2 that are also in the 'nohost' bin from the 'process_host_bam' rule.
    """
    output:
        R1="results/filtered/{sample_id}/{sample_id}_R1.fungi.{paired_strategy}.nohost.kraken.fastq.gz",
        R2="results/filtered/{sample_id}/{sample_id}_R2.fungi.{paired_strategy}.nohost.kraken.fastq.gz",
    input:
        R1_pre=rules.sortmerna.output.R1,
        R2_pre=rules.sortmerna.output.R2,
        #R1fungi_mapped="results/filtered/{sample_id}/{sample_id}_R1.fungi.{paired_strategy}.nohost.fastq.gz",
        #R2fungi_mapped="results/filtered/{sample_id}/{sample_id}_R2.fungi.{paired_strategy}.nohost.fastq.gz",
        R1fungi_mapped=rules.filter_fungal_reads.output.R1,
        R2fungi_mapped=rules.filter_fungal_reads.output.R2,
        R1_prok=expand("results/kraken/{kraken_db}/{{sample_id}}/taxbins/Prokaryota_R1.fastq.gz", kraken_db=config["kraken_db"]),
        R2_prok=expand("results/kraken/{kraken_db}/{{sample_id}}/taxbins/Prokaryota_R2.fastq.gz", kraken_db=config["kraken_db"]),
        R1_fungi=expand("results/kraken/{kraken_db}/{{sample_id}}/taxbins/Fungi_R1.nohost.fastq.gz", kraken_db=config["kraken_db"]),
        R2_fungi=expand("results/kraken/{kraken_db}/{{sample_id}}/taxbins/Fungi_R2.nohost.fastq.gz", kraken_db=config["kraken_db"]),
    log:
        "results/filtered/{sample_id}/{sample_id}.{paired_strategy}.filter_with_kraken.log"
    params:
        tmpdir = "$TMPDIR/{sample_id}.{paired_strategy}",
        kraken_db=config["kraken_db"]
    shell:
        """
        exec &> {log}
        echo "Filtering with Kraken2 db: {params.kraken_db}"
        mkdir -p {params.tmpdir}
        # Remove prokaryota reads from fungi mapped reads
        seqkit grep -v -f <(seqkit seq -n -i {input.R1_prok}) {input.R1fungi_mapped} -o {params.tmpdir}/R1_nohost_noprok.fastq
        seqkit grep -v -f <(seqkit seq -n -i {input.R2_prok}) {input.R2fungi_mapped} -o {params.tmpdir}/R2_nohost_noprok.fastq
        # Take the union of the prokaryota-removed reads and the fungi reads from Kraken2
        seqkit grep -f <(seqkit seq -n -i {input.R1_fungi} {params.tmpdir}/R1_nohost_noprok.fastq) {input.R1_pre} -o {output.R1}
        seqkit grep -f <(seqkit seq -n -i {input.R2_fungi} {params.tmpdir}/R2_nohost_noprok.fastq) {input.R2_pre} -o {output.R2}
        rm -rf {params.tmpdir}
        """
