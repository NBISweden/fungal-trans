localrules:
    strobealign_collect_chunks,
    link_unfiltered,
    link_filtered,
    fastq_stats,

################
## Link fastq ##
################
rule link_unfiltered:
    """
    This rule creates symlinks to the unfiltered fastq files.
    """
    input:
        "results/preprocess/{sample_id}_R{i}.mRNA.fastq.gz"
    output:
        "results/unfiltered/{sample_id}/{sample_id}_R{i}.fastq.gz"
    shell:
        """
        ln -s $(pwd)/{input} $(pwd)/{output}
        """

rule link_filtered:
    """
    This rule creates symlinks to the filtered fastq files.
    """
    input:
        expand("results/star/{{sample_id}}/{paired_strategy}/{{sample_id}}_R{{i}}.fungi.nohost.fastq.gz", paired_strategy=config["paired_strategy"])
    output:
        "results/filtered/{sample_id}/{sample_id}_R{i}.fastq.gz"
    shell:
        """
        ln -s $(pwd)/{input} $(pwd)/{output}
        """

rule fastq_stats:
    """
    This rule calculates statistics for the fastq files.
    """
    input:
        R1="results/{source}/{sample_id}/{sample_id}_R1.fastq.gz",
        R2="results/{source}/{sample_id}/{sample_id}_R2.fastq.gz"
    output:
        "results/{source}/{sample_id}/{sample_id}.stats.tsv"
    wildcard_constraints:
        source="unfiltered|filtered"
    shell:
        """
        seqkit stats -T {input.R1} {input.R2} > {output}
        """

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
    input:
        R1=rules.sortmerna.output.R1,
        R2=rules.sortmerna.output.R2,
        fna=lambda wildcards: strobealign_chunks[wildcards.chunk]
    output:
        paf="results/strobealign/{sample_id}/{chunk}/{sample_id}.fungi.paf.gz",
        both="results/strobealign/{sample_id}/{chunk}/{sample_id}.fungi.both_mapped.gz",
    log:
        "results/strobealign/{sample_id}/{chunk}/{sample_id}.strobealign.log"
    params:
        tmpdir = "$TMPDIR/{sample_id}.{chunk}",
    container: "docker://quay.io/biocontainers/strobealign:0.15.0--h5ca1c30_1"
    threads: 6
    shell:
        """
        mkdir -p {params.tmpdir}
        cat {input.fna} > {params.tmpdir}/ref.fna.gz
        strobealign -v --mcs -x -t {threads} {params.tmpdir}/ref.fna.gz {input.R1} {input.R2} 2>{log} | igzip > {output.paf}
        igzip -c -d {output.paf} | cut -f1 | sort | uniq -d | igzip > {output.both}
        rm -r {params.tmpdir}
        """

rule strobealign_collect_chunks:
    """
    This rule collects the PAF output from the different chunks and concatenates them.
    """
    input:
        paf=expand("results/strobealign/{{sample_id}}/{chunk}/{{sample_id}}.fungi.paf.gz", chunk=strobealign_chunks.keys()),
        both=expand("results/strobealign/{{sample_id}}/{chunk}/{{sample_id}}.fungi.both_mapped.gz", chunk=strobealign_chunks.keys())
    output:
        paf="results/strobealign/{sample_id}/{sample_id}.fungi.paf.gz",
        both="results/strobealign/{sample_id}/{sample_id}.fungi.both_mapped.gz",
    shell:
        """
        cat {input.paf} > {output.paf}
        cat {input.both} > {output.both}
        """

rule process_fungal_bam:
    """
    This rule extracts reads from the preprocessed fastq files that map to fungi.
    Reads with both ends and one or both ends mapped are extracted separately.
    """
    input:
        R1=rules.sortmerna.output.R1,
        R2=rules.sortmerna.output.R2,
        paf="results/strobealign/{sample_id}/{sample_id}.fungi.paf.gz",
        both="results/strobealign/{sample_id}/{sample_id}.fungi.both_mapped.gz",
    output:
        R1one="results/strobealign/{sample_id}/{sample_id}_R1.fungi.one_mapped.fastq.gz",
        R2one="results/strobealign/{sample_id}/{sample_id}_R2.fungi.one_mapped.fastq.gz",
        R1both="results/strobealign/{sample_id}/{sample_id}_R1.fungi.both_mapped.fastq.gz",
        R2both="results/strobealign/{sample_id}/{sample_id}_R2.fungi.both_mapped.fastq.gz",
    shell:
        """
        # Output reads that map with one or two ends to fungi
        seqkit grep -f <(gunzip -c {input.paf} | cut -f1) {input.R1} -o {output.R1one}
        seqkit grep -f <(gunzip -c {input.paf} | cut -f1) {input.R2} -o {output.R2one}
        # Output reads that map with both ends to fungi
        seqkit grep -f <(gunzip -c {input.both}) {input.R1} -o {output.R1both}
        seqkit grep -f <(gunzip -c {input.both}) {input.R2} -o {output.R2both}
        """

if not config["host_fna"]:
    sys.exit("ERROR: No host genome fasta file given for STAR mapping")

def star_build_input(wildcards):
    input = []
    if config["host_fna"]:
        input.append(config["host_fna"])
    if config["host_gff"]:
        input.append(config["host_gff"])
    return input

rule star_build_host:
    input:
        star_build_input,
    output:
        expand("resources/host/{f}",
            f = ["Genome", "SA", "SAindex", "chrLength.txt", "chrName.txt",
                 "chrNameLength.txt", "chrStart.txt", "genomeParameters.txt"])
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
        gunzip -c {input[0]} > {params.tmpdir}/host.fna
        gunzip -c {input[1]} > {params.tmpdir}/host.gff
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
    input:
        db=expand("resources/host/{f}",
            f=["Genome", "SA", "SAindex", "chrLength.txt", "chrName.txt",
               "chrNameLength.txt", "chrStart.txt", "genomeParameters.txt"]),
        R1=rules.sortmerna.output.R1,
        R2=rules.sortmerna.output.R2,
    output:
        "results/star/{sample_id}/{sample_id}.host.bam"
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
    input:
        bam=rules.star_map_host.output[0],
    output:
        R1_nohost="results/star/{sample_id}/{sample_id}_R1.nohost.fastq.gz",
        R2_nohost="results/star/{sample_id}/{sample_id}_R2.nohost.fastq.gz",
        R1_host="results/star/{sample_id}/{sample_id}_R1.host.fastq.gz",
        R2_host="results/star/{sample_id}/{sample_id}_R2.host.fastq.gz"
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
    input:
        R1_nohost=rules.process_host_bam.output.R1_nohost,
        R2_nohost=rules.process_host_bam.output.R2_nohost,
        R1_fungi="results/strobealign/{sample_id}/{sample_id}_R1.fungi.{paired_strategy}.fastq.gz",
        R2_fungi="results/strobealign/{sample_id}/{sample_id}_R2.fungi.{paired_strategy}.fastq.gz",
    output:
        R1="results/filtered/{sample_id}/{paired_strategy}/{sample_id}_R1.fungi.nohost.fastq.gz",
        R2="results/filtered/{sample_id}/{paired_strategy}/{sample_id}_R2.fungi.nohost.fastq.gz"
    shell:
        """
        seqkit grep -f <(seqkit seq -n -i {input.R1_nohost}) {input.R1_fungi} -o {output.R1}
        seqkit grep -f <(seqkit seq -n -i {input.R2_nohost}) {input.R2_fungi} -o {output.R2}
        """

rule remove_prok_add_fungi:
    """
    This rule removes reads classified by Kraken2 as Prokaryota and adds reads classified as Fungi.
    """
    input:
        R1prok="results/kraken/{kraken_db}/{sample_id}/taxbins/Prokaryota_R1.fastq.gz",
        R2prok="results/kraken/{kraken_db}/{sample_id}/taxbins/Prokaryota_R2.fastq.gz",
        R1fungi="results/kraken/{kraken_db}/{sample_id}/taxbins/Fungi_R1.fastq.gz",
        R2fungi="results/kraken/{kraken_db}/{sample_id}/taxbins/Fungi_R2.fastq.gz",
        R1_nohost="results/star/{sample_id}/{sample_id}_R1.nohost.fastq.gz",
        R2_nohost="results/star/{sample_id}/{sample_id}_R2.nohost.fastq.gz",
        R1one="results/strobealign/{sample_id}/{sample_id}_R1.fungi.one_mapped.fastq.gz",
        R2one="results/strobealign/{sample_id}/{sample_id}_R2.fungi.one_mapped.fastq.gz",
        R1both="results/strobealign/{sample_id}/{sample_id}_R1.fungi.both_mapped.fastq.gz",
        R2both="results/strobealign/{sample_id}/{sample_id}_R2.fungi.both_mapped.fastq.gz",
    shell:
        """
        # Reads must be in R*_nohost
        seqkit grep -f <(seqkit seq -n -i {input.R1_nohost}) {input.R1fungi} -o {params.tmpdir}/kraken_fungi
        """