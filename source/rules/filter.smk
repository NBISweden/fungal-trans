localrules:
    get_all_mapped_fungal_refs,
    filter_report,
    link_unfiltered,
    link_filtered

################
## Link fastq ##
################
rule link_unfiltered:
    input:
        "results/preprocess/{sample_id}_R{i}.cut.trim.mRNA.fastq.gz"
    output:
        "results/unfiltered/{sample_id}/{sample_id}_R{i}.fastq.gz"
    shell:
        """
        ln -s $(pwd)/{input} $(pwd)/{output}
        """

rule link_filtered:
    input:
        expand("results/{aligner}/{{sample_id}}/{{sample_id}}_R{{i}}.fungi.nohost.fastq.gz",
               aligner=config["host_aligner"])
    output:
        "results/filtered/{sample_id}/{sample_id}_R{i}.fastq.gz"
    shell:
        """
        ln -s $(pwd)/{input} $(pwd)/{output}
        """

rule fastq_stats:
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


##########################
## Bowtie2/STAR mapping ##
##########################
include: "paired_strategy.smk"

rule bowtie_build_fungi:
    """
    Builds bowtie2 index for fungal transcripts. Because gunzip may report 
    'decompression OK, trailing garbage ignored' when decompressing the transcript file
    and exiting with an error code of 2 we handle this by checking the exit code and only 
    exiting with a 1 if the exit code is exactly 1. Otherwise this rule may fail because
    Snakemake runs in strict mode.
    """
    input:
        rules.concat_transcripts.output
    output:
        expand("resources/fungi/fungi_transcripts.fasta.{index}.bt2l",
               index=range(1,5))
    log:
        bt2="resources/fungi/bowtie2.log",
        shell="resources/fungi/shell.log"
    params:
        fasta = "$TMPDIR/fungi_transcripts/fungi_transcripts.fasta",
        tmpdir = "$TMPDIR/fungi_transcripts",
        outdir = "resources/fungi/"
    threads: 64
    conda: "../../envs/bowtie2.yaml"
    shell:
        """
        set +e
        exec &>{log.shell}
        mkdir -p {params.tmpdir}
        gunzip -c {input} > {params.fasta}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 1
        fi
        bowtie2-build --threads {threads} --large-index {params.fasta} {params.fasta} >{log.bt2} 2>&1
        mv {params.fasta}*.bt2l {params.outdir}/
        rm -r {params.tmpdir}
        """

if not config["host_fna"] and config["host_aligner"]=="star":
    sys.exit("ERROR: No host genome fasta file given for STAR mapping")

def star_build_input(wildcards):
    input = []
    # if host_fna_url is given, put this download under resources/host/host.fna
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
        extra_params = config["star_extra_build_params"],
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    threads: 20
    resources:
        runtime = 60 * 24,
        mem_mb = config["star_limitGenomeGenerateRAM"]*1000, # converts ram in GB to MB
    conda: "../../envs/star.yaml"
    shell:
        """
        exec &>{log}
        mkdir -p {params.tmpdir}
        gunzip -c {input[0]} > {params.tmpdir}/host.fna
        gunzip -c {input[1]} > {params.tmpdir}/host.gff
        STAR --runMode genomeGenerate --genomeDir {params.tmpdir} --genomeFastaFiles {params.tmpdir}/host.fna \
            --sjdbGTFfile {params.tmpdir}/host.gff --runThreadN {threads} --limitGenomeGenerateRAM {params.limitGenomeGenerateRAM} \
            {params.extra_params}
        mv {params.tmpdir}/* {params.outdir}/
        rm -rf {params.tmpdir}
        """


rule bowtie_build_host:
    input:
        "resources/host/host.fna"
    output:
        expand("resources/host/host.fna.{index}.bt2l", index=range(1,5))
    log:
        "resources/host/bowtie2.log"
    resources:
        runtime = 60
    conda: "../../envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            --large-index {input} \
            {input} >{log} 2>&1
        """

rule bowtie_map_fungi:
    """
    Maps preprocessed reads against fungal transcripts.
    
    Reads that do or do not map concordantly in this step are saved directly
    in separate fastq files. If the 'paired_strategy' is set to 'concordant' 
    these files will be used to separate host and fungal reads.
    """
    input:
        R1="results/preprocess/sortmerna/{sample_id}/{sample_id}_R1.cut.trim.mRNA.fastq.gz",
        R2="results/preprocess/sortmerna/{sample_id}/{sample_id}_R2.cut.trim.mRNA.fastq.gz",
        db=expand("resources/fungi/fungi_transcripts.fasta.{index}.bt2l", index=range(1,5))
    output:
        bam="results/bowtie2/{sample_id}/{sample_id}.fungi.bam",
        R1="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.conc.fastq.gz",
        R2="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.conc.fastq.gz",
        R1u="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.noconc.fastq.gz",
        R2u="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.noconc.fastq.gz",
    params:
        prefix = "resources/fungi/fungi_transcripts.fasta",
        al_conc_path = "$TMPDIR/{sample_id}/{sample_id}_R%.fungi.conc.fastq.gz",
        un_conc_path = "$TMPDIR/{sample_id}/{sample_id}_R%.fungi.noconc.fastq.gz",
        R1 = "$TMPDIR/{sample_id}/{sample_id}_R1.fungi.conc.fastq.gz",
        R2 = "$TMPDIR/{sample_id}/{sample_id}_R2.fungi.conc.fastq.gz",
        R1u = "$TMPDIR/{sample_id}/{sample_id}_R1.fungi.noconc.fastq.gz",
        R2u = "$TMPDIR/{sample_id}/{sample_id}_R2.fungi.noconc.fastq.gz",
        tmpdir = "$TMPDIR/{sample_id}",
        temp_bam = "$TMPDIR/{sample_id}/{sample_id}.fungi.bam",
        setting = config["bowtie2_params"]
    log:
        bt2 = "results/bowtie2/{sample_id}/{sample_id}.bowtie2.fungi.log",
        st_sort = "results/bowtie2/{sample_id}/{sample_id}.samtools_sort.fungi.log",
    threads: 10
    conda: "../../envs/bowtie2.yaml"
    shell:
        """
        mkdir -p {params.tmpdir}
        bowtie2 {params.setting} -p {threads} -x {params.prefix} -1 {input.R1} \
            -2 {input.R2} --al-conc-gz {params.al_conc_path} \
            --un-conc-gz {params.un_conc_path} 2> {log.bt2} | samtools sort -n -O BAM - >{params.temp_bam} 2>{log.st_sort} 
        mv {params.temp_bam} {output.bam}
        mv {params.R1} {output.R1} 
        mv {params.R2} {output.R2}
        mv {params.R1u} {output.R1u}
        mv {params.R2u} {output.R2u}
        """

rule get_mapped_fungal_refs:
    input:
        "results/bowtie2/{sample_id}/{sample_id}.fungi.bam"
    output:
        "results/bowtie2/{sample_id}/{sample_id}.fungi.refs"
    resources:
        runtime = 30
    shell:
        """
        samtools view -F 4 -f 64 {input[0]} | cut -f3 | sort -u > {output[0]}
        """

rule get_all_mapped_fungal_refs:
    input:
        expand("results/bowtie2/{sample_id}/{sample_id}.fungi.refs",
               sample_id = samples.keys())
    output:
        "results/collated/bowtie2/fungi.refs.ids"
    shell:
        """
        cat {input} | sort -u > {output[0]}
        """

rule star_map_host:
    """
    Maps putative fungal reads against host with STAR
    """
    input:
        db=expand("resources/host/{f}",
            f=["Genome", "SA", "SAindex", "chrLength.txt", "chrName.txt",
               "chrNameLength.txt", "chrStart.txt", "genomeParameters.txt"]),
        R1="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.fastq.gz",
        R2="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.fastq.gz"
    output:
        "results/star/{sample_id}/{sample_id}.host.bam"
    log:
        star="results/star/{sample_id}/{sample_id}.host.log",
        all="results/star/{sample_id}/map.log",
        stat="results/star/{sample_id}/{sample_id}.Log.final.out"
    params:
        prefix = "$TMPDIR/{sample_id}/{sample_id}.",
        genomedir = lambda wildcards, input: os.path.dirname(input.db[0]),
        temp_bam = "$TMPDIR/{sample_id}/{sample_id}.Aligned.out.bam",
        temp_log = "$TMPDIR/{sample_id}/{sample_id}.Log.out",
        temp_stat = "$TMPDIR/{sample_id}/{sample_id}.Log.final.out",
        setting = config["star_extra_params"]
    conda: "../../envs/star.yaml"
    threads: 20
    shell:
        """
        exec &>{log.all}
        STAR --outFileNamePrefix {params.prefix} --runThreadN {threads} \
            --genomeDir {params.genomedir} --readFilesIn {input.R1} {input.R2} \
            --readFilesCommand 'gunzip -c' --outSAMtype BAM Unsorted \
            --outSAMunmapped Within KeepPairs {params.setting}
        mv {params.temp_bam} {output[0]}
        mv {params.temp_log} {log.star}
        mv {params.temp_stat} {log.stat}
        """

rule bowtie_map_host:
    """
    Maps putative fungal reads against host sequences with bowtie2 
    """
    input:
        db = expand("resources/host/host.fna.{index}.bt2l", index=range(1,5)),
        R1="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.fastq.gz",
        R2="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.fastq.gz"
    output:
        bam="results/bowtie2/{sample_id}/{sample_id}.host.bam",
        R1f="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.host.noconc.fastq.gz",
        R2f="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.host.noconc.fastq.gz",
        R1h="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.host.conc.fastq.gz",
        R2h="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.host.conc.fastq.gz"
    params:
        temp_bam = "$TMPDIR/{sample_id}/{sample_id}.host.bam",
        no_al_path = "$TMPDIR/{sample_id}/{sample_id}_R%.fungi.host.noconc.fastq.gz",
        al_path = "$TMPDIR/{sample_id}/{sample_id}_R%.fungi.host.conc.fastq.gz",
        R1h = "$TMPDIR/{sample_id}/{sample_id}_R1.fungi.host.conc.fastq.gz",
        R2h ="$TMPDIR/{sample_id}/{sample_id}_R2.fungi.host.conc.fastq.gz",
        R1f = "$TMPDIR/{sample_id}/{sample_id}_R1.fungi.host.noconc.fastq.gz",
        R2f = "$TMPDIR/{sample_id}/{sample_id}_R2.fungi.host.noconc.fastq.gz",
        tmpdir = "$TMPDIR/{sample_id}",
        prefix = "resources/host/host.fna",
        setting = config["bowtie2_params"]
    log:
        bt2 = "results/bowtie2/{sample_id}/{sample_id}.host.log",
        samtools = "results/bowtie2/{sample_id}/{sample_id}.samtools.host.log"
    threads: 10
    resources:
        runtime = 60
    shell:
        """
        exec &> {log.samtools}
        mkdir -p {params.tmpdir}
        bowtie2 {params.setting} -p {threads} -x {params.prefix} -1 {input.R1} \
            -2 {input.R2} --al-conc-gz {params.al_path} --un-conc-gz {params.no_al_path} \
            2>{log.bt2} | samtools sort -n -O BAM - > {params.temp_bam} 2>{log.samtools}
        mv {params.temp_bam} {output.bam}
        mv {params.R1f} {output.R1f}
        mv {params.R2f} {output.R2f}
        mv {params.R1h} {output.R1h}
        mv {params.R2h} {output.R2h}
        """

rule host_reads:
    """
    This rule extracts reads from the preprocessed fastq files that are either not mapped to fungi,
    or are marked as putative host reads.
    """
    input:
        R1="results/preprocess/sortmerna/{sample_id}/{sample_id}_R1.cut.trim.mRNA.fastq.gz",
        R2="results/preprocess/sortmerna/{sample_id}/{sample_id}_R2.cut.trim.mRNA.fastq.gz",
        R1_1 = "results/bowtie2/{sample_id}/{sample_id}_R1.nonfungi.fastq.gz",
        R2_1 = "results/bowtie2/{sample_id}/{sample_id}_R2.nonfungi.fastq.gz",
        R1_2 = "results/"+config["host_aligner"]+"/{sample_id}/{sample_id}_R1.fungi.putative-host.fastq.gz",
        R2_2 = "results/"+config["host_aligner"]+"/{sample_id}/{sample_id}_R2.fungi.putative-host.fastq.gz"
    output:
        R1 = "results/host/{sample_id}_R1.host.fastq.gz",
        R2 = "results/host/{sample_id}_R2.host.fastq.gz"
    params:
        tmpids="$TMPDIR/{sample_id}.tmp",
        ids = "$TMPDIR/{sample_id}.ids",
        R1 = "$TMPDIR/{sample_id}_R1.host.fastq.gz",
        R2 = "$TMPDIR/{sample_id}_R2.host.fastq.gz"
    log:
        "results/host/{sample_id}.log"
    threads: 1
    shell:
        """
        exec &>{log}
        touch {params.tmpids}
        seqkit seq -n -i --threads {threads} {input.R1_1} {input.R1_2} {input.R2_1} {input.R2_2} | sort -u > {params.ids}
        seqkit grep -f {params.ids} --immediate-output --threads {threads} -o {params.R1} {input.R1}
        seqkit grep -f {params.ids} --immediate-output --threads {threads} -o {params.R2} {input.R2}
        mv {params.R1} {output.R1}
        mv {params.R2} {output.R2}
        rm {params.ids}
        rm {params.tmpids}
        """

def get_host_logs(config, samples):
    if config["host_aligner"] == "star":
        return expand("results/star/{sample_id}/{sample_id}.host.log", sample_id = samples.keys())
    elif config["host_aligner"] == "bowtie2":
        return expand("results/bowtie2/{sample_id}/{sample_id}.host.log", sample_id = samples.keys())

rule filter_report:
    input:
        hostlogs = get_host_logs(config, samples),
        fungallogs = expand("results/bowtie2/{sample_id}/{sample_id}.bowtie2.fungi.log",
            sample_id = samples.keys())
    output:
        "results/report/filtering/filter_report.html"
    log:
        "results/report/filtering/filter_report.log"
    params:
        tmpdir = "multiqc_filter",
        config = "config/multiqc_filter_config.yaml"
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.hostlogs} {params.tmpdir}
        cp {input.fungallogs} {params.tmpdir}
        multiqc \
            -f -c {params.config} \
            -n filter_report \
            -o results/report/filtering {params.tmpdir}
        rm -r {params.tmpdir}
        """