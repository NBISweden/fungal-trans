localrules:
    multiqc,
    download_fastq,
    download_rRNA_database

#############
## SAMPLES ##
#############
def generate_sra_acc(wildcards):
    try:
        acc = samples[wildcards.sample_id]['accession']
    except KeyError:
        acc = wildcards.sample_id
    return acc

rule download_fastq:
    """Downloads fastq files from a sequence read archive.
    """
    output:
        R1 = os.path.join(config["datadir"], "{sample_id}_1.fastq.gz"),
        R2 = os.path.join(config["datadir"], "{sample_id}_2.fastq.gz")
    log:
        os.path.join(config["datadir"], "{sample_id}.log")
    params:
        tmpdir = "$TMPDIR/{sample_id}",
        outdir = config["datadir"],
        dlparams = config["dlparams"],
        acc=lambda wildcards: generate_sra_acc(wildcards)
    conda: "../../envs/sratools.yaml"
    container: "docker://quay.io/biocontainers/sra-tools:3.1.1--h4304569_2"
    shell:
        """
        mkdir -p {params.tmpdir}
        fastq-dump {params.dlparams} \
            --split-files {params.acc} -O {params.tmpdir} > {log} 2>{log}
        pigz {params.tmpdir}/*_1.fastq > {output.R1}
        pigz {params.tmpdir}/*_2.fastq > {output.R2}
        mv {params.tmpdir}/*_1.fastq.gz {output.R1}
        mv {params.tmpdir}/*_2.fastq.gz {output.R2}
        rm -rf {params.tmpdir}
        """

##############
## TRIMMING ##
##############

rule cutadapt:
    input:
        R1 = lambda wildcards: samples[wildcards.sample_id]["R1"],
        R2 = lambda wildcards: samples[wildcards.sample_id]["R2"]
    output:
        R1 = "results/preprocess/cutadapt/{sample_id}_R1.cut.fastq.gz",
        R2 = "results/preprocess/cutadapt/{sample_id}_R2.cut.fastq.gz",
    log:
        "results/preprocess/{sample_id}.cutadapt.log"
    resources:
        runtime = lambda wildcards, attempt: attempt*20
    params:
        a = config["cutadapt_R1"],
        A = config["cutadapt_R2"],
        cutadapt_error_rate = config["cutadapt_error_rate"],
        extra_settings = config["cutadapt_extra_settings"],
        tmpdir = "$TMPDIR/{sample_id}"
    conda: "../../envs/cutadapt.yaml"
    container: "docker://quay.io/biocontainers/cutadapt:5.0--py310h1fe012e_0"
    threads: 4
    shell:
        """
        mkdir -p {params.tmpdir}
        cutadapt {params.extra_settings} -e {params.cutadapt_error_rate} -j {threads} -a {params.a} -A {params.A} \
        -o {params.tmpdir}/R1.fastq.gz -p {params.tmpdir}/R2.fastq.gz {input.R1} {input.R2} > {log} 2>{log}
        mv {params.tmpdir}/R1.fastq.gz {output.R1}
        mv {params.tmpdir}/R2.fastq.gz {output.R2}
        rm -r {params.tmpdir}
        """

rule trimmomatic:
    input:
        R1 = "results/preprocess/cutadapt/{sample_id}_R1.cut.fastq.gz",
        R2 = "results/preprocess/cutadapt/{sample_id}_R2.cut.fastq.gz"
    output:
        R1 = "results/preprocess/trimmomatic/{sample_id}_R1.cut.trim.fastq.gz",
        R2 = "results/preprocess/trimmomatic/{sample_id}_R2.cut.trim.fastq.gz",
    log:
        "results/preprocess/trimmomatic/{sample_id}.cut.trim.log"
    params:
        jarpath=config["trimmomatic_home"]+"/trimmomatic.jar",
        qc_trim_settings=config["qc_trim_settings"],
        tmpdir = "$TMPDIR/{sample_id}",
        outdir = "results/preprocess",
        R1_temp = "$TMPDIR/{sample_id}/R1.fastq.gz",
        R2_temp = "$TMPDIR/{sample_id}/R2.fastq.gz"
    threads: 4
    conda: "../../envs/trimmomatic.yaml"
    container: "docker://biocontainers/trimmomatic:v0.38dfsg-1-deb_cv1"
    resources: 
        runtime = lambda wildcards, attempt: attempt*20
    shell:
        """
        mkdir -p {params.tmpdir}
        trimmomatic PE -threads {threads} {input.R1} {input.R2} {params.R1_temp} /dev/null \
        {params.R2_temp} /dev/null {params.qc_trim_settings} 2>{log}
        mv {params.R1_temp} {output.R1}
        mv {params.R2_temp} {output.R2}
        rm -r {params.tmpdir}
        """

####################
## rRNA filtering ##
####################

rule download_rRNA_database:
    output:
        default="resources/sortmerna/rRNA_databases_v4/smr_v4.3_default_db.fasta", 
        fast="resources/sortmerna/rRNA_databases_v4/smr_v4.3_fast_db.fasta",
        sensitive="resources/sortmerna/rRNA_databases_v4/smr_v4.3_sensitive_db.fasta",
        sensitive_rfam="resources/sortmerna/rRNA_databases_v4/smr_v4.3_sensitive_db_rfam_seeds.fasta"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        wget -O {params.outdir}/database.tar.gz https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
        tar -xvf {params.outdir}/database.tar.gz -C {params.outdir}
        rm {params.outdir}/database.tar.gz
        """

rule sortmerna_index:
    output:
        touch("resources/sortmerna/rRNA_databases_v4/indexed")
    input:
        rules.download_rRNA_database.output.default
    log:
        "resources/sortmerna/rRNA_databases_v4/index.log"
    params:
        tmpdir="$TMPDIR/sortmerna_index",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    resources:
        runtime = 60
    threads: 2
    conda: "../../envs/sortmerna.yaml"
    container: "docker://quay.io/biocontainers/sortmerna:4.2.0--h9ee0642_1"
    shell:
        """
        mkdir -p {params.tmpdir}
        head -2 {input} > {params.tmpdir}/tmp.fasta
        sortmerna --workdir {params.tmpdir}/ --threads {threads} \
            --idx {params.outdir} --reads {params.tmpdir}/tmp.fasta \
            --ref {input[0]} --task 1 --blast 1 --num_alignments 1 > {log} 2>&1
        rm -r {params.tmpdir}
        """

rule sortmerna:
    input:
        R1=rules.trimmomatic.output.R1,
        R2=rules.trimmomatic.output.R2,
        ref=rules.download_rRNA_database.output.default
    output:
        R1="results/preprocess/sortmerna/{sample_id}/{sample_id}_R1.cut.trim.mRNA.fastq.gz",
        R2="results/preprocess/sortmerna/{sample_id}/{sample_id}_R2.cut.trim.mRNA.fastq.gz",
        R1_rRNA="results/preprocess/sortmerna/{sample_id}/{sample_id}_R1.cut.trim.rRNA.fastq.gz",
        R2_rRNA="results/preprocess/sortmerna/{sample_id}/{sample_id}_R2.cut.trim.rRNA.fastq.gz"
    log:
        run="results/preprocess/sortmerna/{sample_id}/run.log",
        stats="results/preprocess/sortmerna/{sample_id}/{sample_id}.aligned.log",
    params:
        workdir="$TMPDIR/{sample_id}.sortmerna",
        idx_dir=lambda wildcards, input: os.path.dirname(input.ref)
    threads: 4
    conda:
        "../../envs/sortmerna.yaml"
    container: "docker://quay.io/biocontainers/sortmerna:4.2.0--h9ee0642_1"
    shell:
        """
        rm -rf {params.workdir}
        sortmerna --workdir {params.workdir}/ --threads {threads} --idx {params.idx_dir} \
            --ref {input.ref} --reads {input.R1} --reads {input.R2} \
            --paired_in --out2 --fastx --blast 1 --num_alignments 1 --aligned --other -v >{log.run} 2>&1
        gzip {params.workdir}/out/*.fastq
        mv {params.workdir}/out/aligned_fwd.fastq.gz {output.R1_rRNA}
        mv {params.workdir}/out/aligned_rev.fastq.gz {output.R2_rRNA}
        mv {params.workdir}/out/other_fwd.fastq.gz {output.R1}
        mv {params.workdir}/out/other_rev.fastq.gz {output.R2}
        mv {params.workdir}/out/aligned.log {log.stats}
        """

rule fastqc:
    input:
        R1 = "results/preprocess/sortmerna/{sample_id}/{sample_id}_R1.cut.trim.mRNA.fastq.gz",
        R2 = "results/preprocess/sortmerna/{sample_id}/{sample_id}_R2.cut.trim.mRNA.fastq.gz"
    output:
        R1 = "results/preprocess/fastqc/{sample_id}_R1.cut.trim.mRNA_fastqc.zip",
        R2 = "results/preprocess/fastqc/{sample_id}_R2.cut.trim.mRNA_fastqc.zip"
    log: "results/preprocess/fastqc/{sample_id}.fastqc.log"
    params:
        outdir = "results/preprocess/fastqc"
    resources:
        runtime = lambda wildcards, attempt: attempt*20
    conda: "../../envs/fastqc.yaml"
    container: "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    shell:
        """
        fastqc --noextract {input.R1} {input.R2} -o {params.outdir} >{log} 2>&1
        """

rule multiqc:
    input:
        qclogs = expand("results/preprocess/fastqc/{sample_id}_R{i}.cut.trim.mRNA_fastqc.zip", sample_id = samples.keys(), i = [1,2]),
        cutlogs = expand("results/preprocess/cutadapt/{sample_id}.cutadapt.log", sample_id = samples.keys()),
        trimlogs = expand("results/preprocess/trimmomatic/{sample_id}.cut.trim.log", sample_id = samples.keys()),
        sortmernalogs = expand("results/preprocess/sortmerna/{sample_id}/{sample_id}.aligned.log", sample_id = samples.keys())
    output:
        "results/report/preprocess/preprocess_report.html",
        "results/report/preprocess/preprocess_report_data/multiqc.log"
    log:
        "results/report/preprocess/preprocess_report.log"
    params:
        outdir = "results/report/preprocess",
        config = workflow.source_path("../../config/multiqc_config.yaml")
    shell:
        """
        multiqc -f -o {params.outdir} -c {params.config} -n preprocess_report {input.cutlogs} \
            {input.trimlogs} {input.sortmernalogs} {input.qclogs} >{log} 2>&1
        """