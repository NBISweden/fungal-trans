localrules:
    multiqc,
    link_files,
    download_fastq,
    avg_seq_length_bowtie,
    avg_seq_length_filtered,
    avg_seq_length_taxmapper

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
    shell:
        """
        mkdir -p {params.tmpdir}
        fastq-dump {params.dlparams} \
            --split-files {params.acc} -O {params.tmpdir} > {log} 2>{log}
        gzip -c {params.tmpdir}/*_1.fastq > {output.R1}
        gzip -c {params.tmpdir}/*_2.fastq > {output.R2}        
        rm -rf {params.tmpdir}
        """

##################
## READ LENGTHS ##
##################
def calc_read_lengths(input, output, replace):
    import numpy as np
    sample_lengths = {}
    lengths = {}
    for f in input:
        basename = os.path.basename(f)
        sample_id = basename.replace(replace,"")
        for line in shell("seqtk seq -f 0.01 {f} | seqtk comp | cut -f2 | sort | uniq -c", iterable = True):
            line = (line.rstrip()).lstrip()
            items = line.split(" ")
            l = [int(items[1])]*int(items[0])
            try:
                lengths[sample_id] += l
            except KeyError:
                lengths[sample_id] = l
        sample_lengths[sample_id] = np.round(np.mean(lengths[sample_id]),2)
    df = pd.DataFrame(sample_lengths,index=["avg_len"]).T
    df.to_csv(output[0], sep="\t")

rule avg_seq_length_taxmapper:
    input:
        "results/taxmapper/{sample_id}/{sample_id}_R1.cut.trim.filtered.fastq.gz"
    output:
        temp("results/sample_info/{sample_id}.taxmapper_read_lengths.tab")
    run:
        calc_read_lengths(input, output, "_R1.cut.trim.filtered.fastq.gz")

rule avg_seq_length_bowtie:
    input:
        "results/bowtie2/{sample_id}/{sample_id}_R1.fungi.nospruce.fastq.gz"
    output:
        temp("results/sample_info/{sample_id}.bowtie2_read_lengths.tab")
    run:
        calc_read_lengths(input, output, "_R1.fungi.nospruce.fastq.gz")

rule avg_seq_length_filtered:
    input:
        "results/filtered/{sample_id}/{sample_id}_R1.filtered.union.fastq.gz"
    output:
        temp("results/sample_info/{sample_id}.filtered_read_lengths.tab")
    run:
        calc_read_lengths(input, output, "_R1.filtered.union.fastq.gz")

rule link_files:
    input:
        lambda wildcards: samples[wildcards.sample_id][wildcards.pair]
    output:
        "results/intermediate/{sample_id}_{pair}.fastq.gz"
    shell:
        """
        ln -s $(pwd)/{input} $(pwd)/{output}
        """

rule cutadapt:
    input:
        R1 = "results/intermediate/{sample_id}_R1.fastq.gz",
        R2 = "results/intermediate/{sample_id}_R2.fastq.gz"
    output:
        R1 = "results/preprocess/{sample_id}_R1.cut.fastq.gz",
        R2 = "results/preprocess/{sample_id}_R2.cut.fastq.gz",
    log:
        "results/preprocess/{sample_id}.cutadapt.log"
    resources:
        runtime = lambda wildcards, attempt: attempt*20
    params:
        a = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        A = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        settings = "-e 0.2",
        tmpdir = "$TMPDIR/{sample_id}"
    conda: "../../envs/cutadapt.yaml"
    threads: 4
    shell:
        """
        mkdir -p {params.tmpdir}
        cutadapt {params.settings} -j {threads} -a {params.a} -A {params.A} \
        -o {params.tmpdir}/R1.fastq.gz -p {params.tmpdir}/R2.fastq.gz {input.R1} {input.R2} > {log} 2>{log}
        mv {params.tmpdir}/R1.fastq.gz {output.R1}
        mv {params.tmpdir}/R2.fastq.gz {output.R2}
        rm -r {params.tmpdir}
        """

rule trimmomatic:
    input:
        R1 = "results/preprocess/{sample_id}_R1.cut.fastq.gz",
        R2 = "results/preprocess/{sample_id}_R2.cut.fastq.gz"
    output:
        R1 = "results/preprocess/{sample_id}_R1.cut.trim.fastq.gz",
        R2 = "results/preprocess/{sample_id}_R2.cut.trim.fastq.gz",
    log:
        "results/preprocess/{sample_id}.cut.trim.log"
    params:
        jarpath=config["trimmomatic_home"]+"/trimmomatic.jar",
        qc_trim_settings=config["qc_trim_settings"],
        tmpdir = "$TMPDIR/{sample_id}",
        outdir = "results/preprocess",
        R1_temp = "$TMPDIR/{sample_id}/R1.fastq.gz",
        R2_temp = "$TMPDIR/{sample_id}/R2.fastq.gz"
    threads: 4
    conda: "../../envs/trimmomatic.yaml"
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

rule fastqc:
    input:
        R1 = "results/preprocess/{sample_id}_R1.cut.trim.fastq.gz",
        R2 = "results/preprocess/{sample_id}_R2.cut.trim.fastq.gz"
    output:
        R1 = "results/preprocess/{sample_id}_R1.cut.trim_fastqc.zip",
        R2 = "results/preprocess/{sample_id}_R2.cut.trim_fastqc.zip"
    log: "results/preprocess/{sample_id}.fastqc.log"
    params:
        outdir = "results/preprocess"
    resources:
        runtime = lambda wildcards, attempt: attempt*20
    conda: "../../envs/fastqc.yaml"
    shell:
        """
        fastqc --noextract {input.R1} {input.R2} -o {params.outdir} >{log} 2>&1
        """

def reformat_cutadapt_log(f, out1, out2):
    sample = (os.path.basename(f)).replace(".cutadapt.log","")
    with open(f, 'r') as fhin, open(out1, 'w') as fhout1, open(out2, 'w') as fhout2:
        for line in fhin:
            line = line.rstrip()
            if "Command line parameters:" in line:
                items = line.rsplit()
                fhout1.write(" ".join(items[0:-1]+[sample+"_R1"])+"\n")
                fhout2.write(" ".join(items[0:-1]+[sample+"_R2"])+"\n")
            else:
                fhout1.write("{}\n".format(line))
                fhout2.write("{}\n".format(line))

def reformat_trimmomatic_log(f, out1, out2):
    sample = (os.path.basename(f)).replace(".cut.trim.log","")
    with open(f, 'r') as fhin, open(out1, 'w') as fhout1, open(out2, 'w') as fhout2:
        for line in fhin:
            line = line.rstrip()
            fhout1.write("{}\n".format(line.replace("_R2.cut.fastq.gz","_R1.cut.fastq.gz")))
            fhout2.write("{}\n".format(line.replace("_R1.cut.fastq.gz","_R2.cut.fastq.gz")))

rule multiqc:
    input:
        qclogs = expand("results/preprocess/{sample_id}_R{i}.cut.trim_fastqc.zip", sample_id = samples.keys(), i = [1,2]),
        cutlogs = expand("results/preprocess/{sample_id}.cutadapt.log", sample_id = samples.keys()),
        trimlogs = expand("results/preprocess/{sample_id}.cut.trim.log", sample_id = samples.keys())
    output:
        "results/report/preprocess/preprocess_report.html",
        "results/report/preprocess/preprocess_report_data/multiqc.log"
    params:
        dir = "results/preprocess",
        outdir = "results/report/preprocess",
        tmpdir = os.path.join(os.path.expandvars("$TMPDIR"),"multiqc"),
        multiqc_config = "config/multiqc_config.yaml"
    run:
        shell("mkdir -p {params.tmpdir}")
        shell("cp {input} {params.tmpdir}")
        for f in input.cutlogs:
            sample = (os.path.basename(f)).replace(".cutadapt.log","")
            out1 = "{dir}/{sample}_R1.cut.log".format(dir=params.tmpdir, sample=sample)
            out2 = "{dir}/{sample}_R2.cut.log".format(dir=params.tmpdir, sample=sample)
            reformat_cutadapt_log(f, out1, out2)
            shell("rm {params.tmpdir}/{sample}.cutadapt.log")
        for f in input.trimlogs:
            sample = (os.path.basename(f)).replace(".cut.trim.log","")
            out1 = "{dir}/{sample}_R1.cut.trim.log".format(dir=params.tmpdir, sample=sample)
            out2 = "{dir}/{sample}_R2.cut.trim.log".format(dir=params.tmpdir, sample=sample)
            reformat_trimmomatic_log(f, out1, out2)
        shell("multiqc -f -o {params.tmpdir} -c {params.multiqc_config} -n preprocess_report.html {params.tmpdir}")
        shell("rsync -azv {params.tmpdir}/preprocess_report* {params.outdir}")
        shell("rm -r {params.tmpdir}")
