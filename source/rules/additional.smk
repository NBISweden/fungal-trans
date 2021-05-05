localrules: stat_additional

rule map_additional:
    """
    This rule maps preprocessed and filtered reads against a database
    to investigate alignment rate.
    """
    input:
        R1 = "results/filtered/{sample}/{sample}_R1.filtered.union.fastq.gz",
        R2 = "results/filtered/{sample}/{sample}_R2.filtered.union.fastq.gz",
        db = expand("resources/{db}/{db}.fasta.{i}.bt2l",
                    db=config["additional_bt_db"],
                    i=["1","2","3","4","rev.1","rev.2"])
    output:
        temp("results/map/additional/{db}/{sample}/{sample}.bam")
    log:
        "results/map/additional/{db}/{sample}/{sample}.bowtie2.log"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    threads: 10
    params:
        index = "resources/{db}/{db}.fasta",
        temp_bam = "$TMPDIR/{db}.{sample}.bam"
    shell:
        """
        bowtie2 -x {params.index} -1 {input.R1} -2 {input.R2} --very-sensitive \
            -p {threads} 2>{log} | samtools view -b -h - 2>/dev/null | \
            samtools sort -o {params.temp_bam} -O BAM
        mv {params.temp_bam} {output}
        """

rule stat_additional:
    """
    Filters bam file with quality and mapping criteria and runs flagstat
    """
    input:
        "results/map/additional/{db}/{sample}/{sample}.bam"
    output:
        "results/map/additional/{db}/{sample}/{sample}.bam.flagstat"
    shell:
        """
        samtools view -F 12 -h -b -q 10 {input} | samtools \
            flagstat - > {output}
        """


rule multiqc_additional:
    """
    Runs multiqc on the bowtie2 and featurecounts logs
    """
    input:
        fl = expand("results/map/additional/{db}/{sample}/{sample}.bam.flagstat",
                    db=config["additional_bt_db"],
                    sample=samples.keys()),
        bt = expand("results/map/additional/{db}/{sample}/{sample}.bowtie2.log",
                    db=config["additional_bt_db"],
                    sample=samples.keys())
    output:
        html = "results/report/additional/{db}/assign_report.html",
        fc = "results/report/additional/{db}/assign_report_data/multiqc_samtools_flagstat.txt",
        bt = "results/report/additional/{db}/assign_report_data/multiqc_bowtie2.txt"
    log:
        "results/report/additional/{db}/assign_report.log"
    params:
        outdir = lambda w, output: os.path.dirname(output.html),
        multiqc_config = "config/multiqc_config.yaml"
    shadow: True
    shell:
        """
        mkdir {params.outdir}/input
        cp {input} {params.outdir}/input
        multiqc -f -o {params.outdir} -c {params.multiqc_config} \
            -n assign_report.html {params.outdir}/input 
        rm -r {params.outdir}/input
        """
