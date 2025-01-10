localrules: multiqc_map_report, multiqc_map_report_co, gff2bed, infer_experiment

###############################
## Mapping for co-assemblies ##
###############################
rule bowtie_build_co:
    input:
        "results/co-assembly/{assembler}/{assembly}/final.fa"
    output:
        expand("results/map/co-assembly/{{assembler}}/{{assembly}}/final.fa.{index}.bt2l", index=range(1,5))
    params:
        prefix = "results/map/co-assembly/{assembler}/{assembly}/final.fa"
    threads: 1
    resources:
        runtime = 60*4
    shell:
        """
        bowtie2-build --large-index --threads {threads} {input[0]} {params.prefix}
        """

rule bowtie_co:
    input:
        bt = expand("results/map/co-assembly/{{assembler}}/{{assembly}}/final.fa.{index}.bt2l", index=range(1,5)),
        R1 = lambda wildcards: map_dict[wildcards.sample_id]["R1"],
        R2 = lambda wildcards: map_dict[wildcards.sample_id]["R2"]
    output:
        bam = "results/map/co-assembly/{assembler}/{assembly}/{sample_id}.bam"
    log: "results/map/co-assembly/{assembler}/{assembly}/{sample_id}.log"
    threads: 10
    resources:
        runtime = 60*4
    params:
        prefix = "results/map/co-assembly/{assembler}/{assembly}/final.fa",
        tmp_out = "$TMPDIR/{assembly}/{sample_id}.bam",
        tmp_dir = "$TMPDIR/{assembly}",
        setting = config["bowtie2_params"]
    shell:
        """
        mkdir -p {params.tmp_dir}
        bowtie2 {params.setting} -p {threads} -x {params.prefix} \
         -1 {input.R1} -2 {input.R2} 2>{log} | samtools view -h -b - | samtools sort - -o {params.tmp_out}
         mv {params.tmp_out} {output.bam}
        """

rule multiqc_map_report_co:
    input:
        bt_logs = expand("results/map/co-assembly/{{assembler}}/{{assembly}}/{sample_id}.log",
            sample_id = samples.keys()),
        fc_logs = expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/featureCounts/{sample_id}.fc.tab.summary",
            sample_id = samples.keys())
    output:
        "results/report/map/{assembler}.{assembly}_map_report.html"
    params:
        name = "{assembler}.{assembly}_map_report.html",
        tmpdir = "$TMPDIR/{assembly}"
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.bt_logs} {params.tmpdir}
        cp {input.fc_logs} {params.tmpdir}
        multiqc -o results/report/map/ -n {params.name} {params.tmpdir}
        rm -r {params.tmpdir}
        """

#######################################
## Mapping for individual assemblies ##
#######################################
rule bowtie_build:
    input:
        "results/assembly/{assembler}/{source}/{sample_id}/final.fa"
    output:
        expand("results/map/{{assembler}}/{{source}}/{{sample_id}}/final.fa.{index}.bt2l", index=range(1,5))
    params:
        prefix = "results/map/{assembler}/{source}/{sample_id}/final.fa"
    container: "docker://quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5"
    threads: 1
    resources:
        runtime = 60*24
    shell:
        """
        bowtie2-build --large-index --threads {threads} {input[0]} {params.prefix}
        """

rule bowtie:
    input:
        #R1 = "results/assembly/{assembler}/{source}/{sample_id}/{sample_id}_R1.fastq.gz",
        #R2 = "results/assembly/{assembler}/{source}/{sample_id}/{sample_id}_R2.fastq.gz",
        fq = assembly_input,
        bt = expand("results/map/{{assembler}}/{{filter_source}}/{{sample_id}}/final.fa.{index}.bt2l", index=range(1,5))
    output:
        sam = temp("results/map/{assembler}/{filter_source}/{sample_id}/{sample_id}.unsorted.sam")
    log: "results/map/{assembler}/{filter_source}/{sample_id}/bowtie.log"
    params:
        prefix = "results/map/{assembler}/{filter_source}/{sample_id}/final.fa",
        #tmp_out = "$TMPDIR/{sample_id}.bam",
        setting = config["bowtie2_params"]
    container: "docker://quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5"
    threads: 10
    resources:
        runtime = 60
    shell:
        """
        bowtie2 {params.setting} -p {threads} -x {params.prefix} \
         -1 {input.fq[0]} -2 {input.fq[0]} -S {output.sam} 2>{log}
        """

rule samtools_sort:
    input:
        sam = rules.bowtie.output.sam
    output:
        bam = "results/map/{assembler}/{filter_source}/{sample_id}/{sample_id}.bam"
    log: "results/map/{assembler}/{filter_source}/{sample_id}/samtools.log"
    container: "docker://quay.io/biocontainers/samtools:1.21--h96c455f_1"
    shell:
        """
        samtools view {input.sam} -h -b | samtools sort - -o {output.bam}
        """
    
rule gff2bed:
    input:
        gff = "results/annotation/{assembler}/{source}/{sample_id}/frame_selection/final.reformat.gff"
    output:
        bed = "results/annotation/{assembler}/{source}/{sample_id}/frame_selection/final.reformat.bed"
    conda:
        "../../envs/rseqc.yaml"
    shell:
        """
        gff2bed < {input.gff} > {output.bed}
        """

rule infer_experiment:
    input:
        bed = "results/annotation/{assembler}/{source}/{sample_id}/frame_selection/final.reformat.bed",
        bam = "results/map/{assembler}/{source}/{sample_id}/{sample_id}.bam"
    output:
        "results/annotation/{assembler}/{source}/{sample_id}/rseqc/rseqc.out"
    conda:
        "../../envs/rseqc.yaml"
    shell:
        """
        infer_experiment.py -i {input.bam} -r {input.bed}> {output[0]}
        """

rule multiqc_map_report:
    input:
        bt_logs = expand("results/map/{{assembler}}/{{source}}/{sample_id}/bowtie.log",
            sample_id = samples.keys()),
        fc_logs = expand("results/annotation/{{assembler}}/{{source}}/{sample_id}/featureCounts/fc.tab.summary",
            sample_id = samples.keys()),
        rs_logs = expand("results/annotation/{{assembler}}/{{source}}/{sample_id}/rseqc/rseqc.out",
            sample_id = samples.keys())
    output:
        "results/report/map/{assembler}_{source}_map_report.html",
        "results/report/map/{assembler}_{source}_map_report_data/multiqc_bowtie2.txt"
    params:
        dir = "qc",
        outdir = "results/report/map",
        config = "config/multiqc_{assembler}_config.yaml"
    run:
        assembler = wildcards.assembler
        source = wildcards.source
        for sample_id in samples.keys():
            tmp_dir = "{dir}/{sample_id}".format(dir=params.dir, sample_id=sample_id)
            shell("mkdir -p {tmp_dir}")
            bt_log = "results/map/{assembler}/{source}/{sample_id}/bowtie.log".format(
                assembler=assembler, source=source, sample_id=sample_id)
            fc_log = "results/annotation/{assembler}/{source}/{sample_id}/featureCounts/fc.tab.summary".format(
                assembler=assembler, source=source, sample_id=sample_id)
            rs_log = "results/annotation/{assembler}/{source}/{sample_id}/rseqc/rseqc.out".format(
                assembler=assembler, source=source, sample_id=sample_id)
            shell("cp {bt_log} {tmp_dir}/{sample_id}.bowtie.log")
            shell("cp {rs_log} {tmp_dir}/{sample_id}.rseqc.log")
            with open(fc_log, 'r') as fhin, open("{}/{}.fc.log".format(tmp_dir,sample_id), 'w') as fhout:
                for line in fhin:
                    if line[0:6]=="Status":
                        fhout.write("Status\t{}\n".format(sample_id))
                    else:
                        fhout.write(line)
        shell("cd {params.dir}; multiqc -f -c ../{params.config} -n {wildcards.assembler}_{wildcards.source}_map_report -o ../{params.outdir} .; cd -")
        shell("rm -rf {params.dir}")
