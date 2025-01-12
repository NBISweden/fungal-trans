localrules: 
    multiqc_map_report, 
    multiqc_map_report_co, 
    wrap_assembly,
    wrap_assembly_co

###############################
## Mapping for co-assemblies ##
###############################
rule kallisto_index_co:
    input:
        "results/co-assembly/{assembler}/{assembly}/final.fa"
    output:
        "results/map/co-assembly/{assembler}/{assembly}/kallisto_index"
    log:
        "results/map/co-assembly/{assembler}/{assembly}/kallisto_index.log"
    container: "docker://quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1"
    conda: "../../envs/kallisto.yaml"
    threads: 2
    shell:
        """
        kallisto index -t {threads} -i {output} {input} > {log} 2>&1
        """

rule kallisto_quant_co:
    input:
        index = rules.kallisto_index_co.output,
        R1 = lambda wildcards: map_dict[wildcards.sample_id]["R1"],
        R2 = lambda wildcards: map_dict[wildcards.sample_id]["R2"]
    output:
        h5 = "results/map/co-assembly/{assembler}/{assembly}/{sample_id}/abundance.h5",
        tsv = "results/map/co-assembly/{assembler}/{assembly}/{sample_id}/abundance.tsv",
        json = "results/map/co-assembly/{assembler}/{assembly}/{sample_id}/run_info.json",
    log:
        "results/map/co-assembly/{assembler}/{assembly}/{sample_id}/kallisto.log"
    container: "docker://quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1"
    conda: "../../envs/kallisto.yaml"
    params:
        bootstrap = 100,
        outdir = lambda wildcards, output: os.path.dirname(output.h5)
    threads: 2
    shell:
        """
        kallisto quant -b {params.bootstrap} -t {threads} -i {input.index} -o {params.outdir} {input.R1} {input.R2} > {log} 2>&1
        """

rule wrap_assembly_co:
    """
    Subread has a line length limit on fasta files, to be sure we are under this limit we wrap the sequences.
    """
    input:
        "results/co-assembly/{assembler}/{assembly}/final.fa"
    output:
        temp("results/co-assembly/{assembler}/{assembly}/final.wrap.fa")
    log:
        "results/co-assembly/{assembler}/{assembly}/wrap_assembly.log"
    shell:
        """
        seqkit seq -w 60 {input} > {output} 2> {log}
        """

rule subread_index_co:
    input:
        rules.wrap_assembly_co.output
    output:
        expand("results/map/co-assembly/{{assembler}}/{{assembly}}/subread_index.{suff}",
                suff = ["00.b.array", "00.b.tab", "files", "log", "lowinf", "reads"])
    log:
        "results/map/co-assembly/{assembler}/{assembly}/subread_index.log"
    container: "docker://quay.io/biocontainers/subread:2.0.8--h577a1d6_0"
    conda: "../../envs/featurecount.yaml"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        subread-buildindex -o {params.outdir}/subread_index {input} > {log} 2>&1
        """

rule subread_align_co:
    input:
        index = rules.subread_index_co.output,
        R1 = lambda wildcards: map_dict[wildcards.sample_id]["R1"],
        R2 = lambda wildcards: map_dict[wildcards.sample_id]["R2"]
    output:
        "results/map/co-assembly/{assembler}/{assembly}/{sample_id}.bam"
    log:
        "results/map/co-assembly/{assembler}/{assembly}/{sample_id}.subread_align.log"
    container: "docker://quay.io/biocontainers/subread:2.0.8--h577a1d6_0"
    conda: "../../envs/featurecount.yaml"
    params:
        index = lambda wildcards, input: os.path.join(os.path.dirname(input.index[0]), "subread_index"),
    threads: 4
    shell:
        """
        subread-align -T {threads} -sortReadsByCoordinates -r {input.R1} -R {input.R2} -i {params.index} -o {output} -t 0 > {log} 2>&1
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
rule kallisto_index:
    input:
        "results/assembly/{assembler}/{filter_source}/{sample_id}/final.fa"
    output:
        "results/map/{assembler}/{filter_source}/{sample_id}/kallisto_index"
    log:
        "results/map/{assembler}/{filter_source}/{sample_id}/kallisto_index.log"
    container: "docker://quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1"
    conda: "../../envs/kallisto.yaml"
    threads: 2
    shell:
        """
        kallisto index -t {threads} -i {output} {input} > {log} 2>&1
        """

rule kallisto_quant:
    input:
        index = rules.kallisto_index.output,
        fq = assembly_input
    output:
        h5 = "results/map/{assembler}/{filter_source}/{sample_id}/abundance.h5",
        tsv = "results/map/{assembler}/{filter_source}/{sample_id}/abundance.tsv",
        json = "results/map/{assembler}/{filter_source}/{sample_id}/run_info.json",
    log:
        "results/map/{assembler}/{filter_source}/{sample_id}/kallisto.log"
    container: "docker://quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1"
    conda: "../../envs/kallisto.yaml"
    params:
        bootstrap = 100,
        outdir = lambda wildcards, output: os.path.dirname(output.h5),
        settings = config["kallisto_params"]
    threads: 2
    shell:
        """
        kallisto quant {params.settings} -b {params.bootstrap} -t {threads} -i {input.index} -o {params.outdir} {input.fq[0]} {input.fq[1]} > {log} 2>&1
        """

rule wrap_assembly:
    """
    Subread has a line length limit on fasta files, to be sure we are under this limit we wrap the sequences.
    """
    input:
        "results/assembly/{assembler}/{filter_source}/{sample_id}/final.fa"
    output:
        temp("results/assembly/{assembler}/{filter_source}/{sample_id}/final.wrap.fa")
    log:
        "results/assembly/{assembler}/{filter_source}/{sample_id}/wrap_assembly.log"
    shell:
        """
        seqkit seq -w 60 {input} > {output} 2> {log}
        """

rule subread_index:
    input:
        rules.wrap_assembly.output
    output:
        expand("results/map/{{assembler}}/{{filter_source}}/{{sample_id}}/subread_index.{suff}",
                suff = ["00.b.array", "00.b.tab", "files", "log", "lowinf", "reads"])
    log:
        "results/map/{assembler}/{filter_source}/{sample_id}/subread_index.log"
    container: "docker://quay.io/biocontainers/subread:2.0.8--h577a1d6_0"
    conda: "../../envs/featurecount.yaml"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        subread-buildindex -o {params.outdir}/subread_index {input} > {log} 2>&1
        """

rule subread_align:
    input:
        index = rules.subread_index.output,
        fq = assembly_input
    output:
        "results/map/{assembler}/{filter_source}/{sample_id}/{sample_id}.bam"
    log:
        "results/map/{assembler}/{filter_source}/{sample_id}/subread_align.log"
    container: "docker://quay.io/biocontainers/subread:2.0.8--h577a1d6_0"
    conda: "../../envs/featurecount.yaml"
    params:
        index = lambda wildcards, input: os.path.join(os.path.dirname(input.index[0]), "subread_index"),
    threads: 4
    shell:
        """
        subread-align -T {threads} -sortReadsByCoordinates -r {input.fq[0]} -R {input.fq[1]} -i {params.index} -o {output} -t 0 > {log} 2>&1
        """

rule multiqc_map_report:
    input:
        bt_logs = expand("results/map/{{assembler}}/{{filter_source}}/{sample_id}/bowtie.log",
            sample_id = samples.keys()),
        fc_logs = expand("results/annotation/{{assembler}}/{{filter_source}}/{sample_id}/featureCounts/fc.tab.summary",
            sample_id = samples.keys()),
        rs_logs = expand("results/annotation/{{assembler}}/{{filter_source}}/{sample_id}/rseqc/rseqc.out",
            sample_id = samples.keys())
    output:
        "results/report/map/{assembler}_{filter_source}_map_report.html",
        "results/report/map/{assembler}_{filter_source}_map_report_data/multiqc_bowtie2.txt"
    params:
        dir = "qc",
        outdir = "results/report/map",
        config = "config/multiqc_{assembler}_config.yaml"
    run:
        assembler = wildcards.assembler
        source = wildcards.filter_source
        for sample_id in samples.keys():
            tmp_dir = "{dir}/{sample_id}".format(dir=params.dir, sample_id=sample_id)
            shell("mkdir -p {tmp_dir}")
            bt_log = "results/map/{assembler}/{filter_source}/{sample_id}/bowtie.log".format(
                assembler=assembler, source=source, sample_id=sample_id)
            fc_log = "results/annotation/{assembler}/{filter_source}/{sample_id}/featureCounts/fc.tab.summary".format(
                assembler=assembler, source=source, sample_id=sample_id)
            rs_log = "results/annotation/{assembler}/{filter_source}/{sample_id}/rseqc/rseqc.out".format(
                assembler=assembler, source=source, sample_id=sample_id)
            shell("cp {bt_log} {tmp_dir}/{sample_id}.bowtie.log")
            shell("cp {rs_log} {tmp_dir}/{sample_id}.rseqc.log")
            with open(fc_log, 'r') as fhin, open("{}/{}.fc.log".format(tmp_dir,sample_id), 'w') as fhout:
                for line in fhin:
                    if line[0:6]=="Status":
                        fhout.write("Status\t{}\n".format(sample_id))
                    else:
                        fhout.write(line)
        shell("cd {params.dir}; multiqc -f -c ../{params.config} -n {wildcards.assembler}_{wildcards.filter_source}_map_report -o ../{params.outdir} .; cd -")
        shell("rm -rf {params.dir}")
