localrules: 
    multiqc_map_report, 
    multiqc_map_report_co, 
    wrap_assembly,
    wrap_assembly_co,

###############################
## Mapping for co-assemblies ##
###############################
rule kallisto_index_co:
    """
    Build kallisto index for co-assemblies.
    """
    output:
        "results/map/co-assembly/{assembler}/{assembly}/kallisto_index"
    input:
        "results/co-assembly/{assembler}/{assembly}/final.fa"
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
    """
    Quantify reads with kallisto.
    """
    output:
        h5 = "results/map/co-assembly/{assembler}/{assembly}/{sample_id}/kallisto/abundance.h5",
        tsv = "results/map/co-assembly/{assembler}/{assembly}/{sample_id}/kallisto/abundance.tsv",
        json = "results/map/co-assembly/{assembler}/{assembly}/{sample_id}/kallisto/run_info.json",
    input:
        index = rules.kallisto_index_co.output,
        R1 = lambda wildcards: map_dict[wildcards.sample_id]["R1"],
        R2 = lambda wildcards: map_dict[wildcards.sample_id]["R2"]
    log:
        "results/map/co-assembly/{assembler}/{assembly}/{sample_id}/kallisto/kallisto.log"
    container: "docker://quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1"
    conda: "../../envs/kallisto.yaml"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.h5),
        settings = config["kallisto_params"]
    threads: 2
    shell:
        """
        kallisto quant {params.settings} -t {threads} -i {input.index} -o {params.outdir} {input.R1} {input.R2} > {log} 2>&1
        """

rule parse_kallisto_co:
    """
    Parses Kallisto output
    """
    output:
        tsv="results/map/co-assembly/{assembler}/{assembly}/{sample_id}/kallisto/{quant_type}.tsv",
    input:
        tsv=rules.kallisto_quant_co.output.tsv
    run:
        df = pd.read_csv(input.tsv, sep="\t", index_col=0, header=0, comment="#", usecols=[0,3,4], names=["transcript_id", "est_counts", "tpm"])
        df = df.loc[:, [wildcards.quant_type]]
        df.columns=[wildcards.sample_id]
        df.to_csv(output.tsv, sep="\t")

rule collate_kallisto_co:
    """
    Collate Kallisto output
    """
    output:
        tsv="results/collated/co-assembly/{assembler}/{assembly}/abundance/kallisto/{quant_type}.tsv",
    input:
        expand("results/map/co-assembly/{{assembler}}/{{assembly}}/{sample_id}/kallisto/{{quant_type}}.tsv",
            sample_id = samples.keys())
    run:
        df = pd.DataFrame()
        for f in input:
            _df=pd.read_csv(f, sep="\t", index_col=0, header=0)
            df = pd.merge(df, _df, how="outer", left_index=True, right_index=True)
        df.to_csv(output.tsv, sep="\t")

#########################
## RSEM quantification ##
#########################

rule rsem_index_co:
    input:
        "results/co-assembly/{assembler}/{assembly}/final.fa"
    output:
        expand(
            "results/co-assembly/{{assembler}}/{{assembly}}/final.fa.{suff}",
            suff = ["bowtie2.ok", "RSEM.rsem.prepped.ok"])
    log:
        "results/co-assembly/{assembler}/{assembly}/rsem_index.log"
    params:
        trinity_mode = lambda wildcards: "--trinity_mode" if wildcards.assembler == "trinity" else "",
    threads: 10
    conda: "../../envs/trinity.yaml"
    shadow: "shallow"
    container: "docker://trinityrnaseq/trinityrnaseq:2.15.2"
    shell:
        """
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input} --aln_method bowtie2 \
            --est_method RSEM --thread_count {threads} {params.trinity_mode} --prep_reference >{log} 2>&1
        """

rule rsem_map_co:
    input:
        fa="results/co-assembly/{assembler}/{assembly}/final.fa",
        index=rules.rsem_index_co.output,
        R1 = lambda wildcards: map_dict[wildcards.sample_id]["R1"],
        R2 = lambda wildcards: map_dict[wildcards.sample_id]["R2"]
    log:
        "results/map/co-assembly/{assembler}/{assembly}/{sample_id}/rsem.log"
    output:
        genes="results/map/co-assembly/{assembler}/{assembly}/{sample_id}/RSEM/RSEM.genes.results",
        isoforms="results/map/co-assembly/{assembler}/{assembly}/{sample_id}/RSEM/RSEM.isoforms.results",
        isoforms_ok="results/map/co-assembly/{assembler}/{assembly}/{sample_id}/RSEM/RSEM.isoforms.results.ok",
        stat=expand(
            "results/map/co-assembly/{{assembler}}/{{assembly}}/{{sample_id}}/RSEM/RSEM.stat/RSEM.{suff}",
            suff = ["cnt","model","theta"]
        ),
        bam=temp("results/map/co-assembly/{assembler}/{assembly}/{sample_id}/RSEM/bowtie2.bam"),
        rsem_bam=temp("results/map/co-assembly/{assembler}/{assembly}/{sample_id}/RSEM/bowtie2.bam.for_rsem.bam"),
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.genes),
        trinity_mode = lambda wildcards: "--trinity_mode" if wildcards.assembler == "trinity" else "",
    threads: 4
    conda: "../../envs/trinity.yaml"
    container: "docker://trinityrnaseq/trinityrnaseq:2.15.2"
    shell:
        """
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input.fa} --seqType fq \
            --left {input.R1} --right {input.R2} --est_method RSEM {params.trinity_mode} \
            --aln_method bowtie2 --thread_count {threads} --output_dir {params.output_dir} > {log} 2>&1
        """

rule parse_rsem_co:
    input:
        res="results/map/co-assembly/{assembler}/{assembly}/{sample_id}/RSEM/RSEM.{rsem_res}.results"
    output:
        tsv="results/map/co-assembly/{assembler}/{assembly}/{sample_id}/RSEM/{rsem_res}.{quant_type}.tsv"
    run:
        import pandas as pd
        df = pd.read_csv(input.res, sep="\t", index_col=0)
        df = df.loc[:, [wildcards.quant_type]]
        df.columns = [wildcards.sample_id]
        df.to_csv(output.tsv, sep="\t")

rule collate_rsem_co:
    """
    Collate RSEM output
    """
    output:
        tsv="results/collated/co-assembly/{assembler}/{assembly}/abundance/RSEM/{rsem_res}.{quant_type}.tsv",
    input:
        expand("results/map/co-assembly/{{assembler}}/{{assembly}}/{sample_id}/RSEM/{{rsem_res}}.{{quant_type}}.tsv",
            sample_id = samples.keys())
    run:
        df = pd.DataFrame()
        for f in input:
            _df=pd.read_csv(f, sep="\t", index_col=0, header=0)
            df = pd.merge(df, _df, how="outer", left_index=True, right_index=True)
        df.to_csv(output.tsv, sep="\t")

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

rule multiqc_map_report_co:
    """
    Generate a multiqc report for co-assemblies.
    """
    input:
        kallisto_logs = expand("results/map/co-assembly/{{assembler}}/{{assembly}}/{sample_id}/kallisto/kallisto.log",
            sample_id = samples.keys()),
        kallisto = expand("results/map/co-assembly/{{assembler}}/{{assembly}}/{sample_id}/kallisto/abundance.h5",
            sample_id = samples.keys()),
        fc_logs = expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/featureCounts/{sample_id}.fc.tsv.summary",
            sample_id = samples.keys())
    output:
        "results/report/map/{assembler}_{assembly}_map_report.html"
    log:
        "results/report/map/{assembler}_{assembly}_map_report.log"
    params:
        config = workflow.source_path("../../config/multiqc_map_config.yaml"),
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        multiqc -f -c {params.config} -n {wildcards.assembler}_{wildcards.assembly}_map_report \
            -o {params.outdir} {input.kallisto_logs} {input.fc_logs} >{log} 2>&1
        """

#######################################
## Mapping for individual assemblies ##
#######################################
rule rsem_index:
    input:
        "results/assembly/{assembler}/{sample_id}/final.fa"
    output:
        expand(
            "results/assembly/{{assembler}}/{{sample_id}}/final.fa.{suff}",
            suff = ["bowtie2.ok", "RSEM.rsem.prepped.ok"])
    log:
        "results/map/{assembler}/{sample_id}/rsem_index.log"
    params:
        trinity_mode = lambda wildcards: "--trinity_mode" if wildcards.assembler == "trinity" else "",
    threads: 10
    conda: "../../envs/trinity.yaml"
    shadow: "shallow"
    container: "docker://trinityrnaseq/trinityrnaseq:2.15.2"
    shell:
        """
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input} --aln_method bowtie2 \
            --est_method RSEM --thread_count {threads} {params.trinity_mode} --prep_reference >{log} 2>&1
        """

rule rsem_map:
    input:
        fa="results/assembly/{assembler}/{sample_id}/final.fa",
        index=rules.rsem_index.output,
        R1 = lambda wildcards: map_dict[wildcards.sample_id]["R1"],
        R2 = lambda wildcards: map_dict[wildcards.sample_id]["R2"]
    log:
        "results/map/map/{assembler}/{sample_id}/RSEM/rsem.log"
    output:
        genes="results/map/{assembler}/{sample_id}/RSEM/RSEM.genes.results",
        isoforms="results/map/{assembler}/{sample_id}/RSEM/RSEM.isoforms.results",
        isoforms_ok="results/map/{assembler}/{sample_id}/RSEM/RSEM.isoforms.results.ok",
        stat=expand(
            "results/map/{{assembler}}/{{sample_id}}/RSEM/RSEM.stat/RSEM.{suff}",
            suff = ["cnt","model","theta"]
        ),
        bam=temp("results/map/{assembler}/{sample_id}/RSEM/bowtie2.bam"),
        rsem_bam=temp("results/map/{assembler}/{sample_id}/RSEM/bowtie2.bam.for_rsem.bam"),
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.genes),
        trinity_mode = lambda wildcards: "--trinity_mode" if wildcards.assembler == "trinity" else "",
    threads: 4
    conda: "../../envs/trinity.yaml"
    container: "docker://trinityrnaseq/trinityrnaseq:2.15.2"
    shell:
        """
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input.fa} --seqType fq \
            --left {input.R1} --right {input.R2} --est_method RSEM {params.trinity_mode} \
            --aln_method bowtie2 --thread_count {threads} --output_dir {params.output_dir} > {log} 2>&1
        """

rule parse_rsem:
    input:
        res="results/map/{assembler}/{sample_id}/RSEM/RSEM.{rsem_res}.results"
    output:
        tsv="results/map/{assembler}/{sample_id}/RSEM/{rsem_res}.{quant_type}.tsv"
    run:
        import pandas as pd
        df = pd.read_csv(input.res, sep="\t", index_col=0)
        df = df.loc[:, [wildcards.quant_type]]
        df.columns = [wildcards.sample_id]
        df.to_csv(output.tsv, sep="\t")

rule kallisto_index:
    """
    Build kallisto index for individual assemblies.
    """
    output:
        "results/map/{assembler}/{sample_id}/kallisto_index"
    input:
        "results/assembly/{assembler}/{sample_id}/final.fa"
    log:
        "results/map/{assembler}/{sample_id}/kallisto_index.log"
    container: "docker://quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1"
    conda: "../../envs/kallisto.yaml"
    threads: 2
    shell:
        """
        kallisto index -t {threads} -i {output} {input} > {log} 2>&1
        """

rule kallisto_quant:
    """
    Quantify reads with kallisto.
    """
    input:
        index = rules.kallisto_index.output,
        R1 = lambda wildcards: map_dict[wildcards.sample_id]["R1"],
        R2 = lambda wildcards: map_dict[wildcards.sample_id]["R2"]
    output:
        h5 = "results/map/{assembler}/{sample_id}/kallisto/abundance.h5",
        tsv = "results/map/{assembler}/{sample_id}/kallisto/abundance.tsv",
        json = "results/map/{assembler}/{sample_id}/kallisto/run_info.json",
    log:
        "results/map/{assembler}/{sample_id}/kallisto/kallisto.log"
    container: "docker://quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1"
    conda: "../../envs/kallisto.yaml"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.h5),
        settings = config["kallisto_params"]
    threads: 2
    shell:
        """
        kallisto quant {params.settings} -t {threads} -i {input.index} -o {params.outdir} {input.R1} {input.R2} > {log} 2>&1
        """

rule parse_kallisto:
    """
    Parses Kallisto output
    """
    output:
        tsv="results/map/{assembler}/{sample_id}/kallisto/{quant_type}.tsv",
    input:
        tsv=rules.kallisto_quant.output.tsv
    run:
        df = pd.read_csv(input.tsv, sep="\t", index_col=0, header=0, comment="#", usecols=[0,3,4], names=["transcript_id", "est_counts", "tpm"])
        df = df.loc[:, [wildcards.quant_type]]
        df.columns=[wildcards.sample_id]
        df.to_csv(output.tsv, sep="\t")

rule wrap_assembly:
    """
    Subread has a line length limit on fasta files, to be sure we are under this limit we wrap the sequences.
    """
    input:
        "results/assembly/{assembler}/{sample_id}/final.fa"
    output:
        temp("results/assembly/{assembler}/{sample_id}/final.wrap.fa")
    log:
        "results/assembly/{assembler}/{sample_id}/wrap_assembly.log"
    shell:
        """
        seqkit seq -w 60 {input} > {output} 2> {log}
        """

rule multiqc_map_report:
    """
    Generate a multiqc report for individual assemblies.
    """
    input:
        kallisto_logs = expand("results/map/{{assembler}}/{sample_id}/kallisto/kallisto.log",
            sample_id = samples.keys()),
        kallisto = expand("results/map/{{assembler}}/{sample_id}/kallisto/abundance.h5",
            sample_id = samples.keys()),
        fc_logs = expand("results/annotation/{{assembler}}/{sample_id}/featureCounts/fc.tsv.summary",
            sample_id = samples.keys()),
    output:
        "results/report/map/{assembler}_map_report.html",
    log:
        "results/report/map/{assembler}_map_report.log"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        config = workflow.source_path("../../config/multiqc_map_config.yaml")
    shell:
        """
        multiqc -f -c {params.config} -n {wildcards.assembler}_map_report \
            -o {params.outdir} {input.kallisto_logs} {input.kallisto} {input.fc_logs} >{log} 2>&1
        """