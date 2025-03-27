localrules: 
    assembly_stats,
    co_assembly_stats

##############
## Assembly ##
##############

rule transabyss:
    """
    Assemble fungal reads with transabyss at a certain k-mer size.
    """
    output:
        fa = "results/assembly/transabyss/{sample_id}/{k}/{sample_id}.{k}-final.fa",
    input:
        R1=lambda wildcards: map_dict[wildcards.sample_id]["R1"],
        R2=lambda wildcards: map_dict[wildcards.sample_id]["R2"]
    log: "results/assembly/transabyss/{sample_id}/{k}/log"
    conda: "../../envs/transabyss.yaml"
    container: "docker://quay.io/biocontainers/transabyss:2.0.1--pyh864c0ab_7"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        tmpdir = "$TMPDIR/{sample_id}.{k}.transabyss",
        R1 = "$TMPDIR/{sample_id}.{k}.transabyss/R1.fastq",
        R2 = "$TMPDIR/{sample_id}.{k}.transabyss/R2.fastq",
        min_contig_len = config["min_contig_len"]
    threads: 16
    resources:
        runtime = 60 * 24
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.R1} > {params.R1}
        gunzip -c {input.R2} > {params.R2}
        transabyss --pe {params.R1} {params.R2} -k {wildcards.k} --length {params.min_contig_len} \
            --outdir {params.tmpdir} --name {wildcards.sample_id}.{wildcards.k} \
            --threads {threads} >{log} 2>&1
        mv {params.tmpdir}/* {params.outdir}
        rm -rf {params.tmpdir}
        """

rule transabyss_merge:
    """
    Merge transabyss assemblies from different k-mer sizes.
    """
    output:
        fa = "results/assembly/transabyss/{sample_id}/final.fa",
    input:
        fungi_input,
        expand("results/assembly/transabyss/{{sample_id}}/{k}/{{sample_id}}.{k}-final.fa",
            k = config["transabyss_kmers"])
    log:
        "results/logs/transabyss/{sample_id}.merge.log"
    conda: "../../envs/transabyss.yaml"
    container: "docker://quay.io/biocontainers/transabyss:2.0.1--pyh864c0ab_7"
    params:
        mink = min(config["transabyss_kmers"]),
        maxk = max(config["transabyss_kmers"]),
        i = lambda wildcards, input: sorted(input[2:]),
        prefix = ["k{x}.".format(x=x) for x in sorted(config["transabyss_kmers"])],
        tmpout = "$TMPDIR/{sample_id}.ta.merged.fa"
    threads: 16
    resources:
        runtime = 60 * 24
    shell:
        """
        rm -rf {params.tmpout}
        transabyss-merge {params.i} --mink {params.mink} --maxk {params.maxk} \
            --out {params.tmpout} --threads {threads} --force --prefixes {params.prefix} > {log} 2>&1
        mv {params.tmpout} {output.fa}
        """

rule trinity:
    """
    Assemble fungal reads with Trinity.
    """
    output:
        fa="results/assembly/trinity/{sample_id}/final.fa",
    input:
        R1=lambda wildcards: map_dict[wildcards.sample_id]["R1"],
        R2=lambda wildcards: map_dict[wildcards.sample_id]["R2"]
    log: "results/assembly/trinity/{sample_id}/log"
    params:
        min_contig_len = config["min_contig_len"],
        tmpdir="$TMPDIR/{sample_id}.trinity",
        outdir=lambda wildcards, output: os.path.dirname(output.fa),
        out_base = lambda wildcards, output: os.path.basename(output.fa),
        max_mem = lambda wildcards, resources: int(resources.mem_mb /1000)
    threads: 6
    resources:
        runtime = 24 * 60,
        mem_mb = 5000
    conda: "../../envs/trinity.yaml"
    container: "docker://trinityrnaseq/trinityrnaseq:2.15.2"
    shell:
        """
        rm -rf {params.outdir}/*
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}
        Trinity --CPU {threads} --min_contig_length {params.min_contig_len} \
            --output {params.tmpdir} --left {input.R1} --right {input.R2} \
            --seqType fq --max_memory {params.max_mem}G > {log} 2>&1
        mv {params.tmpdir}/* {params.outdir}/
        mv {params.tmpdir}.Trinity.fasta {output.fa}
        rm -rf {params.tmpdir}
        """

rule megahit:
    """
    Assemble fungal reads with Megahit.
    """
    output:
        fa = "results/assembly/megahit/{sample_id}/final.fa"
    input:
        fungi_input
    log: "results/assembly/megahit/{sample_id}/log"
    params:
        min_contig_len = config["min_contig_len"],
        out_dir = "results/assembly/megahit/{sample_id}",
        tmp_dir = "$TMPDIR/megahit/{sample_id}",
        tmp_dir_base = "$TMPDIR/megahit"
    conda: "../../envs/megahit.yaml"
    container: "docker://quay.io/biocontainers/megahit:1.2.9--h43eeafb_5"
    threads: 10
    resources:
        runtime = 60 * 10
    shell:
        """
        mkdir -p {params.tmp_dir_base}
        megahit -1 {input[0]} -2 {input[1]} --prune-level 3 \
            --min-contig-len {params.min_contig_len} -o {params.tmp_dir} \
            -t {threads} > {log} 2>&1
        mv {params.tmp_dir}/final.contigs.fa {output.fa}
        mv {params.tmp_dir}/opt* {params.out_dir}
        rm -rf {params.tmp_dir}
        """

#################
## CO-ASSEMBLY ##
#################

def trinity_partition(wildcards):
    """
    Read memory allocated (in config) to trinity insilico normalization
    multiply by a factor 1.28 and convert to MB
    Then figure out slurm partition to use and return it

    Parameters
    ----------
    wildcards : snakemake.Wildcards
        Wildcards object
    
    Returns
    -------
    str
        Slurm partition to use
    """
    mb = config["insilico_norm_mem"] * 1.28 * 1000
    return slurm_mem_partition(mb)

rule trinity_normalize:
    """
    Normalize reads with Trinity in-silico normalization.
    """
    output:
        R1="results/in-silico-normalization/{assembly}/left.norm.fq.gz",
        R2="results/in-silico-normalization/{assembly}/right.norm.fq.gz",
    input:
        R1 = lambda wildcards: assemblies[wildcards.assembly]["R1"],
        R2 = lambda wildcards: assemblies[wildcards.assembly]["R2"]
    log:
        "results/in-silico-normalization/{assembly}/log"
    params:
        max_cov = config["insilico_norm_max_cov"],
        R1 = lambda wildcards: ",".join(sorted(assemblies[wildcards.assembly]["R1"])),
        R2 = lambda wildcards: ",".join(sorted(assemblies[wildcards.assembly]["R2"])),
        tmpdir="$TMPDIR/{assembly}.norm",
        mem=config["insilico_norm_mem"],
    resources:
        mem_mb = config["insilico_norm_mem"] * 1.28 * 1000
    conda: "../../envs/trinity.yaml"
    container: "docker://trinityrnaseq/trinityrnaseq:2.15.2"
    threads: 10
    shell:
        """
        mkdir -p {params.tmpdir}
        echo -e {params.R1} | tr "," "\n" > {params.tmpdir}/R1.list
        echo -e {params.R2} | tr "," "\n" > {params.tmpdir}/R2.list
        $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM {params.mem}G --max_cov {params.max_cov} \
            --left_list {params.tmpdir}/R1.list --right_list {params.tmpdir}/R2.list \
            --pairs_together --PARALLEL_STATS --CPU {threads} --output {params.tmpdir} --tmp_dir_name out > {log} 2>&1
        gzip -c {params.tmpdir}/left.norm.fq > {params.tmpdir}/left.norm.fq.gz
        gzip -c {params.tmpdir}/right.norm.fq > {params.tmpdir}/right.norm.fq.gz
        mv {params.tmpdir}/left.norm.fq.gz {output.R1}
        mv {params.tmpdir}/right.norm.fq.gz {output.R2}
        rm -r {params.tmpdir}
        """

rule transabyss_co:
    """
    Co-assemble fungal reads with transabyss at a certain k-mer size.
    """
    output:
        "results/co-assembly/transabyss/{assembly}/{k}/{assembly}-final.fa",
    input:
        R1=rules.trinity_normalize.output.R1,
        R2=rules.trinity_normalize.output.R2
    log: "results/co-assembly/transabyss/{assembly}/{k}/log"
    conda: "../../envs/transabyss.yaml"
    container: "docker://quay.io/biocontainers/transabyss:2.0.1--pyh864c0ab_7"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        tmpdir="$TMPDIR/{assembly}.{k}.transabyss",
        min_contig_len=config["min_contig_len"],
        R1="$TMPDIR/{assembly}.{k}.transabyss/R1.fastq",
        R2="$TMPDIR/{assembly}.{k}.transabyss/R2.fastq"
    threads: 16
    resources:
        runtime = 60 * 24
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.R1} > {params.R1}
        gunzip -c {input.R2} > {params.R2}
        transabyss --pe {params.R1} {params.R2} -k {wildcards.k} --length {params.min_contig_len} \
            --outdir {params.tmpdir} --name {wildcards.assembly} \
            --threads {threads} >{log} 2>&1
        mv {params.tmpdir}/* {params.outdir}
        rm -rf {params.tmpdir}
        """

rule transabyss_merge_co:
    """
    Merge transabyss co-assemblies from different k
    """
    output:
        fa = "results/co-assembly/transabyss/{assembly}/final.fa"
    input:
        expand("results/co-assembly/transabyss/{{assembly}}/{k}/{{assembly}}-final.fa",
            k = config["transabyss_kmers"])
    log: "results/co-assembly/transabyss/{assembly}/merge.log"
    conda: "../../envs/transabyss.yaml"
    container: "docker://quay.io/biocontainers/transabyss:2.0.1--pyh864c0ab_7"
    params:
        mink = min(config["transabyss_kmers"]),
        maxk=max(config["transabyss_kmers"]),
        i=lambda wildcards, input: sorted(input),
        prefix=["k{x}.".format(x=x) for x in
                sorted(config["transabyss_kmers"])],
        tmpout="$TMPDIR/{assembly}.ta.merged.fa"
    threads: 16
    resources:
        runtime = 60 * 10
    shell:
        """
        rm -rf {params.tmpout}
        wd=$(pwd)
        transabyss-merge {params.i} --mink {params.mink} --maxk {params.maxk} \
            --out {params.tmpout} --threads {threads} --force --prefixes {params.prefix} > {log} 2>&1
        mv {params.tmpout} {output.fa}
        """

rule trinity_co:
    """
    Co-assemble fungal reads with Trinity.
    """
    output:
        fa = "results/co-assembly/trinity/{assembly}/final.fa"
    input:
        R1=rules.trinity_normalize.output.R1,
        R2=rules.trinity_normalize.output.R2
    log: "results/co-assembly/trinity/{assembly}/log"
    params:
        min_contig_len = config["min_contig_len"],
        tmpdir="$TMPDIR/{assembly}.co.trinity",
        outdir = lambda wildcards, output: os.path.dirname(output.fa),
        out_base = lambda wildcards, output: os.path.basename(output.fa),
        max_mem = lambda wildcards, resources: int(resources.mem_mb /1000)
    threads: 6
    resources:
        runtime = 24 * 60,
        mem_mb = 5000
    conda: "../../envs/trinity.yaml"
    container: "docker://trinityrnaseq/trinityrnaseq:2.15.2"
    shell:
        """
        rm -rf {params.outdir}/*
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}
        Trinity --no_normalize_reads --CPU {threads} --min_contig_length {params.min_contig_len} \
            --output {params.tmpdir} --left {input.R1} --right {input.R2} \
            --seqType fq --max_memory {params.max_mem}G > {log} 2>&1
        mv {params.tmpdir}/* {params.outdir}/
        mv {params.tmpdir}.Trinity.fasta {params.outdir}/{params.out_base}
        rm -rf {params.tmpdir}
        """

rule megahit_co:
    """
    Co-assemble fungal reads with Megahit.
    """
    output:
        fa = "results/co-assembly/megahit/{assembly}/final.fa"
    input:
        R1=rules.trinity_normalize.output.R1,
        R2=rules.trinity_normalize.output.R2
    log: "results/co-assembly/megahit/{assembly}/log"
    params:
        min_contig_len = config["min_contig_len"],
        out_dir = "results/co-assembly/megahit/{assembly}",
        tmp_dir = "$TMPDIR/megahit/{assembly}",
        tmp_dir_base = "$TMPDIR/megahit"
    conda: "../../envs/megahit.yaml"
    container: "docker://quay.io/biocontainers/megahit:1.2.9--h43eeafb_5"
    threads: 20
    resources:
        runtime = 24 * 60
    shell:
        """
        wd=$(pwd)
        mkdir -p {params.tmp_dir_base}
        megahit -1 {input.R1} -2 {input.R2} --prune-level 3 \
            --min-contig-len {params.min_contig_len} -o {params.tmp_dir} \
            -t {threads} > {log} 2>&1
        mv {params.tmp_dir}/final.contigs.fa {output.fa}
        mv {params.tmp_dir}/* {params.out_dir}
        rm -r {params.tmp_dir}
        """

rule assembly_stats:
    """
    Calculate assembly statistics.
    """
    output:
        "results/report/assembly/{assembler}_stats.tsv",
        "results/report/assembly/{assembler}_size_dist.tsv"
    input:
        expand("results/assembly/{{assembler}}/{sample_id}/final.fa",
            sample_id = samples.keys())
    run:
        names = [x.split("/")[-2] for x in input]
        shell("python source/utils/assembly_stats.py -i {input} -n {names} --size-dist-file {output[1]} > {output[0]}")

rule co_assembly_stats:
    """
    Calculate assembly statistics for co-assembly
    """
    output:
        "results/report/co-assembly/{assembler}.{assembly}_assembly_stats.tsv",
        "results/report/co-assembly/{assembler}.{assembly}_assembly_size_dist.tsv"
    input:
        "results/co-assembly/{assembler}/{assembly}/final.fa"
    shell:
        """
        python source/utils/assembly_stats.py -i {input[0]} -n {wildcards.assembler}.{wildcards.assembly} --size-dist-file {output[1]} > {output[0]}
        """