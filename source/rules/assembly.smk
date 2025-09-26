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
    shadow: "shallow"
    container: "docker://quay.io/biocontainers/transabyss:2.0.1--pyh864c0ab_7"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        tmpdir = "{sample_id}.{k}.transabyss",
        R1 = "{sample_id}.{k}.transabyss/R1.fastq",
        R2 = "{sample_id}.{k}.transabyss/R2.fastq",
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
    shadow: "shallow"
    container: "docker://quay.io/biocontainers/transabyss:2.0.1--pyh864c0ab_7"
    params:
        mink = min(config["transabyss_kmers"]),
        maxk = max(config["transabyss_kmers"]),
        i = lambda wildcards, input: sorted(input[2:]),
        prefix = ["k{x}.".format(x=x) for x in sorted(config["transabyss_kmers"])],
        tmpout = "{sample_id}.ta.merged.fa"
    threads: 16
    resources:
        runtime = 60 * 24
    shell:
        """
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
        tmpdir="{sample_id}.trinity",
        outdir=lambda wildcards, output: os.path.dirname(output.fa),
        out_base = lambda wildcards, output: os.path.basename(output.fa),
        max_mem = lambda wildcards, resources: int(resources.mem_mb /1000)
    threads: 6
    resources:
        runtime = 24 * 60,
        mem_mb = 5000
    conda: "../../envs/trinity.yaml"
    shadow: "shallow"
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
        tmp_dir = "megahit/{sample_id}",
        tmp_dir_base = "megahit"
    conda: "../../envs/megahit.yaml"
    shadow: "shallow"
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

import os
from datetime import datetime

def get_or_create_timestamp(path="tmp/.trinity_timestamp"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path):
        with open(path) as f:
            return f.read().strip()
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        with open(path, "w") as f:
            f.write(ts)
        return ts

timestamp = get_or_create_timestamp()
#timestamp = "20250614_214415"

tmpdir = lambda wildcards: os.path.join(os.environ.get("TMPDIR", "/tri"), f"{wildcards.assembly}.co.trinity_{timestamp}")


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
        outdir=lambda wildcards, output: os.path.dirname(output.R1),
        max_mem=config["insilico_norm_mem"],
        tmpdir=tmpdir
    resources:
        mem_mb = int(config["insilico_norm_mem"] * 1.27 * 1000)
    conda: "../../envs/trinity.yaml"
    container: "docker://trinityrnaseq/trinityrnaseq:2.15.2"
    shadow: "shallow"
    threads: 16
    shell:
        """
        # echo -e {params.R1} > R1.list
        # echo -e {params.R2} > R2.list
        Trinity --just_normalize_reads --normalize_by_read_set --seqType fq --normalize_max_read_cov {params.max_cov} \
            --left {params.R1} --right {params.R2} --max_memory {params.max_mem}G \
            --CPU {threads} --output {params.tmpdir} > {log} 2>&1
        gzip -c {params.tmpdir}/insilico_read_normalization_altogether/left.norm.fq > {output.R1}
        gzip -c {params.tmpdir}/insilico_read_normalization_altogether/right.norm.fq > {output.R2}
        #ln -s $(readlink -f {params.tmpdir}/insilico_read_normalization_altogether/left.norm.fq) {output.R1}
        #ln -s $(readlink -f {params.tmpdir}/insilico_read_normalization_altogether/right.norm.fq) {output.R2}
        #touch {output.R1}.done {output.R2}.done
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
        tmpdir="{assembly}.{k}.transabyss",
        min_contig_len=config["min_contig_len"],
        R1="{assembly}.{k}.transabyss/R1.fastq",
        R2="{assembly}.{k}.transabyss/R2.fastq"
    threads: 16
    shadow: "shallow"
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
    shadow: "shallow"
    container: "docker://quay.io/biocontainers/transabyss:2.0.1--pyh864c0ab_7"
    params:
        mink = min(config["transabyss_kmers"]),
        maxk=max(config["transabyss_kmers"]),
        i=lambda wildcards, input: sorted(input),
        prefix=["k{x}.".format(x=x) for x in
                sorted(config["transabyss_kmers"])],
        tmpout="{assembly}.ta.merged.fa"
    threads: 16
    resources:
        runtime = 60 * 10
    shell:
        """
        wd=$(pwd)
        transabyss-merge {params.i} --mink {params.mink} --maxk {params.maxk} \
            --out {params.tmpout} --threads {threads} --force --prefixes {params.prefix} > {log} 2>&1
        mv {params.tmpout} {output.fa}
        """

rule trinity_co_inchworm:
    """
    Run only the Inchworm step of Trinity.
    Produces intermediate assembly files in tmpdir, to be used by Chrysalis.
    """
    input:
        R1=rules.trinity_normalize.output.R1,
        R2=rules.trinity_normalize.output.R2
    output:
        inchworm_finished="results/co-assembly/trinity/{assembly}/trinity.multi-stage/inchworm.DS.fa.finished"
    log:
        "results/co-assembly/trinity/{assembly}/trinity.multi-stage/log_inchworm"
    params:
        min_contig_len=config["min_contig_len"],
        outdir=lambda wildcards, output: os.path.dirname(output.inchworm_finished),
        max_mem = lambda wildcards, resources: int(resources.mem_mb / 1000)
    conda: "../../envs/trinity.yaml"
    container: "docker://trinityrnaseq/trinityrnaseq:2.15.2"
    shell:
        """
        Trinity \
            --no_run_chrysalis --min_kmer_cov 2 \
            --min_contig_length {params.min_contig_len} \
            --max_memory {params.max_mem}G \
            --CPU {threads} \
            --left {input.R1} --right {input.R2} \
            --output {params.outdir} \
            --seqType fq \
            --no_normalize_reads \
            > {log} 2>&1
        """

rule trinity_co_chrysalis:
    """
    Run the Chrysalis step of Trinity, using the output from Inchworm.
    Produces the recursive_trinity.cmds file for downstream steps.
    """
    input:
        inchworm_finished=rules.trinity_co_inchworm.output.inchworm_finished,
        R1=rules.trinity_normalize.output.R1,
        R2=rules.trinity_normalize.output.R2
    output:
        cmds="results/co-assembly/trinity/{assembly}/trinity.multi-stage/recursive_trinity.cmds"
    log:
        "results/co-assembly/trinity/{assembly}/trinity.multi-stage/log_chrysalis"
    params:
        min_contig_len = config["min_contig_len"],
        max_mem = lambda wildcards, resources: int(resources.mem_mb / 1000),
        outdir=lambda wildcards, output: os.path.dirname(output.cmds)
    conda: "../../envs/trinity.yaml"
    container: "docker://trinityrnaseq/trinityrnaseq:2.15.2"
    shell:
        """
        Trinity \
            --no_distributed_trinity_exec \
            --no_bowtie --min_kmer_cov 2 \
            --min_contig_length {params.min_contig_len} \
            --max_memory {params.max_mem}G \
            --CPU {threads} \
            --seqType fq \
            --no_normalize_reads \
            --left {input.R1} --right {input.R2} \
            --output {params.outdir} \
            > {log} 2>&1
        """

checkpoint trinity_co_butterfly_split:
    """
    Prepare for parallel Butterfly execution by splitting commands.
    Splits the cmds file into 999 separate files.
    """
    input:
        rules.trinity_co_chrysalis.output.cmds
    output:
        directory("results/co-assembly/trinity/{assembly}/trinity.multi-stage/parallel_jobs")
    log:
        "results/co-assembly/trinity/{assembly}/trinity.multi-stage/log_split"
    shell:
        """
        mkdir -p {output}
        split -n l/999 -e -d {input} {output}/job_ &>> {log}
        sed -i 's/\r$//' {output}/job_* &>> {log}
        """

rule trinity_co_butterfly_parallel:
    """
    Run Trinity Butterfly commands (split for parallelization).
    """
    input:
        "results/co-assembly/trinity/{assembly}/trinity.multi-stage/parallel_jobs/job_{job_index}"
    output:
        "results/co-assembly/trinity/{assembly}/trinity.multi-stage/parallel_jobs/completed_{job_index}"
    log:
        "results/co-assembly/trinity/{assembly}/trinity.multi-stage/parallel_jobs/log_parallel_{job_index}"
    params:
        tmpdir=tmpdir
    conda: "../../envs/trinity.yaml"
    container: "docker://trinityrnaseq/trinityrnaseq:2.15.2"
    shell:
        """
        bash {input} > {log}
        cp {input} {output}
        """

def trinity_co_completed_jobs(wildcards):
    import os
    from glob import glob
    parallel_dir = checkpoints.trinity_co_butterfly_split.get(**wildcards).output[0]
    job_ids = [os.path.basename(p).replace("job_", "") for p in glob(f"{parallel_dir}/job_*")]
    return [f"{parallel_dir}/completed_{job_index}" for job_index in job_ids]

rule trinity_co_butterfly_merge:
    """
    Merge outputs of parallel Butterfly jobs.
    """
    input:
        jobs=trinity_co_completed_jobs,
        cmds=rules.trinity_co_chrysalis.output.cmds
    output:
        "results/co-assembly/trinity/{assembly}/trinity.multi-stage/recursive_trinity.cmds.completed"
    log:
        "results/co-assembly/trinity/{assembly}/trinity.multi-stage/log_merge"
    shell:
        """
        cat {input.jobs} > {output} 2> {log}
        """

rule trinity_co_final:
    """
    Finalize the Trinity assembly.
    """
    input:
        cmds_completed=rules.trinity_co_butterfly_merge.output[0],
        R1=rules.trinity_normalize.output.R1,
        R2=rules.trinity_normalize.output.R2
    output:
        fa="results/co-assembly/trinity/{assembly}/final.fa",
        map="results/co-assembly/trinity/{assembly}/final.fa.gene_trans_map"
    log:
        "results/co-assembly/trinity/{assembly}/log_final"
    params:
        min_contig_len = config["min_contig_len"],
        max_mem = lambda wildcards, resources: int(resources.mem_mb /1000),
        outdir=lambda wildcards, input: os.path.dirname(input.cmds_completed)
    conda: "../../envs/trinity.yaml"
    container: "docker://trinityrnaseq/trinityrnaseq:2.15.2"
    shell:
        """
        Trinity \
            --max_memory {params.max_mem}G \
            --CPU {threads} \
            --no_bowtie --min_kmer_cov 2 \
            --output {params.outdir} \
            --left {input.R1} --right {input.R2} \
            --min_contig_length {params.min_contig_len} \
            --seqType fq \
            --no_normalize_reads \
            --full_cleanup \
            > {log} 2>&1
        mv {params.outdir}.Trinity.fasta {output.fa}
        mv {params.outdir}.Trinity.fasta.gene_trans_map {output.map}
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
        tmp_dir = "megahit/{assembly}",
        tmp_dir_base = "megahit"
    conda: "../../envs/megahit.yaml"
    shadow: "shallow"
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