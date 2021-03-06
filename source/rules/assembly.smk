localrules: assembly_stats
##############
## Assembly ##
##############

def assembly_input(wildcards):
    suffices = {'unfiltered': 'cut.trim.fastq.gz',
                'filtered': 'filtered.union.fastq.gz',
                'taxmapper': 'cut.trim.filtered.fastq.gz',
                'bowtie2': 'fungi.nohost.fastq.gz'}
    R1 = 'results/{source}/{sample}/{sample}_R1.{suff}'.format(source=wildcards.filter_source,
        sample=wildcards.sample_id,suff=suffices[wildcards.filter_source])
    R2 = 'results/{source}/{sample}/{sample}_R2.{suff}'.format(source=wildcards.filter_source,
        sample=wildcards.sample_id,suff=suffices[wildcards.filter_source])
    return [R1, R2]

rule transabyss:
    input:
        assembly_input
    output:
        fa = "results/assembly/transabyss/{filter_source}/{sample_id}/{sample_id}.{k}-final.fa",
    log: "results/assembly/transabyss/{filter_source}/{sample_id}/{k}.log"
    conda:
        "../../envs/transabyss.yaml"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        tmpdir = "$TMPDIR/{sample_id}.{k}.transabyss",
        R1 = "$TMPDIR/{sample_id}.{k}.transabyss/R1.fastq",
        R2 = "$TMPDIR/{sample_id}.{k}.transabyss/R2.fastq"
    threads: 16
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 150
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input[0]} > {params.R1}
        gunzip -c {input[1]} > {params.R2}
        transabyss --pe {params.R1} {params.R2} -k {wildcards.k} \
            --outdir {params.tmpdir} --name {wildcards.sample_id}.{wildcards.k} \
            --threads {threads} >{log} 2>&1
        mv {params.tmpdir}/* {params.outdir}
        rm -rf {params.tmpdir}
        """

rule transabyss_merge:
    input:
        assembly_input,
        expand("results/assembly/transabyss/{{filter_source}}/{{sample_id}}/{{sample_id}}.{k}-final.fa",
            k = config["transabyss_kmers"])
    output:
        fa = "results/assembly/transabyss/{filter_source}/{sample_id}/final.fa",
        R1 = "results/assembly/transabyss/{filter_source}/{sample_id}/{sample_id}_R1.fastq.gz",
        R2= "results/assembly/transabyss/{filter_source}/{sample_id}/{sample_id}_R2.fastq.gz"
    log:
        "results/logs/transabyss/{sample_id}.{filter_source}.merge.log"
    conda:
        "../../envs/transabyss.yaml"
    params:
        mink = min(config["transabyss_kmers"]),
        maxk = max(config["transabyss_kmers"]),
        i = lambda wildcards, input: sorted(input[2:]),
        prefix = ["k{x}.".format(x=x) for x in sorted(config["transabyss_kmers"])],
        tmpout = "$TMPDIR/{sample_id}.ta.merged.fa"
    threads: 16
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 150
    shell:
        """
        rm -rf {params.tmpout}
        wd=$(pwd)
        # Link input into output directory
        ln -s $wd/{input[0]} $wd/{output.R1}
        ln -s $wd/{input[1]} $wd/{output.R2}
        transabyss-merge {params.i} --mink {params.mink} --maxk {params.maxk} \
            --out {params.tmpout} --threads {threads} --force --prefixes {params.prefix} > {log} 2>&1
        mv {params.tmpout} {output.fa}
        """


rule trinity:
    input:
        assembly_input
    output:
        fa="results/assembly/trinity/{filter_source}/{sample_id}/final.fa",
        R1="results/assembly/trinity/{filter_source}/{sample_id}/{sample_id}_R1.fastq.gz",
        R2="results/assembly/trinity/{filter_source}/{sample_id}/{sample_id}_R2.fastq.gz"
    log: "results/assembly/trinity/{filter_source}/{sample_id}/log"
    params:
        min_contig_len = config["min_contig_len"],
        tmpdir = "$TMPDIR/{sample_id}.trinity",
        outdir = "results/assembly/trinity/{filter_source}/{sample_id}",
        cpumem = config["mem_per_cpu"],
        out_base= lambda wildcards,output: os.path.basename(output.fa)
    threads: 10
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 10
    conda:
        "../../envs/trinity.yaml"
    shell:
        """
        rm -rf {params.tmpdir}
        rm -rf {params.outdir}/*
        mkdir -p {params.tmpdir}
        wd=$(pwd)
        max_mem=$(({params.cpumem} * {threads}))
        # Link input into output directory
        ln -s $wd/{input[0]} $wd/{output.R1}
        ln -s $wd/{input[1]} $wd/{output.R2}
        Trinity --CPU {threads} --min_contig_length {params.min_contig_len} \
            --output {params.tmpdir} --left {input[0]} --right {input[1]} \
            --seqType fq --max_memory ${{max_mem}}G > {log} 2>&1
        mv {params.tmpdir}/* {params.outdir}/
        cd {params.outdir}
        ln -s Trinity.fasta {params.out_base}
        rm -rf {params.tmpdir}
        """

rule megahit:
    input:
        assembly_input
    output:
        fa = "results/assembly/megahit/{filter_source}/{sample_id}/final.fa",
        R1 = "results/assembly/megahit/{filter_source}/{sample_id}/{sample_id}_R1.fastq.gz",
        R2 = "results/assembly/megahit/{filter_source}/{sample_id}/{sample_id}_R2.fastq.gz"
    log: "results/assembly/megahit/{filter_source}/{sample_id}/log"
    params:
        min_contig_len = config["min_contig_len"],
        out_dir = "results/assembly/megahit/{filter_source}/{sample_id}",
        tmp_dir = "$TMPDIR/megahit/{sample_id}",
        tmp_dir_base = "$TMPDIR/megahit"
    conda: "../../envs/megahit.yaml"
    threads: 10
    resources:
        runtime = lambda wildcards, attempt: attempt*attempt*30
    shell:
        """
        wd=$(pwd)
        # Link input into output directory
        ln -s $wd/{input[0]} $wd/{output.R1}
        ln -s $wd/{input[1]} $wd/{output.R2}
        mkdir -p {params.tmp_dir_base}
        megahit -1 {input[0]} -2 {input[1]} --prune-level 3 \
            --min-contig-len {params.min_contig_len} -o {params.tmp_dir} \
            -t {threads} > {log} 2>&1
        mv {params.tmp_dir}/final.contigs.fa {output.fa}
        mv {params.tmp_dir}/opt* {params.out_dir}
        rm -rf {params.tmp_dir}
        """

rule fastuniq_for_assembly:
    input:
        R1 = lambda wildcards: assemblies[wildcards.assembly]["R1"],
        R2 = lambda wildcards: assemblies[wildcards.assembly]["R2"]
    output:
        "results/fastuniq/{assembly}/R1.fastuniq.gz",
        "results/fastuniq/{assembly}/R2.fastuniq.gz"
    params:
        filelist = os.path.expandvars("$TMPDIR/fastuniq/{assembly}/filelist"),
        tmpdir = os.path.expandvars("$TMPDIR/fastuniq/{assembly}"),
        R1 = os.path.expandvars("$TMPDIR/fastuniq/{assembly}/R1.fastuniq"),
        R2 = os.path.expandvars("$TMPDIR/fastuniq/{assembly}/R2.fastuniq")
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*10,
        mem_mb = lambda wildcards, attempt: attempt**2*128000
    threads: 20
    run:
        shell("mkdir -p {params.tmpdir}")
        files = []
        # Unzip and create filelist for fastuniq
        for i, r1 in enumerate(input.R1):
            basename1 = os.path.basename(r1).rstrip(".gz")
            shell("gunzip -c {r1} > {params.tmpdir}/{basename1}")
            r2 = input.R2[i]
            basename2 = os.path.basename(r2).rstrip(".gz")
            shell("gunzip -c {r2} > {params.tmpdir}/{basename2}")
            files+=[os.path.join(params.tmpdir,basename1),os.path.join(params.tmpdir, basename2)]
        # Write the filelist
        with open(params.filelist, 'w') as fhout:
            fhout.write("{}".format("\n".join(files)))
        shell("fastuniq -i {params.filelist} -o {params.R1} -p {params.R2} -c 1")
        shell("gzip {params.R1}")
        shell("gzip {params.R2}")
        shell("mv {params.R1}.gz {output[0]}")
        shell("mv {params.R2}.gz {output[1]}")
        shell("rm -r {params.tmpdir}")

rule transabyss_co:
    input:
        R1="results/fastuniq/{assembly}/R1.fastuniq.gz",
        R2="results/fastuniq/{assembly}/R2.fastuniq.gz"
    output:
        "results/co-assembly/transabyss/{assembly}/{assembly}.{k}-final.fa",
    log: "results/co-assembly/transabyss/{assembly}/{k}.log"
    conda:
        "../../envs/transabyss.yaml"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        tmpdir="$TMPDIR/{assembly}.{k}.transabyss",
        R1="$TMPDIR/{assembly}.{k}.transabyss/R1.fastq",
        R2="$TMPDIR/{assembly}.{k}.transabyss/R2.fastq"
    threads: 16
    resources:
        runtime = lambda wildcards,attempt: attempt ** 2 * 60 * 150
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.R1} > {params.R1}
        gunzip -c {input.R2} > {params.R2}
        transabyss --pe {params.R1} {params.R2} -k {wildcards.k} \
            --outdir {params.tmpdir} --name {wildcards.assembly}.{wildcards.k} \
            --threads {threads} >{log} 2>&1
        mv {params.tmpdir}/* {params.outdir}
        rm -rf {params.tmpdir}
        """

rule transabyss_merge_co:
    input:
        expand("results/co-assembly/transabyss/{{assembly}}/{{assembly}}.{k}-final.fa",
            k = config["transabyss_kmers"])
    output:
        fa = "results/co-assembly/transabyss/{assembly}/final.fa"
    log: "results/co-assembly/transabyss/{assembly}/merge.log"
    conda:
        "../../envs/transabyss.yaml"
    params:
        mink = min(config["transabyss_kmers"]),
        maxk=max(config["transabyss_kmers"]),
        i=lambda wildcards, input: sorted(input),
        prefix=["k{x}.".format(x=x) for x in
                sorted(config["transabyss_kmers"])],
        tmpout="$TMPDIR/{assembly}.ta.merged.fa"
    threads: 16
    resources:
        runtime = lambda wildcards,attempt: attempt ** 2 * 60 * 10
    shell:
        """
        rm -rf {params.tmpout}
        wd=$(pwd)
        transabyss-merge {params.i} --mink {params.mink} --maxk {params.maxk} \
            --out {params.tmpout} --threads {threads} --force --prefixes {params.prefix} > {log} 2>&1
        mv {params.tmpout} {output.fa}
        """

rule trinity_co:
    input:
        R1="results/fastuniq/{assembly}/R1.fastuniq.gz",
        R2="results/fastuniq/{assembly}/R2.fastuniq.gz"
    output:
        fa = "results/co-assembly/trinity/{assembly}/final.fa"
    log: "results/co-assembly/trinity/{assembly}/log"
    params:
        min_contig_len = config["min_contig_len"],
        tmpdir="$TMPDIR/{assembly}.co.trinity",
        outdir="results/co-assembly/trinity/{assembly}",
        out_base = lambda wildcards, output: os.path.basename(output.fa),
        cpumem = config["mem_per_cpu"]
    threads: 10
    resources:
        runtime = lambda wildcards,attempt: attempt ** 2 * 60 * 240
    conda:
        "../../envs/trinity.yaml"
    shell:
        """
        rm -rf {params.outdir}/*
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}
        wd=$(pwd)
        max_mem=$(({params.cpumem} * {threads}))
        Trinity --CPU {threads} --min_contig_length {params.min_contig_len} \
            --output {params.tmpdir} --left {input.R1} --right {input.R2} \
            --seqType fq --max_memory ${{max_mem}}G > {log} 2>&1
        mv {params.tmpdir}/* {params.outdir}/
        cd {params.outdir}
        ln -s Trinity.fasta {params.out_base}
        rm -rf {params.tmpdir}
        """

rule megahit_co:
    input:
        R1 = "results/fastuniq/{assembly}/R1.fastuniq.gz",
        R2 = "results/fastuniq/{assembly}/R2.fastuniq.gz"
    output:
        fa = "results/co-assembly/megahit/{assembly}/final.fa"
    log: "results/co-assembly/megahit/{assembly}/log"
    params:
        min_contig_len = config["min_contig_len"],
        out_dir = "results/co-assembly/megahit/{assembly}",
        tmp_dir = "$TMPDIR/megahit/{assembly}",
        tmp_dir_base = "$TMPDIR/megahit"
    conda: "../../envs/megahit.yaml"
    threads: 20
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*240
    shell:
        """
        wd=$(pwd)
        mkdir -p {params.tmp_dir_base}
        megahit -1 {input.R1} -2 {input.R2} --prune-level 3 \
            --min-contig-len {params.min_contig_len} -o {params.tmp_dir} \
            -t {threads} > {log} 2>&1
        mv {params.tmp_dir}/final.contigs.fa {output.fa}
        mv {params.tmp_dir}/opts.txt {params.out_dir}
        rm -r {params.tmp_dir}
        """

rule assembly_stats:
    input:
        expand("results/assembly/{{assembler}}/{{source}}/{sample_id}/final.fa",
            sample_id = samples.keys())
    output:
        "results/report/assembly/{source}_{assembler}_stats.tsv",
        "results/report/assembly/{source}_{assembler}_size_dist.tsv"
    run:
        names = [x.split("/")[-2] for x in input]
        shell("python source/utils/assembly_stats.py -i {input} -n {names} --size-dist-file {output[1]} > {output[0]}")

rule co_assembly_stats:
    input:
        "results/co-assembly/{assembler}/{assembly}/final.fa"
    output:
        "results/report/co-assembly/{assembler}.{assembly}_assembly_stats.tsv",
        "results/report/co-assembly/{assembler}.{assembly}_assembly_size_dist.tsv"
    shell:
        """
        python source/utils/assembly_stats.py -i {input[0]} -n {wildcards.assembler}.{wildcards.assembly} --size-dist-file {output[1]} > {output[0]}
        """