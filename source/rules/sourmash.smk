###############################################
## Sourmash analysis of reads and assemblies ##
###############################################
localrules:
    sourmash_compute_bowtie,
    sourmash_compute_taxmapper,
    sourmash_compute_filtered,
    sourmash_compute_assembly

rule sourmash_compute_bowtie:
    input:
        R1="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.nospruce.fastq.gz",
        R2="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.nospruce.fastq.gz"
    output:
        "results/sourmash/{sample_id}/{sample_id}.bowtie.sig"
    conda: "../../envs/sourmash.yaml"
    shell:
        """
        sourmash compute --scaled 100 --merge {wildcards.sample_id} -o {output[0]} -k 31 {input.R1} {input.R2}
        """

rule sourmash_compute_taxmapper:
    input:
        R1 = "results/taxmapper/{sample_id}/{sample_id}_R1.cut.trim.filtered.fastq.gz",
        R2 = "results/taxmapper/{sample_id}/{sample_id}_R2.cut.trim.filtered.fastq.gz"
    output:
        "results/sourmash/{sample_id}/{sample_id}.taxmapper.sig"
    conda: "../../envs/sourmash.yaml"
    shell:
        """
        sourmash compute --scaled 100 --merge {wildcards.sample_id} -o {output[0]} -k 31 {input.R1} {input.R2}
        """

rule sourmash_compute_filtered:
    input:
        R1 = "results/filtered/{sample_id}/{sample_id}_R1.filtered.union.fastq.gz",
        R2 = "results/filtered/{sample_id}/{sample_id}_R2.filtered.union.fastq.gz"
    output:
        "results/sourmash/{sample_id}/{sample_id}.filtered.sig"
    conda: "../../envs/sourmash.yaml"
    shell:
        """
        sourmash compute --scaled 100 --merge {wildcards.sample_id} -o {output[0]} -k 31 {input.R1} {input.R2}
        """

rule sourmash_compute_assembly:
    input:
        fa = "results/assembly/{assembler}/{source}/{sample_id}/final.fa"
    output:
        sig = "results/assembly/{assembler}/{source}/{sample_id}/final.sig"
    conda: "../../envs/sourmash.yaml"
    shell:
        """
        sourmash compute --scaled 100 --merge {wildcards.assembler}.{wildcards.sample_id} -o {output.sig} -k 31 {input.fa}
        """

rule sourmash_compare:
    input:
        expand("results/sourmash/{sample_id}/{sample_id}.{{source}}.sig", sample_id = samples.keys())
    output:
        "results/report/sourmash/{source}_sample.dist",
        "results/report/sourmash/{source}_sample.dist.csv",
        "results/report/sourmash/{source}_sample.dist.labels.txt"
    conda: "../../envs/sourmash.yaml"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*10
    shell:
        """
        sourmash compare {input} --csv {output[1]} -o {output[0]}
        """

rule sourmash_compare_assembly:
    input:
        expand("results/assembly/{{assembler}}/{{source}}/{sample_id}/final.sig", sample_id=samples.keys())
    output:
        "results/report/sourmash/{source}_assembly_{assembler}.dist",
        "results/report/sourmash/{source}_assembly_{assembler}.dist.csv",
        "results/report/sourmash/{source}_assembly_{assembler}.dist.labels.txt"
    conda: "../../envs/sourmash.yaml"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*10
    shell:
        """
        sourmash compare {input} --csv {output[1]} -o {output[0]}
        """
