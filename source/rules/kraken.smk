localrules: download_kraken_db, extract_kraken_reads, kraken_nohost

rule download_kraken_db:
    output:
        expand("resources/kraken/{f}.k2d", f = ["hash","opts","taxo"])
    params:
        hash = "https://zenodo.org/record/3259407/files/hash.k2d?download=1",
        info = "https://zenodo.org/record/3259407/files/kraken.info?download=1",
        opts = "https://zenodo.org/record/3259407/files/opts.k2d?download=1",
        taxo = "https://zenodo.org/record/3259407/files/taxo.k2d?download=1"
    shell:
        """
        curl -L -o resources/kraken/hash.k2d {params.hash}
        curl -L -o resources/kraken/kraken.info {params.info}
        curl -L -o resources/kraken/opts.k2d {params.opts}
        curl -L -o resources/kraken/taxo.k2d {params.taxo}
        """

rule preload_kraken_db:
    """
    Preload the Kraken database into memory
    """
    output:
        temp(touch("resources/kraken/{kraken_db}.preload"))
    input:
        expand("{kraken_db_dir}/{{kraken_db}}/{f}.k2d", f = ["hash","opts","taxo"], kraken_db_dir=config["kraken_db_dir"])
    params:
        db = lambda wildcards: os.path.join(config["kraken_db_dir"], wildcards.kraken_db)
    resources:
        mem_mb=lambda wildcards, input: max(1.15 * int(os.stat(input[0]).st_size / 1024**2), 888),
    group: "kraken2"
    shell:
        """
        mkdir -p /dev/shm/{wildcards.kraken_db}
        cp {input} /dev/shm/{wildcards.kraken_db}
        """

rule run_kraken:
    """
    Run Kraken2 on the preprocessed reads
    """
    output:
        "results/kraken/{kraken_db}/{sample_id}/{sample_id}.out.gz",
        "results/kraken/{kraken_db}/{sample_id}/{sample_id}.kreport",
        "results/kraken/{kraken_db}/{sample_id}/{sample_id}.unclassified_1.fastq.gz",
        "results/kraken/{kraken_db}/{sample_id}/{sample_id}.unclassified_2.fastq.gz",
    input:
        R1 = rules.sortmerna.output.R1,
        R2 = rules.sortmerna.output.R2,
        db=rules.preload_kraken_db.output
    log:
        "results/kraken/{kraken_db}/{sample_id}/{sample_id}.log"
    threads: 4
    params:
        db = lambda wildcards: wildcards.kraken_db,
        tmp = "{sample_id}.out",
        unc_out = lambda wildcards, output: os.path.join(os.path.dirname(output[0]), wildcards.sample_id+ ".unclassified#.fastq.gz"),
    conda: "../../envs/kraken.yaml"
    container: "docker://quay.io/biocontainers/kraken2:2.1.3--pl5321h077b44d_4"
    group: "kraken2"
    shadow: "shallow"
    shell:
        """
        kraken2 --db /dev/shm/{params.db} --memory-mapping --output {params.tmp} --report {output[1]} --gzip-compressed \
            --threads {threads} --unclassified-out {params.unc_out} --paired {input.R1} {input.R2} > {log} 2>&1
        gzip {params.tmp}
        mv {params.tmp}.gz {output[0]}
        """
    
rule kraken:
    output:
        temp(touch(expand("results/kraken/{kraken_db}.done", kraken_db=config["kraken_db"])))
    input:
        expand("results/kraken/{kraken_db}/{sample_id}/{sample_id}.out.gz", sample_id=samples.keys(), kraken_db=config["kraken_db"])
    params:
        kraken_db=config["kraken_db"]
    group: "kraken2"
    shell:
        """
        rm -rf /dev/shm/{params.kraken_db}
        """

rule extract_kraken_reads:
    """
    Extract reads for different taxonomic bins
    """
    output:
        R1="results/kraken/{kraken_db}/{sample_id}/taxbins/{taxname}_R1.fastq.gz",
        R2="results/kraken/{kraken_db}/{sample_id}/taxbins/{taxname}_R2.fastq.gz",
    input:
        kraken=rules.run_kraken.output[0],
        report=rules.run_kraken.output[1],
        R1=rules.sortmerna.output.R1,
        R2=rules.sortmerna.output.R2,
    log:
        "results/kraken/{kraken_db}/{sample_id}/taxbins/{taxname}.extract_reads.log"
    params:
        tmpdir="{taxname}.{sample_id}",
        taxids=lambda wildcards: " ".join([str(x) for x in config["taxmap"][wildcards.taxname]]),
        R1="{taxname}.{sample_id}_R1.fastq",
        R2="{taxname}.{sample_id}_R2.fastq",
    conda: "../../envs/kraken-tools.yaml" 
    shadow: "shallow"
    container: "docker://quay.io/biocontainers/krakentools:1.1--pyh5e36f6f_0"
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.kraken} > {params.tmpdir}/kraken.out
        extract_kraken_reads.py -k {params.tmpdir}/kraken.out -s1 {input.R1} -s2 {input.R2} \
            --fastq-output -r {input.report} --include-children -t {params.taxids} -o {params.R1} -o2 {params.R2} >{log} 2>&1
        gzip {params.R1} {params.R2}
        mv {params.R1}.gz {output.R1}
        mv {params.R2}.gz {output.R2}
        rm -rf {params.tmpdir}
        """

rule kraken_nohost:
    """
    Remove host reads from Kraken taxbins
    """
    output:
        R1="results/kraken/{kraken_db}/{sample_id}/taxbins/{taxname}_R1.nohost.fastq.gz",
        R2="results/kraken/{kraken_db}/{sample_id}/taxbins/{taxname}_R2.nohost.fastq.gz",
    input:
        R1=rules.extract_kraken_reads.output.R1,
        R2=rules.extract_kraken_reads.output.R2,
        R1nohost="results/filtered/{sample_id}/{sample_id}_R1.nohost.fastq.gz",
        R2nohost="results/filtered/{sample_id}/{sample_id}_R2.nohost.fastq.gz",
    log:
        "results/kraken/{kraken_db}/{sample_id}/taxbins/{taxname}.extract_reads.log"
    shell:
        """
        exec &> {log}
        seqkit grep -f <(seqkit seq -n -i {input.R1nohost}) {input.R1} -o {output.R1}
        seqkit grep -f <(seqkit seq -n -i {input.R2nohost}) {input.R2} -o {output.R2}
        """
