localrules: download_kraken_db, extract_kraken_reads

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

rule run_kraken:
    input:
        R1 = "results/preprocess/{sample_id}_R1.cut.trim.fastq.gz",
        R2 = "results/preprocess/{sample_id}_R2.cut.trim.fastq.gz",
        db = expand("resources/kraken/{f}.k2d", f = ["hash","opts","taxo"])
    output:
        "results/kraken/{sample_id}.out.gz",
        "results/kraken/{sample_id}.kreport",
        "results/kraken/{sample_id}.unclassified_1.fastq.gz",
        "results/kraken/{sample_id}.unclassified_2.fastq.gz",
    log:
        "results/kraken/{sample_id}.log"
    threads: 20
    params:
        db = "resources/kraken",
        tmp = "$TMPDIR/{sample_id}.out",
        unc_out = "results/kraken/{sample_id}.unclassified#.fastq.gz",
    conda: "../../envs/kraken.yaml"
    shell:
        """
        kraken2 --db {params.db} --output {params.tmp} --report {output[1]} --gzip-compressed \
        --threads {threads} --unclassified-out {params.unc_out} --paired {input.R1} {input.R2} > {log} 2>&1
        gzip {params.tmp}
        mv {params.tmp}.gz {output[0]}
        """

rule extract_kraken_reads:
    """
    Extract reads for different taxonomic bins
    """
    output:
        R1="results/taxbins/{taxname}/{sample_id}_R1.fastq.gz",
        R2="results/taxbins/{taxname}/{sample_id}_R2.fastq.gz",
    input:
        kraken=rules.run_kraken.output[0],
        report=rules.run_kraken.output[1],
        R1="results/preprocess/{sample_id}_R1.cut.trim.fastq.gz",
        R2="results/preprocess/{sample_id}_R2.cut.trim.fastq.gz",
    log:
        "results/taxbins/{taxname}/extract_kraken_reads.{sample_id}.log"
    params:
        tmpdir="$TMPDIR/{taxname}.{sample_id}",
        taxid=lambda wildcards: config["taxmap"][wildcards.taxname]
    conda: "../../envs/kraken-tools.yaml" 
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.kraken} > {params.tmpdir}/kraken.out
        extract_kraken_reads.py -k {params.tmpdir}/kraken.out -s1 {input.R1} -s2 {input.R2} --fastq-output -r {input.report} --include-children -t {wildcards.taxid} -o {output.R1} -o2 {output.R2} >{log} 2>&1
        rm -rf {params.tmpdir}
        """