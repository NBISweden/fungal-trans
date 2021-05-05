localrules: download_kraken_db

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
        "results/kraken/{sample_id}.kreport"
    threads: 20
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*1
    params:
        db = "resources/kraken",
        tmp = "$TMPDIR/{sample_id}.out"
    conda: "../../envs/kraken.yaml"
    shell:
        """
        kraken2 --db {params.db} --output {params.tmp} --report {output[1]} --gzip-compressed \
        --threads {threads} --paired {input.R1} {input.R2}
        gzip {params.tmp}
        mv {params.tmp}.gz {output[0]}
        """
