localrules:
    download_taxdump,
    download_taxmapper,
    download_eggnog,
    download_dbCAN,
    press_dbCAN,
    get_kegg_files,
    download_refseq_db,
    init_jgi,
    download_jgi_transcripts,
    download_jgi_proteins,
    concat_proteins,
    mmseqs_filter_fungalDB

# eggnog

rule download_eggnog:
    """
    Downloads the eggnog database and proteins
    """
    output:
        expand("{emapper_db_dir}/{f}", emapper_db_dir=config["emaper_db_dir"], f = ["eggnog.db","eggnog_proteins.dmnd"])
    params:
        data_dir = "resources/eggnog"
    conda: "../../envs/emapper.yaml"
    shell:
        """
        mkdir -p {params.data_dir}
        # Download eggnog.db
        download_eggnog_data.py --data_dir resources/eggnog -y
        """

## Fungal transcripts ##
rule init_jgi:
    """
    Creates a cookie file for JGI authentication using user-defined credentials
    """
    output:
        cookies="resources/JGI/cookies"
    log:
        "resources/JGI/init_jgi.log"
    shadow: "minimal"
    params:
        user_name = config["jgi_user"],
        password = config["jgi_password"]
    retries: 3
    shell:
        """
        curl 'https://signon.jgi.doe.gov/signon/create' \
            --data-urlencode 'login={params.user_name}' \
            --data-urlencode 'password={params.password}' -c {output.cookies} > /dev/null 2>{log}
        """

rule download_jgi_transcripts:
    """
    For each 'portal', download a transcripts file
    """
    output:
        transcripts=temp(touch("resources/JGI/genomes/{portal}.transcripts.fna.gz")),
        #proteins=temp(touch("resources/JGI/genomes/{portal}.proteins.faa.gz"))
    input:
        cookies=rules.init_jgi.output.cookies,
    log:
        "resources/JGI/genomes/download_jgi_transcripts.{portal}.log"
    retries: 3
    shell:
        """
        python source/utils/download_jgi_transcripts.py -p {wildcards.portal} -c {input.cookies} -o {output.transcripts} 2>{log}
        """

rule filter_transcripts:
    """
    Filters transcripts to remove sequences shorter than 300 bp
    Also renames sequences so as to be compatible with downstream mapping
    """
    output:
        touch("resources/JGI/genomes/{portal}.transcripts.filt.fna.gz")
    input:
        rules.download_jgi_transcripts.output.transcripts
    log:
        "resources/JGI/genomes/filter_transcripts.{portal}.log"
    shell:
        """
        exec &>{log}
        if [ -s {input} ];
        then
            seqkit seq -m 300 {input} | seqkit replace -p .+ -r "{wildcards.portal}_{{nr}}" | gzip -c > {output}
        fi
        """

rule download_jgi_proteins:
    """
    For each 'portal', download the filtered proteins. Also output a mapping file
    of protein accession to taxid.
    """
    output:
        proteins=temp(touch("resources/JGI/genomes/{portal}.proteins.faa.gz")),
        mapping=temp(touch("resources/JGI/genomes/{portal}.mapping.tsv"))
    input:
        cookies=rules.init_jgi.output.cookies,
    log:
        "resources/JGI/genomes/download_jgi_proteins.{portal}.log"
    params:
        taxid = lambda wildcards: extra_genomes[wildcards.portal]
    retries: 3
    shell:
        """
        python source/utils/download_jgi_transcripts.py -p {wildcards.portal} -c {input.cookies} \
            --protein_out {output.proteins} --taxidmap {output.mapping} --taxid {params.taxid} 2>{log}
        """

rule concat_proteins:
    """
    Concatenates all non-zero proteins downloaded for genomes listed in 'extra_genomes' file.
    """
    output:
        faa="resources/extra_genomes/proteins.faa.gz",
        mapping="resources/extra_genomes/mapping.tsv"
    input:
        faa=expand(rules.download_jgi_proteins.output.proteins, portal=list(extra_genomes.keys())),
        mapping=expand(rules.download_jgi_proteins.output.mapping, portal=list(extra_genomes.keys()))
    log:
        "resources/extra_genomes/concat_proteins.log"
    params:
        min_len = 100
    run:
        from Bio.SeqIO import parse
        import gzip as gz
        with gz.open(output.faa, "wt") as fhout:
            for f in input.faa:
                with gz.open(f, 'rt') as fhin:
                    for record in parse(fhin, "fasta"):
                        if len(record.seq) >= params.min_len:
                            newid = (record.id).replace("|", ".")
                            fhout.write(f">{newid}\n{record.seq}\n")
        with open(output.mapping, "w") as fhout:
            for f in input.mapping:
                with open(f) as fhin:
                    for line in fhin:
                        seqid, taxid = line.strip().split("\t")
                        newid = seqid.replace("|", ".")
                        fhout.write(f"{newid}\t{taxid}\n")

rule mmseqs_extract_fungalDB:
    """
    Uses the filtertaxseqdb command to extract fungal sequences from the official mmseqs2 database
    """
    output:
        db="resources/mmseqs2/fungi-{mmseqs_db}"
    input:
        db=os.path.join(config["mmseqs_db_dir"],"{mmseqs_db}"),
    log:
        "resources/mmseqs2/extract_fungal_{mmseqs_db}_mmseqsDB.log"
    conda: "../../envs/mmseqs.yaml"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    threads: 1
    shell:
        """
        mmseqs filtertaxseqdb {input.db} {output.db} --taxon-list 4751 --threads {threads} > {log} 2>&1
        """

rule mmseqs_convert2fasta_fungalDB:
    """
    Converts the extracted mmseqs2 database to fasta format
    """
    output:
        fasta="resources/mmseqs2/fungi-{mmseqs_db}.fasta"
    input:
        db=rules.mmseqs_extract_fungalDB.output.db
    conda: "../../envs/mmseqs.yaml"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    threads: 1
    shell:
        """
        mmseqs convert2fasta {input.db} {output.fasta}
        """

rule mmseqs_filter_fungalDB:
    """
    Filters the extracted fungal database to remove sequences shorter than a minimum length
    """
    output:
        fasta="resources/mmseqs2/filtered-fungi-{mmseqs_db}.fasta"
    input:
        rules.mmseqs_convert2fasta_fungalDB.output.fasta
    log:
        "resources/mmseqs2/filter_fungi_{mmseqs_db}_mmseqsDB.log"
    params:
        min_len = 100
    shell:
        """
        seqkit seq -m {params.min_len} {input} > {output.fasta}
        """

rule mmseqs_createseqdb:
    """
    Creates a sequence database for the concatenated proteins
    """
    output:
        db="resources/mmseqs2/combined-fungi-{mmseqs_db}",
        db_files=expand("resources/mmseqs2/combined-fungi-{{mmseqs_db}}{ext}", 
            ext=[".dbtype","_h","_h.dbtype","_h.index",".index",".lookup",".source"])
    input:
        jgi_fasta=rules.concat_proteins.output.faa,
        mmseqs_fasta=rules.mmseqs_filter_fungalDB.output.fasta
    log:
        "resources/mmseqs2/create-fungi-{mmseqs_db}-seqdb.log"
    conda: "../../envs/mmseqs.yaml"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    threads: 1
    shell:
        """
        mmseqs createdb {input.jgi_fasta} {input.mmseqs_fasta} {output.db} --dbtype 1 > {log} 2>&1
        """

rule mmseqs_create_taxidmap:
    """
    Create taxid mapping file for the official mmseqs2 database
    """
    output:
        tsv="resources/mmseqs2/{mmseqs_db}.taxidmap.tsv"
    input:
        mapfile=os.path.join(config["mmseqs_db_dir"], "{mmseqs_db}_mapping"),
        lookupfile=os.path.join(config["mmseqs_db_dir"], "{mmseqs_db}.lookup"),
    threads: 1
    run:
        import pandas as pd
        protmap = pd.read_csv(input.lookupfile, 
                                sep="\t", 
                                index_col=0, 
                                header=None, 
                                usecols=[0,1], 
                                names=["id","accession"])
        taxmap = pd.read_csv(input.mapfile,
                                sep="\t", 
                                index_col=0, 
                                header=None, 
                                usecols=[0,1], 
                                names=["id", "taxid"])
        pd.merge(protmap, taxmap, left_index=True, right_index=True).to_csv(output.tsv, sep="\t", header = False, index=False)


rule download_taxdump:
    """
    Downloads the NCBI taxonomy dump
    """
    output:
        expand("resources/ncbi-taxdump/{pref}.dmp", pref=["nodes","delnodes","gencode","merged","names"]),
    log:
        "resources/ncbi-taxdump/download_taxdump.log"
    params:
        outdir = lambda wc, output: os.path.dirname(output[0])
    shell:
        """
        wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O - 2>{log} | tar -xz -C {params.outdir} 
        """

rule mmseqs_createtaxdb:
    """
    Creates a taxonomy database for the concatenated protein database
    """
    output:
        db_files=expand("resources/mmseqs2/combined-fungi-{{mmseqs_db}}{ext}",  ext=["_taxonomy","_mapping"])
    input:
        db=rules.mmseqs_createseqdb.output.db,
        mapfiles=[rules.mmseqs_create_taxidmap.output.tsv, rules.concat_proteins.output.mapping],
        taxdump=rules.download_taxdump.output
    log:
        "resources/mmseqs2/createtaxdb-combined-fungi-{mmseqs_db}.log"
    conda: "../../envs/mmseqs.yaml"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    params:
        taxdump = lambda wc, input: os.path.dirname(input.taxdump[0]),
        tmpdir = os.environ.get("TMPDIR", "scratch"),
        mapfile = lambda wc: f"resources/mmseqs2/combined-fungi-{wc.mmseqs_db}.mapping.tsv"
    threads: 1
    shell:
        """
        cat {input.mapfiles} > {params.mapfile}
        mmseqs createtaxdb {input.db} {params.tmpdir} --ncbi-tax-dump {params.taxdump} --tax-mapping-file {params.mapfile} --threads {threads} > {log} 2>&1
        """

## Protein databases ##
rule download_refseq_db:
    output:
        "resources/refseq/{domain}/refseq_{domain}.faa",
        "resources/refseq/{domain}/refseq_{domain}.version"
    params:
        dir = "resources/refseq/{domain}"
    shell:
        """
        mkdir -p {params.dir}
        wget -P {params.dir} ftp://ftp.ncbi.nlm.nih.gov/refseq/release/{wildcards.domain}/*protein.faa.gz
        gunzip -c {params.dir}/*.faa.gz > {output[0]}
        rm {params.dir}/*.faa.gz
        wget -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER > {output[1]}
        """

## KEGG info ##
rule get_kegg_files:
    """
    Downloads KEGG files
    """
    output:
        expand("resources/kegg/{f}",
            f = ["kegg_ec2pathways.tsv","kegg_ko2ec.tsv",
                 "kegg_ko2pathways.tsv","kegg_kos.tsv","kegg_modules.tsv","kegg_pathways.tsv"])
    params:
        dldir = "resources/kegg",
        src = workflow.source_path("../utils/eggnog-parser.py"),
    shell:
        """
        python {params.src} download {params.dldir}
        """
