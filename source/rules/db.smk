localrules:
    ncbi_taxonomy,
    download_taxdump,
    download_taxmapper,
    download_eggnog,
    download_dbCAN,
    press_dbCAN,
    get_kegg_files,
    download_refseq_db,
    download_host,
    init_jgi,
    download_jgi_transcripts,
    download_jgi_proteins,
    concat_transcripts,
    concat_proteins,
    mmseqs_createseqdb,
    prepare_diamond_JGI,
    mmseqs_createtaxdb

wildcard_constraints:
    db = config["mmseqs_db"],

## NCBI Taxonomy ##
rule ncbi_taxonomy:
    output:
        "resources/taxonomy/taxonomy.sqlite"
    run:
        from ete3 import NCBITaxa
        shell("touch {output[0]}")
        ncbi_taxa = NCBITaxa(output[0])

rule build_blobdb:
    input:
        names = "resources/taxonomy/names.dmp",
        nodes = "resources/taxonomy/nodes.dmp"
    output:
        flag = "resources/taxonomy/blobdb.done"
    conda: "../../envs/blobtools.yaml"
    shell:
        """
        blobtools nodesdb --nodes {input.nodes} --names {input.names}
        """

## Taxmapper ##
rule download_taxmapper:
    output:
        "resources/taxmapper/databases/taxonomy/meta_database.fasta",
        "resources/taxmapper/databases/taxonomy/meta_database.fa.gz"
    params:
        dbdir = "resources/taxmapper",
        tmpdir = "$TMPDIR/taxmapper"
    shell:
        """
        mkdir -p {params.tmpdir} {params.dbdir}
        wget -O {params.tmpdir}/supplement.zip https://bitbucket.org/dbeisser/taxmapper_supplement/get/supplement.zip
        unzip -o -d {params.tmpdir} {params.tmpdir}/supplement.zip
        rm -rf {params.dbdir}/*
        mv -f {params.tmpdir}/dbeisser-taxmapper_supplement-854f5f60158a/* {params.dbdir}/
        gunzip -c {output[1]} > {output[0]}
        rm -r {params.tmpdir}
        """

rule format_taxmapper:
    input:
        "resources/taxmapper/databases/taxonomy/meta_database.fasta"
    output:
        "resources/taxmapper/databases/taxonomy/meta_database.db"
    conda:
        "../../envs/taxmapper.yaml"
    resources:
        mem_mb=4096
    shell:
        """
        prerapsearch -d {input[0]} -n {output[0]}
        """

rule download_eggnog:
    output:
        expand("resources/eggnog/{f}", f = ["eggnog.db","eggnog_proteins.dmnd"])
    params:
        data_dir = "resources/eggnog"
    conda: "../../envs/emapper.yaml"
    shell:
        """
        mkdir -p {params.data_dir}
        # Download eggnog.db
        download_eggnog_data.py --data_dir resources/eggnog -y
        """

## dbCAN ##
rule download_dbCAN:
    output:
        "resources/dbCAN/dbCAN-fam-HMMs.txt"
    shell:
        """
        curl -L -o {output[0]} http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V6.txt
        """

rule press_dbCAN:
    input:
        "resources/dbCAN/dbCAN-fam-HMMs.txt"
    output:
        expand("resources/dbCAN/dbCAN-fam-HMMs.txt.h3{suffix}", suffix = ["f","i","m","p"])
    conda: "../../envs/hmmer.yaml"
    shell:
        """
        hmmpress {input[0]}
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
    For each 'portal', download the filtered transcripts 
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

rule concat_transcripts:
    """
    Concatenates all non-zero transcripts downloaded
    """
    output:
        "resources/fungi/fungi_transcripts.fasta.gz"
    input:
        expand(rules.download_jgi_transcripts.output.transcripts, portal=genomes.index.tolist())
    log:
        "resources/fungi/concat_transcripts.log"
    params:
        tmpfile = "$TMPDIR/fungi_transcripts.fasta.gz"
    shell:
        """
        for f in {input};
        do
            if [-s $f];
            then
                cat $f >> {params.tmpfile}
            fi
        done
        mv {params.tmpfile} {output} 
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
    run:
        from Bio.SeqIO import parse
        import gzip as gz
        with gz.open(output.faa, "wt") as fhout:
            for f in input.faa:
                with gz.open(f, 'rt') as fhin:
                    for record in parse(fhin, "fasta"):
                        fhout.write(f">{record.id}\n{record.seq}\n")
        with open(output.mapping, "w") as fhout:
            for f in input.mapping:
                with open(f) as fhin:
                    for line in fhin:
                        fhout.write(line)

rule mmseqs_extract_fungalDB:
    """
    Uses the filtertaxseqdb command to extract fungal sequences from the official mmseqs2 database
    """
    input:
        db=os.path.join(config["mmseqs_db_dir"],"{db}"),
    output:
        db="resources/mmseqs2/fungi-{db}"
    log:
        "resources/mmseqs2/extract_fungal_{db}_mmseqsDB.log"
    threads: 4
    shell:
        """
        mmseqs filtertaxseqdb {input.db} {output.db} --taxon-list 4751 --threads {threads} > {log} 2>&1
        """

rule mmseqs_convert2fasta_fungalDB:
    """
    Converts the extracted mmseqs2 database to fasta format
    """
    input:
        db=rules.mmseqs_extract_fungalDB.output.db
    output:
        fasta="resources/mmseqs2/fungi-{db}.fasta"
    conda: "../../envs/mmseqs.yaml"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    shell:
        """
        mmseqs convert2fasta {input.db} {output.fasta}
        """

rule mmseqs_createseqdb:
    """
    Creates a sequence database for the concatenated proteins
    """
    output:
        db="resources/mmseqs2/combined-fungi-{db}",
        db_files=expand("resources/mmseqs2/combined-fungi-{{db}}{ext}", 
            ext=[".dbtype","_h","_h.dbtype","_h.index",".index",".lookup",".source"])
    input:
        jgi_fasta=rules.concat_proteins.output.faa,
        mmseqs_fasta=rules.mmseqs_convert2fasta_fungalDB.output.fasta
    log:
        "resources/mmseqs2/create-fungi-{db}-seqdb.log"
    conda: "../../envs/mmseqs.yaml"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    shell:
        """
        mmseqs createdb {input.jgi_fasta} {input.mmseqs_fasta} {output.db} --dbtype 1 > {log} 2>&1
        """

rule mmseqs_create_taxidmap:
    """
    Create taxid mapping file for the official mmseqs2 database
    """
    input:
        mapfile=os.path.join(config["mmseqs_db_dir"], "{db}_mapping"),
        lookupfile=os.path.join(config["mmseqs_db_dir"], "{db}.lookup"),
    output:
        tsv="resources/mmseqs2/{db}.taxidmap.tsv"
    resources:
        mem_mb = 90000
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
        pd.merge(protmap, taxmap, left_index=True, right_index=True).to_csv(output.tsv, sep="\t", index=False)


rule download_taxdump:
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
        db_files=expand("resources/mmseqs2/combined-fungi-{{db}}{ext}",  ext=["_taxonomy","_mapping"])
    input:
        db=rules.mmseqs_createseqdb.output.db,
        mapfiles=[rules.mmseqs_create_taxidmap.output.tsv, rules.concat_proteins.output.mapping],
        taxdump=rules.download_taxdump.output
    log:
        "resources/mmseqs2/createtaxdb-combined-fungi-{db}.log"
    conda: "../../envs/mmseqs.yaml"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    params:
        taxdump = lambda wc, input: os.path.dirname(input.taxdump[0]),
        tmpdir = os.environ.get("TMPDIR", "scratch"),
        mapfile = lambda wc: f"resources/mmseqs2/combined-fungi-{wc.db}.mapping.tsv"
    #shadow: "minimal"
    threads: 100
    shell:
        """
        cat {input.mapfiles} > {params.mapfile}
        mmseqs createtaxdb {input.db} {params.tmpdir} --ncbi-tax-dump {params.taxdump} --tax-mapping-file {params.mapfile} --threads {threads} > {log} 2>&1
        """

rule cluster_transcripts:
    """
    Clusters transcripts using vsearch
    """
    output:
        "resources/fungi/fungi_transcripts.clustered.fasta"
    input:
        rules.concat_transcripts.output
    log:
        "resources/fungi/cluster_transcripts.log"
    params:
        tmpdir = "$TMPDIR/fungi",
        id = 1.0
    conda:
        "../../envs/vsearch.yaml"
    threads: 10
    shell:
        """
        mkdir -p {params.tmpdir}
        pigz -p {threads} -c -d {input} > {params.tmpdir}/fungi_transcripts.fasta
        vsearch --cluster_fast {params.tmpdir}/fungi_transcripts.fasta \
            --centroids {params.tmpdir}/fungi_transcripts.clustered.fasta --id {params.id} \
            --log {log} --threads {threads}
        mv {params.tmpdir}/fungi_transcripts.clustered.fasta {output}
        rm -r {params.tmpdir}
        """

## Host data ##
rule download_host:
    """
    Downloads host sequences for filtering 
    """
    output:
        "resources/host/host.{suff}"
    log:
        "resources/host/host.{suff}.download.log"
    params:
        url = lambda wildcards: config["host_"+wildcards.suff+"_url"]
    shell:
        """
        exec &>{log}
        curl -L -o {output[0]}.gz {params.url}
        gunzip -c {output[0]}.gz > {output[0]}
        rm {output[0]}.gz
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

rule prepare_diamond_JGI:
    input:
        zip = "resources/JGI/fungi_proteins_download.zip",
        list = "resources/JGI/fungi_proteins_download.protfile1",
        gzip = "resources/taxmapper/databases/taxonomy/meta_database.fa.gz",
        extra = "resources/refseq/hygrophorus/HygMG78.faa",
        extra_id = "resources/refseq/hygrophorus/HygMG78.idmap.gz",
        tmmap = "resources/taxmapper/taxmapper_taxonomy_taxids.tsv",
        jgimap = "resources/JGI/taxonomy_taxids_newids.tsv",
        jgiidmap = "resources/JGI/JGI_idmap.tsv"
    output:
        "resources/diamond/taxonmap.gz",
        "resources/diamond/jgi.map.gz",
        "resources/diamond/tm.map.gz",
        "resources/diamond/refseq.map.gz",
        "resources/diamond/fasta.gz"
    params:
        src = "source/utils/reformat_prot_fasta.py"
    shell:
        """
        # Create temporary directory
        mkdir -p $TMPDIR/fasta
        rm -rf $TMPDIR/fasta/*
        # Add taxmapper proteins
        gunzip -c {input.gzip} | python {params.src} {input.tmmap} --accessionlength 14 > $TMPDIR/fasta/db.fasta 2>$TMPDIR/fasta/tm.map
        # Add JGI proteins
        files=$(cat {input.list})
        for f in $files;
        do
            genome_id=$(echo $f | cut -f1 -d '/')
            new_id=$(grep -w $genome_id {input.jgiidmap} | cut -f1)
            unzip -p {input.zip} $f | gunzip -c | python {params.src} {input.jgimap} --accessionlength 14 --genome_id $new_id >> $TMPDIR/fasta/db.fasta 2>>$TMPDIR/fasta/jgi.map
        done
        # Add extra Hygro proteins
        cat {input.extra} >> $TMPDIR/fasta/db.fasta
        gunzip -c {input.extra_id} >> $TMPDIR/fasta/refseq.map
        # Make taxonmap file
        echo -e "accession\taccession.version\ttaxid\tgi" > $TMPDIR/fasta/taxonmap
        cut -f1 $TMPDIR/fasta/tm.map $TMPDIR/fasta/jgi.map $TMPDIR/fasta/refseq.map > $TMPDIR/fasta/taxonmap.1
        cut -f1,3 $TMPDIR/fasta/tm.map $TMPDIR/fasta/jgi.map > $TMPDIR/fasta/taxonmap.2
        paste $TMPDIR/fasta/taxonmap.1 $TMPDIR/fasta/taxonmap.2 $TMPDIR/fasta/taxonmap.1 >> $TMPDIR/fasta/taxonmap
        gzip $TMPDIR/fasta/taxonmap
        gzip $TMPDIR/fasta/jgi.map
        gzip $TMPDIR/fasta/tm.map
        gzip $TMPDIR/fasta/refseq.map
        gzip $TMPDIR/fasta/db.fasta
        mv $TMPDIR/fasta/taxonmap.gz {output[0]}
        mv $TMPDIR/fasta/jgi.map.gz {output[1]}
        mv $TMPDIR/fasta/tm.map.gz {output[2]}
        mv $TMPDIR/fasta/refseq.map.gz {output[3]}
        mv $TMPDIR/fasta/db.fasta.gz {output[4]}
        rm -r $TMPDIR/fasta
        """

rule download_diamond_files:
    output:
        taxonmap = "resources/diamond/taxonmap.gz",
        fasta = "resources/diamond/fasta.gz",
        nodes = "resources/diamond/nodes.dmp"
    params:
        diamond_taxonmap_url = config["diamond_taxonmap_url"],
        diamond_fasta_url = config["diamond_fasta_url"],
        diamond_nodes_url = config["diamond_nodes_url"]
    shell:
        """
        curl -L -o {output.taxonmap} {params.diamond_taxonmap_url}
        curl -L -o {output.fasta} {params.diamond_fasta_url}
        curl -L -o {output.nodes}.gz {params.diamond_nodes_url}
        gunzip {output.nodes}.gz
        """
ruleorder: download_diamond_files > prepare_diamond_JGI

rule build_diamond_JGI:
    input:
        taxonmap = "resources/diamond/taxonmap.gz",
        fasta = "resources/diamond/fasta.gz",
        nodes = "resources/diamond/nodes.dmp"
    output:
        db = "resources/diamond/diamond.dmnd"
    conda: "../../envs/diamond.yaml"
    params:
        tmpdir = "$TMPDIR/diamond"
    threads: 4
    resources:
        runtime = 60*5
    shell:
        """
        # Create temporary directory
        mkdir -p $TMPDIR
        # Create the database
        zcat {input.fasta} | diamond makedb -d $TMPDIR/diamond -p {threads} \
         --taxonmap {input.taxonmap} --taxonnodes {input.nodes}
        # Move to output
        mv $TMPDIR/diamond.dmnd {output.db}
        """

rule build_diamond_JGI_legacy:
    input:
        taxonmap = "resources/diamond/taxonmap.gz",
        fasta = "resources/diamond/fasta.gz",
        nodes = "resources/diamond/nodes.dmp"
    output:
        db = "resources/diamond_legacy/diamond.dmnd"
    conda: "../../envs/contigtax.yaml"
    params:
        tmpdir = "$TMPDIR/diamond"
    threads: 4
    resources:
        runtime = 60*5
    shell:
        """
        # Create temporary directory
        mkdir -p $TMPDIR
        # Create the database
        zcat {input.fasta} | diamond makedb -d $TMPDIR/diamond -p {threads} \
         --taxonmap {input.taxonmap} --taxonnodes {input.nodes}
        # Move to output
        mv $TMPDIR/diamond.dmnd {output.db}
        """

## KEGG info ##
rule get_kegg_files:
    output:
        expand("resources/kegg/{f}",
            f = ["kegg_ec2pathways.tsv","kegg_ko2ec.tsv",
                 "kegg_ko2pathways.tsv","kegg_kos.tsv","kegg_modules.tsv","kegg_pathways.tsv"])
    params:
        dldir = "resources/kegg",
        src = "source/utils/eggnog-parser.py"
    shell:
        """
        python {params.src} download {params.dldir}
        """
