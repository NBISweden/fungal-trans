localrules:
    collate_dbcan,
    collate_taxonomy,
    dbcan_parse,
    eggnog_merge_and_sum,
    normalize_featurecount,
    parse_eggnog,
    quantify_eggnog,
    sum_dbcan,
    sum_taxonomy,

##########################################
##          SINGLE ASSEMBLIES           ##
##########################################
# Below are rules for annotation of single assemblies

#####################
## FRAME SELECTION ##
#####################
rule transdecoder_longorfs:
    """
    Generate longest ORFs from assembly with transdecoder
    """
    input:
        fa = "results/assembly/{assembler}/{filter_source}/{sample_id}/final.fa"
    output:
        expand("results/annotation/{{assembler}}/{{filter_source}}/{{sample_id}}/transdecoder/final.fa.transdecoder_dir/longest_orfs.{suff}",
            suff = ["pep","gff3","cds"])
    log:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/transdecoder/longest_orfs.log"
    container: "docker://trinityrnaseq/transdecoder:5.7.1"
    params:
        output_dir = lambda wildcards: f"results/annotation/{wildcards.assembler}/{wildcards.filter_source}/{wildcards.sample_id}/transdecoder"
    resources:
        tasks= 1,
        runtime = 240
    shell:
        """
        TransDecoder.LongOrfs -t {input.fa} -O {params.output_dir} > {log} 2>&1
        """

rule mmseqs_createquerydb_longorfs:
    """
    Create mmseqs database from longest ORFs
    """
    input:
        query = rules.transdecoder_longorfs.output[0],
    output:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/transdecoder/longorfs-queryDB"
    log:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/transdecoder/longorfs-queryDB.log"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    conda: "../../envs/mmseqs.yaml"
    shell:
        """
        mmseqs createdb {input} {output} > {log} 2>&1
        """

rule mmseqs_firstpass_taxonomy:
    """
    Run mmseqs taxonomy on longest ORFs, also store the alignment results
    for use with transdecoder
    """
    input:
        query = rules.mmseqs_createquerydb_longorfs.output,
        target = os.path.join(config["mmseqs_db_dir"], "{td_db}"),
    output:
        tax = expand("results/annotation/{{assembler}}/{{filter_source}}/{{sample_id}}/transdecoder/{{td_db}}-taxaDB.{suff}", suff = ["index","dbtype"]),
        aln = expand("results/annotation/{{assembler}}/{{filter_source}}/{{sample_id}}/transdecoder/{{td_db}}-taxaDB_aln/first.{suff}", suff = ["index","dbtype"]),
    log:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/transdecoder/mmseqs_firstpass_taxonomy.{td_db}.log"
    threads: 10
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    conda: "../../envs/mmseqs.yaml"
    params:
        tmp = lambda wildcards: f"{os.environ.get("TMPDIR", "scratch")}/mmseqs_firstpass_taxonomy.{wildcards.assembler}.{wildcards.sample_id}.{wildcards.td_db}", 
        split_memory_limit = lambda wildcards, resources: int(resources.mem_mb*.8),
        output = lambda wildcards, output: f"{os.path.dirname(output[0])}/{wildcards.td_db}-taxaDB",
        ranks = "superkingdom,kingdom,phylum,class,order,family,genus,species",
        aln_dir = lambda wildcards, output: os.path.dirname(output.aln[0])
    resources:
        mem_mb = 3600,
        tasks = 1
    shell:
        """
        mkdir -p {params.tmp}
        mmseqs taxonomy {input.query} {input.target} {params.output} {params.tmp} -e 1e-5 -a 1 \
            --lca-mode 3 --tax-output-mode 0 --lca-ranks {params.ranks} --tax-lineage 1 \
            --split-memory-limit {params.split_memory_limit}M --threads {resources.tasks} > {log} 2>&1
        mv {params.tmp}/latest/first.* {params.aln_dir}
        rm -r {params.tmp}
        """

rule mmseqs_convertali:
    """
    Convert mmseqs alignment to m8 format for use with transdecoder
    """
    input:
        alignment = rules.mmseqs_firstpass_taxonomy.output.aln,
        target = os.path.join(config["mmseqs_db_dir"], "{td_db}"),
        query = rules.mmseqs_createquerydb_longorfs.output,
    output:
        m8 = "results/annotation/{assembler}/{filter_source}/{sample_id}/transdecoder/{td_db}.m8"
    log:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/transdecoder/{td_db}.convertali.log"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    conda: "../../envs/mmseqs.yaml"
    params:
        alignment = lambda wildcards, input: os.path.splitext(input.alignment[0])[0]
    shell:
        """
        mmseqs convertalis {input.query} {input.target} {params.alignment} {output.m8} --threads 1 --format-mode 0 > {log} 2>&1
        """

rule transdecoder_predict:
    """
    Run transdecoder predict using homology search results for ORF retention criteria
    """
    input:
        fa = "results/assembly/{assembler}/{filter_source}/{sample_id}/final.fa",
        m8 = expand("results/annotation/{{assembler}}/{{filter_source}}/{{sample_id}}/transdecoder/{td_db}.m8", td_db = config["transdecoder_homology_db"])
    output:
        expand("results/annotation/{{assembler}}/{{filter_source}}/{{sample_id}}/transdecoder/final.fa.transdecoder.{suff}", 
            suff = ["pep","gff3","cds","bed"])
    log:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/transdecoder/predict.log"
    container: "docker://trinityrnaseq/transdecoder:5.7.1"
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        TransDecoder.Predict -t {input.fa} --retain_blastp_hits {input.m8} -O {params.output_dir} > {log} 2>&1
        """

rule mmseqs_createtsv_first:
    """
    Create tsv file from mmseqs taxonomy output
    """
    input:
        query = rules.mmseqs_createquerydb_longorfs.output,
        result = rules.mmseqs_firstpass_taxonomy.output.tax
    output:
        tsv = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{td_db}/firstpass.tsv"
    log:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{td_db}/mmseqs_createtsv_first.log"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    conda: "../../envs/mmseqs.yaml"
    params:
        result = lambda wildcards, input: f"{os.path.dirname(input.result[0])}/{wildcards.td_db}-taxaDB",
    threads: 1
    resources:
        mem_mb = 1000,
        tasks = 1
    shell:
        """
        mmseqs createtsv {input.query} {params.result} {output.tsv} --threads {resources.tasks} > {log} 2>&1
        """

rule parse_mmseqs_first:
    """
    Parse mmseqs taxonomy output, adds columns with ranks
    """
    input:
        tsv = rules.mmseqs_createtsv_first.output,
    output:
        tsv = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{td_db}/firstpass.parsed.tsv"
    log:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{td_db}/firstpass.parsed.log"
    params:
        script = workflow.source_path("../utils/parse_mmseqs.py"),
        ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    shell:
        """
        python {params.script} -i {input.tsv} -o {output.tsv} -r {params.ranks} > {log} 2>&1
        """

rule firstpass_fungal_proteins:
    """
    Extract proteins from first taxonomy pass that are classified as fungi
    """
    input:
        parsed = rules.parse_mmseqs_first.output.tsv,
        fa = rules.transdecoder_predict.output[0]
    output:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{td_db}/firstpass.fungal.faa"
    params:
        idlist = lambda wildcards, input: f"{os.path.dirname(input.parsed)}/fungal.ids"
    run:
        import pandas as pd
        df = pd.read_csv(input.parsed, sep="\t", header=0, index_col=0)
        fungi = df.loc[df["kingdom"]=="Fungi"].index
        with open(params.idlist, "w") as f:
            for i in fungi:
                f.write(i + "\n")
        shell("seqkit grep -f {params.idlist} {input.fa} > {output}")
        os.remove(params.idlist)

rule mmseqs_secondpass_taxonomy:
    """
    Run mmseqs taxonomy on fungal proteins
    """
    input:
        query = expand("results/annotation/{{assembler}}/{{filter_source}}/{{sample_id}}/taxonomy/{td_db}/firstpass.fungal.faa",td_db = config["transdecoder_homology_db"]),
        target = get_mmseq_taxdb,
    output:
        expand("results/annotation/{{assembler}}/{{filter_source}}/{{sample_id}}/taxonomy/{mmseqs_db}/secondpass-taxresult_{suff}", suff = ["lca.tsv", "report","tophit_aln","tophit_report"], mmseqs_db = config["mmseqs_db"])
    log:
        expand("results/annotation/{{assembler}}/{{filter_source}}/{{sample_id}}/taxonomy/{mmseqs_db}/mmseqs_secondpass_taxonomy.log", mmseqs_db = config["mmseqs_db"])
    params:
        output = lambda wildcards, output: f"{os.path.dirname(output[0])}/secondpass-taxresult",
        tmp = lambda wildcards: f"{os.environ.get("TMPDIR", "scratch")}/mmseqs_secondpass_taxonomy_co.{wildcards.assembler}.{wildcards.sample_id}", 
        ranks = "superkingdom,kingdom,phylum,class,order,family,genus,species",
        split_memory_limit = lambda wildcards, resources: int(resources.mem_mb*.8),
        target = lambda wildcards, input: (input.target).replace("_taxonomy", "")
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    conda: "../../envs/mmseqs.yaml"
    threads: 10
    resources:
        mem_mb = 3600,
        tasks = 1
    shell:
        """
        mkdir -p {params.tmp}
        mmseqs easy-taxonomy {input.query} {params.target} {params.output} {params.tmp} \
            --lca-ranks {params.ranks} --lca-mode 3 --tax-lineage 1 --split-memory-limit {params.split_memory_limit}M \
            --threads {resources.tasks} > {log} 2>&1
        """

rule parse_mmseqs_second:
    """
    Parse mmseqs taxonomy output from second run
    """
    input:
        taxresult = rules.mmseqs_secondpass_taxonomy.output,
    output:
        tsv = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{mmseqs_db}/secondpass.parsed.tsv"
    log:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{mmseqs_db}/secondpass.parsed.log"
    params:
        script = workflow.source_path("../utils/parse_mmseqs.py"),
        tsv = lambda wildcards, input: os.path.join(os.path.dirname(input.taxresult[0]), "secondpass-taxresult_lca.tsv"),
        ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    shell:
        """
        python {params.script} -i {params.tsv} -o {output.tsv} -r {params.ranks} > {log} 2>&1
        """

rule secondpass_fungal_proteins:
    input:
        parsed = expand("results/annotation/{{assembler}}/{{filter_source}}/{{sample_id}}/taxonomy/{mmseqs_db}/secondpass.parsed.tsv", mmseqs_db = config["mmseqs_db"]),
        genecall = rules.transdecoder_predict.output,
    output:
        expand("results/annotation/{{assembler}}/{{filter_source}}/{{sample_id}}/genecall/fungal.{suff}", 
            suff = ["faa","gff3","cds","bed","taxonomy.tsv"])
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        indir = lambda wildcards, input: os.path.dirname(input.genecall[0]),
    shadow: "minimal"
    run:
        write_fungal_proteins(input.parsed, params.indir, params.outdir)

###################
## READ COUNTING ##
###################
rule featurecount:
    input:
        gff = "results/annotation/{assembler}/{filter_source}/{sample_id}/frame_selection/final.reformat.gff",
        bam = "results/map/{assembler}/{filter_source}/{sample_id}/{sample_id}.bam"
    output:
        cnt = "results/annotation/{assembler}/{filter_source}/{sample_id}/featureCounts/fc.tab",
        summary = "results/annotation/{assembler}/{filter_source}/{sample_id}/featureCounts/fc.tab.summary"
    resources:
        runtime = lambda wildcards, attempt: attempt*attempt*20
    conda: "../../envs/featurecount.yaml"
    params:
        setting = config["fc_params"]
    threads: 4
    shell:
        """
        featureCounts -T {threads} {params.setting} -a {input.gff} -o {output.cnt} {input.bam}
        """

rule normalize_featurecount:
    input:
        fc="results/annotation/{assembler}/{filter_source}/{sample_id}/featureCounts/fc.tab",
        stats="results/{filter_source}/{sample_id}/{sample_id}.stats.tsv"
    output:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/featureCounts/fc.tpm.tab",
        "results/annotation/{assembler}/{filter_source}/{sample_id}/featureCounts/fc.raw.tab"
    params:
        s = "{sample_id}",
        script = "source/utils/featureCountsTPM.py"
    run:
        df = pd.csv(input.stats[0], sep="\t")
        rl = df.avg_len.mean()
        shell("python {params.script} --rl {rl} -i {input[0]} -o {output[0]} --rc {output[1]} --sampleName {params.s}")
        
############################
## DIAMOND BLAST SEARCHES ##
############################
rule diamond_similarity_search:
    input:
        faa = "results/annotation/{assembler}/{filter_source}/{sample_id}/frame_selection/final.reformat.faa",
        db = "resources/diamond/{db}.dmnd"
    output:
        blast_out = "results/annotation/{assembler}/{filter_source}/{sample_id}/similarity_search/blastp_{db}.out"
    params:
        tmp_out = "$TMPDIR/{sample_id}/{assembler}/{filter_source}/blastp_{db}.out",
        tmp_dir = "$TMPDIR/{sample_id}/{assembler}/{filter_source}"
    resources:
        runtime=lambda wildcards, attempt: attempt**3*5*60
    threads: 20
    shell:
        """
        mkdir -p {params.tmp_dir}
        mkdir -p /dev/shm/$SLURM_JOB_ID
        diamond blastp -d {input.db} -c 1 -b 20 --query-cover 50.000000 --subject-cover 50.000000 --evalue 0.000010 \
        --more-sensitive --top 3 -q {input.faa} -o {params.tmp_out} -p {threads} --tmpdir /dev/shm/$SLURM_JOB_ID \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle
        mv {params.tmp_out} {output.blast_out}
        conda deactivate
        """

###################
## EGGNOG-MAPPER ##
###################
rule emapper_search:
    input:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/frame_selection/final.reformat.faa",
        "resources/eggnog/eggnog.db"
    output:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/annotation_results.emapper.seed_orthologs",
    params:
        out_base = "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/annotation_results",
        resource_dir = "resources/eggnog",
        out = "annotation_results",
        tmpdir = "$TMPDIR/eggnog/{assembler}/{filter_source}/{sample_id}",
        tmp_out = "$TMPDIR/eggnog/{assembler}/{filter_source}/{sample_id}/annotation_results",
        flags = "-m diamond --no_annot --no_file_comments",
        outdir = "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG"
    log: "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/emapper.log"
    threads: 20
    shadow: "shallow"
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60
    conda: "../../envs/emapper.yaml"
    shell:
        """
        mkdir -p {params.tmpdir}
        mkdir -p {params.outdir}
        emapper.py {params.flags} --cpu {threads} -i {input[0]} -o {params.out} --temp_dir {params.tmpdir} \
        --output_dir {params.tmpdir} --data_dir {params.resource_dir} 2>{log}
        mv {params.tmp_out}.emapper.seed_orthologs {output[0]}
        """

rule emapper_annotate_hits:
    input:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/annotation_results.emapper.seed_orthologs",
        "resources/eggnog/eggnog.db"
    output:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/annotation_results.emapper.annotations"
    params:
        resource_dir = "resources/eggnog",
        tmpdir = "$TMPDIR/eggnog/{assembler}/{filter_source}/{sample_id}",
        out = "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/annotation_results",
        flags = "--no_file_comments"
    log: "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/emapper.annotate.log"
    threads: 10
    shadow: "minimal"
    conda: "../../envs/emapper.yaml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 60
    shell:
        """
        #Copy eggnog.db to /dev/shm
        mkdir -p /dev/shm/$SLURM_JOB_ID
        cp {input[1]} /dev/shm/$SLURM_JOB_ID
        #Run annotation of hits table
        emapper.py {params.flags} --cpu {threads} --annotate_hits_table {input[0]} -o {params.out} --data_dir /dev/shm/$SLURM_JOB_ID --usemem 2>{log}
        # Clean up
        rm -rf /dev/shm/$SLURM_JOB_ID
        """

##################
## PARSE EGGNOG ##
##################
rule parse_eggnog:
    input:
        f = "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/annotation_results.emapper.annotations",
        db = expand("resources/kegg/{f}",
            f = ["kegg_ec2pathways.tsv","kegg_ko2ec.tsv",
                 "kegg_ko2pathways.tsv","kegg_kos.tsv","kegg_modules.tsv","kegg_pathways.tsv"])
    output:
        expand("results/annotation/{{assembler}}/{{filter_source}}/{{sample_id}}/eggNOG/{db}.parsed.tab",
            db = ["enzymes","kos","modules","pathways","tc","cazy"])
    log: "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/parser.log"
    params:
        src = "source/utils/eggnog-parser.py",
        dldir = "resources/kegg",
        outdir = "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/"
    shell:
        """
        python {params.src} parse {params.dldir} {input.f} {params.outdir} 2>{log}
        """

rule quantify_eggnog:
    input:
        abundance = "results/annotation/{assembler}/{filter_source}/{sample_id}/featureCounts/fc.{fc}.tab",
        parsed = "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/{db}.parsed.tab"
    output:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/{db}.{fc}.tsv"
    params:
        src = "source/utils/eggnog-parser.py"
    shell:
        """
        python {params.src} quantify {input.abundance} {input.parsed} {output[0]}
        """

rule quantify_eggnog_normalized:
    """Normalize modules/pathways further by the number of KEGG orthologs belonging to each module/pathway"""
    input:
        abundance = "results/annotation/{assembler}/{filter_source}/{sample_id}/featureCounts/fc.{fc}.tab",
        parsed = "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/{db}.parsed.tab",
        norm = "resources/kegg/kegg_ko2{db}.tsv"
    output:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/eggNOG/{db}.norm.{fc}.tsv"
    params:
        src = "source/utils/eggnog-parser.py"
    shell:
        """
        python {params.src} quantify --normalize {input.norm} {input.abundance} {input.parsed} {output[0]}
        """

rule eggnog_merge_and_sum:
    input:
        expand("results/annotation/{{assembler}}/{{filter_source}}/{sample_id}/eggNOG/{{db}}.{{fc}}.tsv",
            sample_id = samples.keys())
    output:
        "results/collated/{assembler}/{filter_source}/eggNOG/{db}.{fc}.tsv"
    params:
        src = "source/utils/eggnog-parser.py"
    shell:
        """
        python {params.src} merge --sum {input} {output[0]}
        """

###########
## DBCAN ##
###########
rule dbcan_scan:
    input:
        fasta = "results/annotation/{assembler}/{filter_source}/{sample_id}/frame_selection/final.reformat.faa",
        db = "resources/dbCAN/dbCAN-fam-HMMs.txt",
        pressed = expand("resources/dbCAN/dbCAN-fam-HMMs.txt.{suffix}", suffix = ["h3f","h3i","h3m","h3p"])
    output:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/dbCAN/dbCAN.out.dm",
        "results/annotation/{assembler}/{filter_source}/{sample_id}/dbCAN/dbCAN.out"
    threads: 1
    conda: "../../envs/hmmer.yaml"
    resources:
        runtime=lambda wildcards, attempt: attempt**2*30
    shell:
        """
        hmmscan --cpu {threads} --domtblout {output[0]} {input.db} {input.fasta} > {output[1]}
        """

rule dbcan_parse:
    input:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/dbCAN/dbCAN.out.dm",
        "source/external/hmmscan-parser.sh"
    output:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/dbCAN/dbCAN.parsed.tsv"
    shell:
        """
        sh {input[1]} {input[0]} > {output[0]}
        """

rule sum_dbcan:
    input:
        tpm = "results/annotation/{assembler}/{filter_source}/{sample_id}/featureCounts/fc.tpm.tab",
        raw = "results/annotation/{assembler}/{filter_source}/{sample_id}/featureCounts/fc.raw.tab",
        dbcan = "results/annotation/{assembler}/{filter_source}/{sample_id}/dbCAN/dbCAN.parsed.tsv"
    output:
        tpm = "results/annotation/{assembler}/{filter_source}/{sample_id}/dbCAN/dbCAN.tpm.tsv",
        raw = "results/annotation/{assembler}/{filter_source}/{sample_id}/dbCAN/dbCAN.raw.tsv"
    params:
        evalue = config["dbCAN_eval"],
        coverage = config["dbCAN_cov"]
    run:
        tpm = pd.read_table(input.tpm, index_col=0)
        raw = pd.read_table(input.raw, index_col=0)
        dbcan = pd.read_table(input.dbcan, header=None, names=["HMM","HMM length","Query","Query length","evalue","hmmstart","hmmend","querystart","queryend","coverage"],
            index_col=2, dtype={"evalue": float, "cov": float})
        l = len(dbcan.index)
        dbcan = dbcan.loc[(dbcan.evalue<float(params.evalue))&(dbcan.coverage>float(params.coverage))]
        print("Filtered {}/{} records from dbCAN".format(l-len(dbcan.index),l))
        dbcan_tpm = pd.merge(tpm, dbcan, right_index=True, left_index=True, how="left")
        dbcan_raw = pd.merge(raw, dbcan, right_index=True, left_index=True, how="left")
        dbcan_tpm_sum = dbcan_tpm.loc[:,["HMM",wildcards.sample_id]].groupby("HMM").sum()
        dbcan_raw_sum = dbcan_raw.loc[:,["HMM",wildcards.sample_id]].groupby("HMM").sum()
        dbcan_tpm_sum.to_csv(output.tpm, sep="\t")
        dbcan_raw_sum.to_csv(output.raw, sep="\t")

rule collate_dbcan:
    input:
        expand("results/annotation/{{assembler}}/{{filter_source}}/{sample_id}/dbCAN/dbCAN.{{fc}}.tsv",
            sample_id = samples.keys())
    output:
        "results/collated/{assembler}/{filter_source}/dbCAN/dbCAN.{fc}.tsv"
    run:
        df = pd.DataFrame()
        for f in input:
            sample_id = f.split("/")[-3]
            _df = pd.read_table(f, index_col=0)
            df = pd.merge(df, _df, left_index=True, right_index=True, how="outer")
        df.to_csv(output[0], sep="\t")

##########################
## TAXONOMIC ANNOTATION ##
##########################

rule sum_taxonomy:
    input:
        gene_tax = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/gene_taxonomy.tsv",
        tpm = "results/annotation/{assembler}/{filter_source}/{sample_id}/featureCounts/fc.tpm.tab",
        raw = "results/annotation/{assembler}/{filter_source}/{sample_id}/featureCounts/fc.raw.tab"
    output:
        tpm = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/taxonomy.tpm.tsv",
        raw = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/taxonomy.raw.tsv"
    run:
        ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
        gene_tax = pd.read_table(input.gene_tax, header=0, sep="\t", index_col=0)
        tpm = pd.read_table(input.tpm, header=0, sep="\t", index_col=0)
        raw = pd.read_table(input.raw, header=0, sep="\t", index_col=0)
        gene_tax_tpm = pd.merge(gene_tax, tpm, left_index = True, right_index = True)
        species_tpm = gene_tax_tpm.groupby(ranks).sum().reset_index()
        gene_tax_raw = pd.merge(gene_tax, raw, left_index = True, right_index = True)
        species_raw = gene_tax_raw.groupby(ranks).sum().reset_index()
        # Write results
        species_tpm.to_csv(output.tpm, sep="\t", index=False)
        species_raw.to_csv(output.raw, sep="\t", index=False)


def merge_tax(files):
    df = pd.DataFrame()
    for f in files:
        _df = pd.read_table(f, header=0, sep="\t")
        ranks = list(_df.columns[0:-1])
        # Set unique index
        index = []
        for j in _df.index:
             index.append(";".join(_df.loc[j,ranks].values))
        _df.index = index
        _df.drop(ranks, axis=1, inplace=True)
        df = pd.merge(df, _df, left_index=True, right_index=True, how="outer")
    cols = {}
    # Reset columns
    for item in df.index:
        names = item.split(";")
        for i, name in enumerate(names):
            rank = ranks[i]
            try:
                cols[rank].append(name)
            except KeyError:
                cols[rank] = [name]
    df.index = list(range(0,len(df)))
    df = pd.merge(pd.DataFrame(cols)[ranks], df, left_index=True, right_index=True)
    df.fillna(0, inplace=True)
    return df

rule collate_taxonomy:
    input:
        tpm = expand("results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/taxonomy.tpm.tsv",
            assembler = config["assembler"], filter_source = config["read_source"], sample_id = samples.keys()),
        raw = expand("results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/taxonomy.raw.tsv",
            assembler = config["assembler"], filter_source = config["read_source"], sample_id = samples.keys())
    output:
        tpm = "results/collated/{assembler}/{filter_source}/taxonomy/taxonomy.tpm.tsv",
        raw = "results/collated/{assembler}/{filter_source}/taxonomy/taxonomy.raw.tsv"
    run:
        tpm = merge_tax(input.tpm)
        raw = merge_tax(input.raw)
        tpm.to_csv(output.tpm, sep="\t", index=False)
        raw.to_csv(output.raw, sep="\t", index=False)
