localrules:
    mmseqs_convertali_co,
    parse_mmseqs_first_co,
    parse_mmseqs_second_co,
    firstpass_fungal_proteins_co,
    secondpass_fungal_proteins_co,
    collate_featurecount_co,
    dbcan_parse_co,
    normalize_featurecount_co,
    parse_eggnog_co,
    sum_dbcan_co,
    sum_taxonomy_co,
    taxonomy_featurecount_co,
    dbcan_tax_annotations,
    eggnog_tax_annotations,

##########################################
##              CO ASSEMBLIES           ##
##########################################
# Below are rules for annotation of co-assemblies

###################################
## FRAME SELECTION CO-ASSEMBLIES ##
###################################
rule transdecoder_longorfs_co:
    """
    Generate longest ORFs from co-assembly with transdecoder
    """
    input:
        fa = "results/co-assembly/{assembler}/{assembly}/final.fa",
    output:
        expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/transdecoder/final.fa.transdecoder_dir/longest_orfs.{suff}",
            suff = ["pep","gff3","cds"])
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/longest_orfs.log"
    container: "docker://trinityrnaseq/transdecoder:5.7.1"
    params:
        output_dir = lambda wildcards: f"results/annotation/co-assembly/{wildcards.assembler}/{wildcards.assembly}/transdecoder"
    resources:
        tasks= 1,
        runtime = 240
    shell:
        """
        TransDecoder.LongOrfs -t {input.fa} -O {params.output_dir} > {log} 2>&1
        """

rule mmseqs_createquerydb_longorfs_co:
    """
    Create mmseqs database from longest ORFs
    """
    input:
        query = rules.transdecoder_longorfs_co.output[0],
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/longorfs-queryDB"
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/longorfs-queryDB.log"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    conda: "../../envs/mmseqs.yaml"
    shell:
        """
        mmseqs createdb {input} {output} > {log} 2>&1
        """

rule mmseqs_firstpass_taxonomy_co:
    """
    Run mmseqs taxonomy on longest ORFs, also store the alignment results
    for use with transdecoder
    """
    input:
        query = rules.mmseqs_createquerydb_longorfs_co.output,
        target = os.path.join(config["mmseqs_db_dir"], "{td_db}")
    output:
        tax = expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/transdecoder/{{td_db}}-taxaDB.{suff}", suff = ["index","dbtype"]),
        aln = expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/transdecoder/{{td_db}}-taxaDB_aln/first.{suff}", suff = ["index","dbtype"]),
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/mmseqs_firstpass_taxonomy_co.{td_db}.log"
    threads: 10
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    conda: "../../envs/mmseqs.yaml"
    params:
        tmp = lambda wildcards: f"{os.environ.get("TMPDIR", "scratch")}/mmseqs_firstpass_taxonomy_co.{wildcards.assembler}.{wildcards.assembly}.{wildcards.td_db}", 
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
            --split-memory-limit {params.split_memory_limit}M --threads {threads} > {log} 2>&1
        mv {params.tmp}/latest/first.* {params.aln_dir}
        rm -r {params.tmp}
        """

rule mmseqs_convertali_co:
    """
    Convert mmseqs alignment to m8 format for transdecoder
    """
    input:
        alignment = rules.mmseqs_firstpass_taxonomy_co.output.aln,
        target = os.path.join(config["mmseqs_db_dir"], "{td_db}"),
        query = rules.mmseqs_createquerydb_longorfs_co.output
    output:
        m8 = "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/{td_db}.m8"
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/{td_db}.convertali.log"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    conda: "../../envs/mmseqs.yaml"
    params:
        alignment = lambda wildcards, input: os.path.splitext(input.alignment[0])[0]
    shell:
        """
        mmseqs convertalis {input.query} {input.target} {params.alignment} {output.m8} --threads 1 --format-mode 0 > {log} 2>&1
        """

rule transdecoder_predict_co:
    """
    Run transdecoder predict using homology search results for ORF retention criteria
    """
    input:
        fa = "results/co-assembly/{assembler}/{assembly}/final.fa",
        m8 = expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/transdecoder/{td_db}.m8",
            td_db = config["transdecoder_homology_db"])
    output:
        expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/transdecoder/final.fa.transdecoder.{suff}",
            suff = ["pep","gff3","cds","bed"])
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/predict.log"
    container: "docker://trinityrnaseq/transdecoder:5.7.1"
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        TransDecoder.Predict -t {input.fa} --retain_blastp_hits {input.m8} -O {params.output_dir} > {log} 2>&1
        """

rule mmseqs_createtsv_first_co:
    """
    Create tsv file from mmseqs taxonomy output
    """
    input:
        query = rules.mmseqs_createquerydb_longorfs_co.output,
        result = rules.mmseqs_firstpass_taxonomy_co.output.tax
    output:
        tsv = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{td_db}/firstpass.tsv"
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{td_db}/mmseqs_createtsv_first_co.log"
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
        mmseqs createtsv {input.query} {params.result} {output.tsv} --threads {threads} > {log} 2>&1
        """

rule parse_mmseqs_first_co:
    """
    Parse mmseqs taxonomy output, adds columns with ranks
    """
    input:
        tsv = rules.mmseqs_createtsv_first_co.output,
    output:
        tsv = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{td_db}/firstpass.parsed.tsv"
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{td_db}/firstpass.parsed.log"
    params:
        script = workflow.source_path("../utils/parse_mmseqs.py"),
        ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    shell:
        """
        python {params.script} -i {input.tsv} -o {output.tsv} -r {params.ranks} > {log} 2>&1
        """

rule firstpass_fungal_proteins_co:
    """
    Extract proteins from first taxonomy pass that are classified as fungi
    """
    input:
        parsed = rules.parse_mmseqs_first_co.output.tsv,
        fa = rules.transdecoder_predict_co.output[0]
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{td_db}/firstpass.fungal.faa"
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

rule mmseqs_secondpass_taxonomy_co:
    """
    Run mmseqs taxonomy on fungal proteins
    """
    input:
        query = expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/taxonomy/{td_db}/firstpass.fungal.faa", td_db = config["transdecoder_homology_db"]),
        target = get_mmseq_taxdb
    output:
        expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/taxonomy/{mmseqs_db}/secondpass-taxresult_{suff}", suff = ["lca.tsv", "report","tophit_aln","tophit_report"], mmseqs_db = config["mmseqs_db"])
    log:
        expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/taxonomy/{mmseqs_db}/mmseqs_secondpass_taxonomy_co.log", mmseqs_db = config["mmseqs_db"])
    params:
        output = lambda wildcards, output: f"{os.path.dirname(output[0])}/secondpass-taxresult",
        tmp = lambda wildcards: f"{os.environ.get("TMPDIR", "scratch")}/mmseqs_secondpass_taxonomy_co.{wildcards.assembler}.{wildcards.assembly}", 
        ranks = "superkingdom,kingdom,phylum,class,order,family,genus,species",
        split_memory_limit = lambda wildcards, resources: int(resources.mem_mb*.8),
        target = lambda wildcards, input: (input.target).replace("_taxonomy", "")
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    conda: "../../envs/mmseqs.yaml"
    threads: 10
    shell:
        """
        mkdir -p {params.tmp}
        mmseqs easy-taxonomy {input.query} {params.target} {params.output} {params.tmp} \
            --lca-ranks {params.ranks} --lca-mode 3 --tax-lineage 1 --split-memory-limit {params.split_memory_limit}M \
            --threads {threads} > {log} 2>&1
        """

rule parse_mmseqs_second_co:
    """
    Parse mmseqs taxonomy output from second run
    """
    input:
        taxresult = rules.mmseqs_secondpass_taxonomy_co.output,
    output:
        tsv = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{mmseqs_db}/secondpass.parsed.tsv"
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{mmseqs_db}/secondpass.parsed.log"
    params:
        script = workflow.source_path("../utils/parse_mmseqs.py"),
        tsv = lambda wildcards, input: os.path.join(os.path.dirname(input.taxresult[0]), "secondpass-taxresult_lca.tsv"),
        ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    shell:
        """
        python {params.script} -i {params.tsv} -o {output.tsv} -r {params.ranks} > {log} 2>&1
        """

rule secondpass_fungal_proteins_co:
    input:
        parsed = expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/taxonomy/{mmseqs_db}/secondpass.parsed.tsv", mmseqs_db = config["mmseqs_db"]),
        genecall = rules.transdecoder_predict_co.output,
    output:
        expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/genecall/fungal.{suff}",
            suff = ["faa","gff3","cds","bed","taxonomy.tsv"])
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        indir = lambda wildcards, input: os.path.dirname(input.genecall[0]),
    shadow: "minimal"
    run:
        write_fungal_proteins(input.parsed[0], params.indir, params.outdir)

#################################
## READ COUNTING CO-ASSEMBLIES ##
#################################
rule featurecount_co:
    input:
        gff = "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/final.fa.transdecoder.gff3",
        bam = "results/map/co-assembly/{assembler}/{assembly}/{sample_id}.bam"
    output:
        cnt = "results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.fc.tab",
        summary = "results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.fc.tab.summary"
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.fc.log"
    resources:
        runtime = 30,
        tasks = 1,
        mem_mb = 1000
    conda: "../../envs/featurecount.yaml"
    params:
        setting = config["fc_params"]
    container: "docker://quay.io/biocontainers/subread:2.0.8--h577a1d6_0"
    threads: 4
    shell:
        """
        featureCounts -T {threads} {params.setting} -a {input.gff} -o {output.cnt} {input.bam} > {log} 2>&1
        """

rule normalize_featurecount_co:
    input:
        fc="results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.fc.tab",
        stats=expand("results/{source}/{{sample_id}}/{{sample_id}}.stats.tsv", source="filtered" if config["filter_reads"] else "unfiltered")
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.tpm.tab",
        "results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.raw.tab"
    params:
        s = "{sample_id}",
        script = "source/utils/featureCountsTPM.py"
    run:
        df = pd.read_csv(input.stats[0], sep="\t")
        rl = df.avg_len.mean()
        shell("python {params.script} --rl {rl} -i {input.fc} -o {output[0]} --rc {output[1]} --sampleName {params.s}")

rule collate_featurecount_co:
    input:
        expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/featureCounts/{sample_id}.{{fc}}.tab",
            sample_id = samples.keys())
    output:
        "results/collated/co-assembly/{assembler}/{assembly}/abundance/{assembly}.{fc}.tsv"
    run:
        df = pd.DataFrame()
        for f in input:
            _df = pd.read_table(f, sep="\t", header=0, index_col=0)
            df = pd.merge(df, _df, right_index=True, left_index=True, how="outer")
        df.to_csv(output[0], sep="\t", index=True, header=True)

#################################
## EGGNOG-MAPPER CO-ASSEMBLIES ##
#################################
rule emapper_search_co:
    input:
        "results/annotation/co-assembly/{assembler}/{assembly}/genecall/fungal.faa",
        "resources/eggnog/eggnog.db"
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/annotation_results.emapper.seed_orthologs"
    params:
        out_base = "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/annotation_results",
        resource_dir = "resources/eggnog",
        out = "annotation_results",
        tmpdir = "$TMPDIR/eggnog/{assembly}",
        tmp_out = "$TMPDIR/eggnog/{assembly}/annotation_results",
        flags = "-m diamond --no_annot --no_file_comments",
        outdir = "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/"
    log: "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/emapper.log"
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

rule emapper_annotate_hits_co:
    input:
        "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/annotation_results.emapper.seed_orthologs",
        "resources/eggnog/eggnog.db",
        "resources/eggnog/eggnog_proteins.dmnd"
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/annotation_results.emapper.annotations"
    params:
        resource_dir = "resources/eggnog",
        tmpdir = "$TMPDIR/eggnog/{assembly}/",
        out = "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/annotation_results",
        flags = "--no_file_comments"
    log: "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/emapper.annotate.log"
    threads: 10
    shadow: "minimal"
    conda: "../../envs/emapper.yaml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 60
    shell:
        """
        #Copy eggnog.db to /dev/shm
        mkdir -p /dev/shm/$SLURM_JOB_ID
        cp {input[1]} {input[2]} /dev/shm/$SLURM_JOB_ID
        #Run annotation of hits table
        emapper.py {params.flags} --cpu {threads} --annotate_hits_table {input[0]} -o {params.out} \
        --data_dir /dev/shm/$SLURM_JOB_ID --usemem 2>{log}
        # Clean up
        rm -rf /dev/shm/$SLURM_JOB_ID
        """

################################
## PARSE EGGNOG CO-ASSEMBLIES ##
################################
rule parse_eggnog_co:
    input:
        f = "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/annotation_results.emapper.annotations",
        db = expand("resources/kegg/{f}",
            f = ["kegg_ec2pathways.tsv","kegg_ko2ec.tsv",
                 "kegg_ko2pathways.tsv","kegg_kos.tsv","kegg_modules.tsv","kegg_pathways.tsv"])
    log: "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/parser.log"
    output:
        expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/eggNOG/{db}.parsed.tab",
            db = ["enzymes","ko","modules","pathways","tc","cazy"])
    params:
        src = "source/utils/eggnog-parser.py",
        dldir = "resources/kegg",
        outdir = "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/"
    shell:
        """
        python {params.src} parse {params.dldir} {input.f} {params.outdir} 2>{log}
        """

rule quantify_eggnog_co:
    input:
        abundance = expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/featureCounts/{sample_id}.{{fc}}.tab",
            sample_id = samples.keys()),
        parsed = "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/{db}.parsed.tab",
    output:
        "results/collated/co-assembly/{assembler}/{assembly}/eggNOG/{db}.{fc}.tsv"
    params:
        src = "source/utils/eggnog-parser.py",
        tmpdir = os.path.join(os.path.expandvars("$TMPDIR"), "{assembly}", "{db}")
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        mkdir -p {params.tmpdir}
        for f in {input.abundance};
        do
            base=$(basename $f)
            sample=$(echo -e $base | sed "s/.{wildcards.fc}.tab//g")
            python {params.src} quantify $f {input.parsed} {params.tmpdir}/$sample.tsv
        done
        python {params.src} merge --sum {params.tmpdir}/*.tsv {output[0]}
        rm -rf {params.tmpdir}
        """

rule quantify_eggnog_normalized_co:
    """
    Normalize modules/pathways further by the number of KEGG orthologs belonging to each module/pathway
    """
    input:
        abundance = expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/featureCounts/{sample_id}.{{fc}}.tab",
            sample_id = samples.keys()),
        parsed = "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/{db}.parsed.tab",
        norm = "resources/kegg/kegg_ko2{db}.tsv"
    output:
        "results/collated/co-assembly/{assembler}/{assembly}/eggNOG/{db}.norm.{fc}.tsv"
    params:
        src = "source/utils/eggnog-parser.py",
        tmpdir = os.path.join(os.path.expandvars("$TMPDIR"), "{assembly}", "norm", "{db}")
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        mkdir -p {params.tmpdir}
        for f in {input.abundance};
        do
            base=$(basename $f)
            sample=$(echo -e $base | sed "s/.{wildcards.fc}.tab//g")
            python {params.src} quantify --normalize {input.norm} $f {input.parsed} {params.tmpdir}/$sample.tsv
        done
        python {params.src} merge --sum {params.tmpdir}/*.tsv {output[0]}
        rm -rf {params.tmpdir}
        """

rule eggnog_tax_annotations:
    """Collate gene annotations and abundances for a certain rank:taxon combination"""
    input:
        gene_tax = "results/annotation/co-assembly/{assembler}/{assembly}/genecall/fungal.taxonomy.tsv",
        parsed = "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/{db}.parsed.tab",
        abundance = "results/collated/co-assembly/{assembler}/{assembly}/abundance/{assembly}.{fc}.tsv"
    output:
        expand("results/collated/co-assembly/{{assembler}}/{{assembly}}/eggNOG_taxonomy/{tax_rank}.{tax_name}.{{db}}.{{fc}}.tsv",
               tax_rank = config["tax_rank"], tax_name = config["tax_name"])
    run:
        import pandas as pd
        abundance = pd.read_csv(input.abundance, header=0, sep="\t", index_col=0)
        parsed = pd.read_csv(input.parsed, header=0, sep="\t", index_col=0)
        gene_tax = pd.read_csv(input.gene_tax, header=0, sep="\t", index_col=0)
        features = len(set(parsed.loc[:,parsed.columns[0]]))
        orig_orfs = abundance.shape[0]
        try:
            orfs = list(set(gene_tax.loc[gene_tax[config["tax_rank"]]==config["tax_name"]].index))
        except KeyError:
            orfs = []
        if len(orfs) == 0:
            shell("touch {output}")
        else:
            abundance = abundance.loc[set(orfs).intersection(set(abundance.index))]
            parsed = parsed.loc[set(orfs).intersection(set(parsed.index))]
            tax_features = len(set(parsed.loc[:,parsed.columns[0]]))
            print("{}/{} features, {}/{} orfs matched".format(tax_features, features, len(orfs), orig_orfs))
            df = pd.merge(parsed, abundance, left_index=True, right_index=True)
            df_sum = df.groupby(list(parsed.columns)).sum().reset_index()
            df_sum.to_csv(output[0], sep="\t", index=False, header=True)

##########################
## TAXONOMIC ANNOTATION ##
##########################

rule sum_taxonomy_co:
    input:
        gene_tax = "results/annotation/co-assembly/{assembler}/{assembly}/genecall/fungal.taxonomy.tsv",
        abundance = "results/collated/co-assembly/{assembler}/{assembly}/abundance/{assembly}.{fc}.tsv"
    output:
        "results/collated/co-assembly/{assembler}/{assembly}/taxonomy/taxonomy.{fc}.tsv"
    run:
        ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
        gene_tax = pd.read_table(input.gene_tax, header=0, sep="\t", index_col=0)
        abundance = pd.read_table(input.abundance, header=0, sep="\t", index_col=0)
        gene_tax_abundance = pd.merge(gene_tax, abundance, left_index = True, right_index = True)
        species_abundance = gene_tax_abundance.groupby(ranks).sum().reset_index()
        # Write results
        species_abundance.to_csv(output[0], sep="\t", index=False)

#########################
## DBCAN CO-ASSEMBLIES ##
#########################
rule dbcan_scan_co:
    input:
        fasta = "results/annotation/co-assembly/{assembler}/{assembly}/genecall/fungal.faa",
        db = "resources/dbCAN/dbCAN-fam-HMMs.txt",
        pressed = expand("resources/dbCAN/dbCAN-fam-HMMs.txt.{suffix}", suffix = ["h3f","h3i","h3m","h3p"])
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/dbCAN/dbCAN.out.dm",
        "results/annotation/co-assembly/{assembler}/{assembly}/dbCAN/dbCAN.out"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    conda: "../../envs/hmmer.yaml"
    shell:
        """
        hmmscan --cpu {threads} --domtblout {output[0]} {input.db} {input.fasta} > {output[1]}
        """

rule dbcan_parse_co:
    input:
        "results/annotation/co-assembly/{assembler}/{assembly}/dbCAN/dbCAN.out.dm",
        "source/external/hmmscan-parser.sh"
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/dbCAN/dbCAN.parsed.tsv"
    shell:
        """
        sh {input[1]} {input[0]} > {output[0]}
        """

rule sum_dbcan_co:
    input:
        abundance = "results/collated/co-assembly/{assembler}/{assembly}/abundance/{assembly}.{fc}.tsv",
        dbcan = "results/annotation/co-assembly/{assembler}/{assembly}/dbCAN/dbCAN.parsed.tsv"
    output:
        "results/collated/co-assembly/{assembler}/{assembly}/dbCAN/dbCAN.{fc}.tsv"
    params:
        evalue = config["dbCAN_eval"],
        coverage = config["dbCAN_cov"]
    run:
        abundance = pd.read_table(input.abundance, index_col=0)
        samples = list(abundance.columns)
        dbcan = pd.read_table(input.dbcan, header=None, names=["HMM","HMM length","Query","Query length","evalue",
        "hmmstart","hmmend","querystart","queryend","coverage"],
            index_col=2, dtype={"evalue": float, "cov": float})
        l = len(dbcan.index)
        dbcan = dbcan.loc[(dbcan.evalue<float(params.evalue))&(dbcan.coverage>float(params.coverage))]
        print("Filtered {}/{} records from dbCAN".format(l-len(dbcan.index),l))
        dbcan_abundance = pd.merge(abundance, dbcan, right_index=True, left_index=True, how="left")
        dbcan_sum = dbcan_abundance.loc[:,["HMM"]+samples].groupby("HMM").sum()
        dbcan_sum.to_csv(output[0], sep="\t", index=True, header=True)

rule dbcan_tax_annotations:
    input:
        abundance = "results/collated/co-assembly/{assembler}/{assembly}/abundance/{assembly}.{fc}.tsv",
        parsed = "results/annotation/co-assembly/{assembler}/{assembly}/dbCAN/dbCAN.parsed.tsv",
        gene_tax = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/gene_taxonomy.tsv"
    output:
        expand("results/collated/co-assembly/{{assembler}}/{{assembly}}/dbCAN_taxonomy/{rank}.{taxon}.dbCAN.{{fc}}.tsv",
               rank = config["tax_rank"], taxon = config["tax_name"])
    run:
        import pandas as pd
        abundance = pd.read_csv(input.abundance, header=0, sep="\t", index_col=0)
        parsed = pd.read_csv(input.parsed, header=0, sep="\t", index_col=0)
        gene_tax = pd.read_csv(input.gene_tax, header=0, sep="\t", index_col=0)
        try:
            orfs = list(set(gene_tax.loc[gene_tax[config["tax_rank"]]==config["tax_name"]].index))
        except KeyError:
            orfs = []
        if len(orfs) == 0:
            shell("touch {output}")
        else:
            abundance = abundance.loc[set(orfs).intersection(set(abundance.index))]
            parsed = parsed.loc[set(orfs).intersection(set(parsed.index))]
            df = pd.merge(parsed, abundance, left_index=True, right_index=True)
            df_sum = df.groupby(list(parsed.columns)).sum().reset_index()
            df_sum.to_csv(output[0], sep="\t", index=False, header=True)
