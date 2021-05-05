localrules:
    collate_featurecount_co,
    create_blob_co,
    dbcan_parse_co,
    normalize_featurecount_co,
    parse_eggnog_co,
    plot_blob_co,
    reformat_fasta_co,
    sum_dbcan_co,
    sum_taxonomy_co,
    taxonomy_featurecount_co,
    dbcan_tax_annotations,
    eggnog_tax_annotations,
    transfer_taxonomy_co

##########################################
##              CO ASSEMBLIES           ##
##########################################
# Below are rules for annotation of co-assemblies

###################################
## FRAME SELECTION CO-ASSEMBLIES ##
###################################
rule frame_selection_co:
    input:
        fa = "results/co-assembly/{assembler}/{assembly}/final.fa",
    output:
        expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/frame_selection/final.{suffix}",
            suffix = ["gff","faa","fnn"])
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/frame_selection/gms.log"
    params:
        tmpdir = "$TMPDIR/{assembly}/",
        outdir = "results/annotation/co-assembly/{assembler}/{assembly}/frame_selection"
    threads: 2
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*2
    shell:
        """
        wd=$(pwd)
        mkdir -p {params.tmpdir}
        mkdir -p {params.outdir}
        cd {params.tmpdir}
        ln -s $wd/{input.fa} final.fa
        gmst.pl --output final --format GFF --fnn --faa final.fa
        rm final.fa
        mv final $wd/{params.outdir}/final.gff
        mv {params.tmpdir}/* $wd/{params.outdir}
        cd $wd
        rm -r {params.tmpdir}
        """

rule reformat_fasta_co:
    input:
        faa = "results/annotation/co-assembly/{assembler}/{assembly}/frame_selection/final.faa"
    output:
        faa = "results/annotation/co-assembly/{assembler}/{assembly}/frame_selection/final.reformat.faa",
        gff = "results/annotation/co-assembly/{assembler}/{assembly}/frame_selection/final.reformat.gff"
    params:
        assembly = "{assembly}"
    shell:
        """
        python source/utils/reformat_fasta.py {input.faa} {output.faa} {output.gff} {params.assembly}
        """

#################################
## READ COUNTING CO-ASSEMBLIES ##
#################################
rule featurecount_co:
    input:
        gff = "results/annotation/co-assembly/{assembler}/{assembly}/frame_selection/final.reformat.gff",
        bam = "results/map/co-assembly/{assembler}/{assembly}/{sample_id}.bam"
    output:
        cnt = "results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.fc.tab",
        summary = "results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.fc.tab.summary"
    resources:
        runtime = lambda wildcards, attempt: attempt*attempt*20
    params:
        setting = config["fc_params"]
    threads: 4
    shell:
        """
        featureCounts -T {threads} {params.setting} -a {input.gff} -o {output.cnt} {input.bam}
        """

rule normalize_featurecount_co:
    input:
        "results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.fc.tab",
        expand("results/sample_info/{{sample_id}}.{source}_read_lengths.tab",
            source = config["read_source"])
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.tpm.tab",
        "results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.raw.tab"
    params:
        s = "{sample_id}",
        script = "source/utils/featureCountsTPM.py"
    run:
        df = pd.read_table(input[1], index_col=0)
        rl = df.loc[wildcards.sample_id,"avg_len"]
        shell("python {params.script} --rl {rl} -i {input[0]} -o {output[0]} --rc {output[1]} --sampleName {params.s}")

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

rule taxonomy_featurecount_co:
    input:
        abundance = "results/collated/co-assembly/{assembler}/{assembly}/abundance/{assembly}.{fc}.tsv",
        gene_tax = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/gene_taxonomy.tsv",
    output:
        expand("results/collated/co-assembly/{{assembler}}/{{assembly}}/abundance_taxonomy/{tax_rank}.{tax_name}.{{assembly}}.{{fc}}.tsv",
               tax_rank = config["tax_rank"], tax_name = config["tax_name"])
    run:
        abundance = pd.read_csv(input.abundance, sep="\t", index_col=0, header=0)
        orig_genes = abundance.shape[0]
        gene_tax = pd.read_csv(input.gene_tax, sep="\t", index_col=0, header=0)
        # Filter to taxonomy
        try:
            orfs = list(set(gene_tax.loc[gene_tax[config["tax_rank"]]==config["tax_name"]].index))
        except KeyError:
            orfs = []
        if len(orfs) == 0:
            shell("touch {output}")
        abundance = abundance.loc[set(orfs).intersection(set(abundance.index))]
        print("{}/{} genes kept".format(abundance.shape[0], orig_genes))
        abundance.to_csv(output[0], sep="\t", index=True)

#################################
## EGGNOG-MAPPER CO-ASSEMBLIES ##
#################################
rule emapper_search_co:
    input:
        "results/annotation/co-assembly/{assembler}/{assembly}/frame_selection/final.reformat.faa",
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
        "resources/eggnog/eggnog.db"
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
        cp {input[1]} /dev/shm/$SLURM_JOB_ID
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
            db = ["enzymes","kos","modules","pathways","tc","cazy"])
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
        tmpdir = os.path.join(os.path.expandvars("$TMPDIR"), "{assembly}")
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
    """Normalize modules/pathways further by the number of KEGG orthologs belonging to each module/pathway"""
    input:
        abundance = expand("results/annotation/co-assembly/{{assembler}}/{{assembly}}/featureCounts/{sample_id}.{{fc}}.tab",
            sample_id = samples.keys()),
        parsed = "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/{db}.parsed.tab",
        norm = "resources/kegg/kegg_ko2{db}.tsv"
    output:
        "results/collated/co-assembly/{assembler}/{assembly}/eggNOG/{db}.norm.{fc}.tsv"
    params:
        src = "source/utils/eggnog-parser.py",
        tmpdir = os.path.join(os.path.expandvars("$TMPDIR"), "{assembly}", "norm")
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
        gene_tax = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/gene_taxonomy.tsv",
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

rule contigtax_search_co:
    """
    Runs blastx of contigs against protein database and reports best LCA
    """
    input:
        db = "resources/diamond/diamond.dmnd",
        fa = "results/co-assembly/{assembler}/{assembly}/final.fa"
    output:
        diamond = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/diamond.tsv.gz",
        logfile = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/diamond.log"
    params:
        tmp_out = "$TMPDIR/{assembly}/diamond.tsv",
        tmp_log = "$TMPDIR/{assembly}/diamond.log",
        tmp_dir = "$TMPDIR/{assembly}"
    threads: 20
    conda: "../../envs/contigtax.yaml"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*10
    shell:
        """
        mkdir -p {params.tmp_dir}
        contigtax search -m blastx -p {threads} --top 10 -e 0.001 -b 20 -c 1 \
            -t {params.tmp_dir} {input.fa} {input.db} {params.tmp_out}
        mv {params.tmp_out}.gz {output.diamond}
        mv {params.tmp_log} {output.logfile}
        """

rule assign_taxonomy_co:
    input:
        diamond = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/diamond.tsv.gz",
        taxdb = ancient("resources/taxonomy/taxonomy.sqlite")
    output:
        tax = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/contig_taxonomy.tsv"
    log: "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/taxonomy.log"
    params:
        taxdir = "resources/taxonomy",
        reportranks = "superkingdom kingdom phylum class order family genus species"
    threads: 10
    conda: "../../envs/contigtax.yaml"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    shell:
        """
        contigtax assign -t {params.taxdir} --reportranks {params.reportranks} -p {threads} \
            -T 5 {input.diamond} {output.tax} 2>{log}
        """

rule create_blob_co:
    input:
        hit = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/contig_hits.tsv",
        bam = "results/map/co-assembly/{assembler}/{assembly}/{sample_id}.bam",
        fa = "results/co-assembly/{assembler}/{assembly}/final.fa",
        names = "resources/taxonomy/names.dmp",
        nodes = "resources/taxonomy/nodes.dmp"
    output:
        blob = temp("results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/blobtools/{sample_id}.blobDB.json")
    params:
        prefix = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/blobtools/{sample_id}"
    conda: "../../envs/blobtools.yaml"
    shell:
        """
        blobtools create --nodes {input.nodes} --names {input.names} -i {input.fa} -b {input.bam} -t {input.hit} \ 
            --title {wildcards.sample_id} -o {params.prefix}
        """

rule plot_blob_co:
    input:
        blob = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/blobtools/{sample_id}.blobDB.json"
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/blobtools/{sample_id}.bestsum.{rank}.p{num}.span.100.exclude_other.blobplot.bam0.png",
        "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/blobtools/{sample_id}.bestsum.{rank}.p{num}.span.100.exclude_other.blobplot.read_cov.bam0.png",
        "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/blobtools/{sample_id}.bestsum.{rank}.p{num}.span.100.exclude_other.blobplot.stats.txt"
    params:
        prefix = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/blobtools/",
    conda: "../../envs/blobtools.yaml"
    shell:
        """
        blobtools blobplot -i {input.blob} --sort_first 'no-hit,undef' -o {params.prefix} -p {wildcards.num} \
            -r {wildcards.rank} --exclude other
        """

rule transfer_taxonomy_co:
    input:
        contig_tax = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/contig_taxonomy.tsv",
        gff = "results/annotation/co-assembly/{assembler}/{assembly}/frame_selection/final.reformat.gff"
    output:
        gene_tax = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/gene_taxonomy.tsv"
    run:
        gff = pd.read_table(input.gff, usecols = [0,8], names=["contig","gene"], header=None, index_col=1)
        contig_tax = pd.read_table(input.contig_tax, header=0, index_col=0)
        gff.rename(index=lambda x: x.split("=")[-1], inplace=True)
        gene_tax = pd.merge(gff, contig_tax, left_on="contig", right_index=True, how="left")
        gene_tax.drop("contig", axis=1, inplace=True)
        gene_tax.fillna("Unclassified", inplace=True)
        gene_tax.to_csv(output.gene_tax, sep="\t")

rule sum_taxonomy_co:
    input:
        gene_tax = "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/gene_taxonomy.tsv",
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
        fasta = "results/annotation/co-assembly/{assembler}/{assembly}/frame_selection/final.reformat.faa",
        db = "resources/dbCAN/dbCAN-fam-HMMs.txt",
        pressed = expand("resources/dbCAN/dbCAN-fam-HMMs.txt.{suffix}", suffix = ["h3f","h3i","h3m","h3p"])
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/dbCAN/dbCAN.out.dm",
        "results/annotation/co-assembly/{assembler}/{assembly}/dbCAN/dbCAN.out"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
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
