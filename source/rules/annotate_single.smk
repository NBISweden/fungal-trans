localrules:
    collate_dbcan,
    collate_taxonomy,
    create_blob,
    dbcan_parse,
    eggnog_merge_and_sum,
    normalize_featurecount,
    parse_eggnog,
    plot_blob,
    quantify_eggnog,
    reformat_fasta,
    sum_dbcan,
    sum_taxonomy,
    transfer_taxonomy

##########################################
##          SINGLE ASSEMBLIES           ##
##########################################
# Below are rules for annotation of single assemblies

#####################
## FRAME SELECTION ##
#####################
rule frame_selection:
    input:
        fa = "results/assembly/{assembler}/{filter_source}/{sample_id}/final.fa"
    output:
        expand("results/annotation/{{assembler}}/{{filter_source}}/{{sample_id}}/frame_selection/final.{suffix}",
            suffix = ["gff","faa","fnn"])
    log:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/frame_selection/gms.log"
    params:
        tmpdir = "$TMPDIR/{assembler}/{filter_source}/{sample_id}",
        outdir = "results/annotation/{assembler}/{filter_source}/{sample_id}/frame_selection"
    resources:
        runtime=lambda wildcards, attempt: attempt*attempt*20
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

rule reformat_fasta:
    input:
        faa = "results/annotation/{assembler}/{filter_source}/{sample_id}/frame_selection/final.faa"
    output:
        faa = "results/annotation/{assembler}/{filter_source}/{sample_id}/frame_selection/final.reformat.faa",
        gff = "results/annotation/{assembler}/{filter_source}/{sample_id}/frame_selection/final.reformat.gff"
    params:
        sample_id = "{sample_id}"
    shell:
        """
        python source/utils/reformat_fasta.py {input.faa} {output.faa} {output.gff} {params.sample_id}
        """

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

rule contigtax_search:
    """
    Runs blastx of contigs against protein database and reports best LCA
    """
    input:
        db = "resources/diamond/diamond.dmnd",
        fa = "results/assembly/{assembler}/{filter_source}/{sample_id}/final.fa"
    output:
        diamond = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/diamond.tsv.gz",
        logfile = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/diamond.log"
    params:
        tmp_out = "$TMPDIR/{sample_id}/diamond.tsv",
        tmp_log = "$TMPDIR/{sample_id}/diamond.log",
        tmp_dir = "$TMPDIR/{sample_id}"
    threads: 20
    conda: "../../envs/contigtax.yaml"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    shell:
        """
        mkdir -p {params.tmp_dir}
        contigtax search -m blastx -p {threads} --top 10 -e 0.001 -b 20 -c 1 \
            -t {params.tmp_dir} {input.fa} {input.db} {params.tmp_out}
        mv {params.tmp_out}.gz {output.diamond}
        mv {params.tmp_log} {output.logfile}
        """

rule assign_taxonomy:
    """Use contigtax package to assign queries"""
    input:
        diamond = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/diamond.tsv.gz",
        taxdb = ancient("resources/taxonomy/taxonomy.sqlite")
    output:
        tax = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/contig_taxonomy.tsv"
    log: "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/taxonomy.log"
    threads: 10
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    params:
        taxdir = "resources/taxonomy",
        reportranks = "superkingdom kingdom phylum class order family genus species"
    conda: "../../envs/contigtax.yaml"
    shell:
        """
        contigtax assign -t {params.taxdir} --reportranks {params.reportranks} -p {threads} \
        -T 5 {input.diamond} {output.tax} 2>{log}
        """

rule create_blob:
    input:
        hit = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/contig_hits.tsv",
        bam = "results/map/{assembler}/{filter_source}/{sample_id}/{sample_id}.bam",
        fa = "results/assembly/{assembler}/{filter_source}/{sample_id}/final.fa",
        names = "resources/taxonomy/names.dmp",
        nodes = "resources/taxonomy/nodes.dmp"
    output:
        blob = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{sample_id}.blobDB.json"
    params:
        prefix = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{sample_id}"
    conda: "../../envs/blobtools.yaml"
    shell:
        """
        blobtools create --nodes {input.nodes} --names {input.names} -i {input.fa} -b {input.bam} -t {input.hit} \
         --title {wildcards.sample_id} -o {params.prefix}
        """

rule plot_blob:
    input:
        blob = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{sample_id}.blobDB.json"
    output:
        "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{sample_id}.bestsum.{rank}.p{num}.span.100.exclude_other.blobplot.bam0.png",
        "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{sample_id}.bestsum.{rank}.p{num}.span.100.exclude_other.blobplot.read_cov.bam0.png",
        "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/{sample_id}.bestsum.{rank}.p{num}.span.100.exclude_other.blobplot.stats.txt"
    params:
        prefix = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/"
    conda: "../../envs/blobtools.yaml"
    shell:
        """
        blobtools blobplot -i {input.blob} --sort_first 'no-hit,undef' -o {params.prefix} -p {wildcards.num} \
        -r {wildcards.rank} --exclude other
        """

rule transfer_taxonomy:
    input:
        contig_tax = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/contig_taxonomy.tsv",
        gff = "results/annotation/{assembler}/{filter_source}/{sample_id}/frame_selection/final.reformat.gff"
    output:
        gene_tax = "results/annotation/{assembler}/{filter_source}/{sample_id}/taxonomy/gene_taxonomy.tsv"
    run:
        gff = pd.read_table(input.gff, usecols = [0,8], names=["contig","gene"], header=None, index_col=1)
        contig_tax = pd.read_table(input.contig_tax, header=0, index_col=0)
        gff.rename(index=lambda x: x.split("=")[-1], inplace=True)
        gene_tax = pd.merge(gff, contig_tax, left_on="contig", right_index=True, how="left")
        gene_tax.drop("contig", axis=1, inplace=True)
        gene_tax.fillna("Unclassified", inplace=True)
        gene_tax.to_csv(output.gene_tax, sep="\t")

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
