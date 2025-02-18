localrules:
    mmseqs_convertali,
    parse_mmseqs_first,
    firstpass_fungal_proteins,
    secondpass_fungal_proteins,
    collate_taxonomy,
    eggnog_merge_and_sum,
    parse_featurecounts,
    parse_eggnog,
    quantify_eggnog,
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
    output:
        expand("results/annotation/{{assembler}}/{{sample_id}}/transdecoder/final.fa.transdecoder_dir/longest_orfs.{suff}",
            suff = ["pep","gff3","cds"])
    input:
        fa = "results/assembly/{assembler}/{sample_id}/final.fa"
    log:
        "results/annotation/{assembler}/{sample_id}/transdecoder/longest_orfs.log"
    container: "docker://trinityrnaseq/transdecoder:5.7.1"
    params:
        output_dir = lambda wildcards: f"results/annotation/{wildcards.assembler}/{wildcards.sample_id}/transdecoder"
    shell:
        """
        TransDecoder.LongOrfs -t {input.fa} -O {params.output_dir} > {log} 2>&1
        """

rule mmseqs_createquerydb_longorfs:
    """
    Create mmseqs database from longest ORFs
    """
    output:
        "results/annotation/{assembler}/{sample_id}/transdecoder/longorfs-queryDB"
    input:
        query = rules.transdecoder_longorfs.output[0],
    log:
        "results/annotation/{assembler}/{sample_id}/transdecoder/longorfs-queryDB.log"
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
    output:
        tax = expand("results/annotation/{{assembler}}/{{sample_id}}/transdecoder/{{td_db}}-taxaDB.{suff}", suff = ["index","dbtype"]),
        aln = expand("results/annotation/{{assembler}}/{{sample_id}}/transdecoder/{{td_db}}-taxaDB_aln/first.{suff}", suff = ["index","dbtype"]),
    input:
        query = rules.mmseqs_createquerydb_longorfs.output,
        target = os.path.join(config["mmseqs_db_dir"], "{td_db}"),
    log:
        "results/annotation/{assembler}/{sample_id}/transdecoder/mmseqs_firstpass_taxonomy.{td_db}.log"
    threads: 4
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    conda: "../../envs/mmseqs.yaml"
    params:
        tmp = lambda wildcards: f"{os.environ.get('TMPDIR', 'scratch')}/mmseqs_firstpass_taxonomy.{wildcards.assembler}.{wildcards.sample_id}.{wildcards.td_db}", 
        split_memory_limit = lambda wildcards, resources: int(resources.mem_mb*.8),
        output = lambda wildcards, output: f"{os.path.dirname(output[0])}/{wildcards.td_db}-taxaDB",
        ranks = "superkingdom,kingdom,phylum,class,order,family,genus,species",
        aln_dir = lambda wildcards, output: os.path.dirname(output.aln[0])
    shell:
        """
        mkdir -p {params.tmp}
        mmseqs taxonomy {input.query} {input.target} {params.output} {params.tmp} -e 1e-5 -a 1 \
            --lca-mode 3 --tax-output-mode 0 --lca-ranks {params.ranks} --tax-lineage 1 \
            --split-memory-limit {params.split_memory_limit}M --threads {threads} > {log} 2>&1
        mv {params.tmp}/latest/first.* {params.aln_dir}
        rm -r {params.tmp}
        """

rule mmseqs_convertali:
    """
    Convert mmseqs alignment to m8 format for use with transdecoder
    """
    output:
        m8 = "results/annotation/{assembler}/{sample_id}/transdecoder/{td_db}.m8"
    input:
        alignment = rules.mmseqs_firstpass_taxonomy.output.aln,
        target = os.path.join(config["mmseqs_db_dir"], "{td_db}"),
        query = rules.mmseqs_createquerydb_longorfs.output,
    log:
        "results/annotation/{assembler}/{sample_id}/transdecoder/{td_db}.convertali.log"
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
    output:
        expand("results/annotation/{{assembler}}/{{sample_id}}/transdecoder/final.fa.transdecoder.{suff}", 
            suff = ["pep","gff3","cds","bed"])
    input:
        fa = "results/assembly/{assembler}/{sample_id}/final.fa",
        m8 = expand("results/annotation/{{assembler}}/{{sample_id}}/transdecoder/{td_db}.m8", td_db = config["transdecoder_homology_db"])
    log:
        "results/annotation/{assembler}/{sample_id}/transdecoder/predict.log"
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
    output:
        tsv = "results/annotation/{assembler}/{sample_id}/taxonomy/{td_db}/firstpass.tsv"
    input:
        query = rules.mmseqs_createquerydb_longorfs.output,
        result = rules.mmseqs_firstpass_taxonomy.output.tax
    log:
        "results/annotation/{assembler}/{sample_id}/taxonomy/{td_db}/mmseqs_createtsv_first.log"
    container: "docker://quay.io/biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0"
    conda: "../../envs/mmseqs.yaml"
    params:
        result = lambda wildcards, input: f"{os.path.dirname(input.result[0])}/{wildcards.td_db}-taxaDB",
    threads: 1
    shell:
        """
        mmseqs createtsv {input.query} {params.result} {output.tsv} --threads {threads} > {log} 2>&1
        """

rule parse_mmseqs_first:
    """
    Parse mmseqs taxonomy output, adds columns with ranks
    """
    output:
        tsv = "results/annotation/{assembler}/{sample_id}/taxonomy/{td_db}/firstpass.parsed.tsv"
    input:
        tsv = rules.mmseqs_createtsv_first.output,
    log:
        "results/annotation/{assembler}/{sample_id}/taxonomy/{td_db}/firstpass.parsed.log"
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
    output:
        "results/annotation/{assembler}/{sample_id}/taxonomy/{td_db}/firstpass.fungal.faa"
    input:
        parsed = rules.parse_mmseqs_first.output.tsv,
        fa = rules.transdecoder_predict.output[0]
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
    output:
        expand("results/annotation/{{assembler}}/{{sample_id}}/taxonomy/{mmseqs_db}/secondpass-taxresult_{suff}", suff = ["lca.tsv", "report","tophit_aln","tophit_report"], mmseqs_db = config["mmseqs_db"])
    input:
        query = expand("results/annotation/{{assembler}}/{{sample_id}}/taxonomy/{td_db}/firstpass.fungal.faa",td_db = config["transdecoder_homology_db"]),
        target = get_mmseq_taxdb,
    log:
        expand("results/annotation/{{assembler}}/{{sample_id}}/taxonomy/{mmseqs_db}/mmseqs_secondpass_taxonomy.log", mmseqs_db = config["mmseqs_db"])
    params:
        output = lambda wildcards, output: f"{os.path.dirname(output[0])}/secondpass-taxresult",
        tmp = lambda wildcards: f"{os.environ.get('TMPDIR', 'scratch')}/mmseqs_secondpass_taxonomy_co.{wildcards.assembler}.{wildcards.sample_id}", 
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

rule parse_mmseqs_second:
    """
    Parse mmseqs taxonomy output from second run
    """
    output:
        tsv = "results/annotation/{assembler}/{sample_id}/taxonomy/{mmseqs_db}/secondpass.parsed.tsv"
    input:
        taxresult = rules.mmseqs_secondpass_taxonomy.output,
    log:
        "results/annotation/{assembler}/{sample_id}/taxonomy/{mmseqs_db}/secondpass.parsed.log"
    params:
        script = workflow.source_path("../utils/parse_mmseqs.py"),
        tsv = lambda wildcards, input: os.path.join(os.path.dirname(input.taxresult[0]), "secondpass-taxresult_lca.tsv"),
        ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    shell:
        """
        python {params.script} -i {params.tsv} -o {output.tsv} -r {params.ranks} > {log} 2>&1
        """

rule secondpass_fungal_proteins:
    output:
        expand("results/annotation/{{assembler}}/{{sample_id}}/genecall/fungal.{suff}", 
            suff = ["faa","gff3","cds","bed","taxonomy.tsv"])
    input:
        parsed = expand("results/annotation/{{assembler}}/{{sample_id}}/taxonomy/{mmseqs_db}/secondpass.parsed.tsv", mmseqs_db = config["mmseqs_db"]),
        genecall = rules.transdecoder_predict.output,
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        indir = lambda wildcards, input: os.path.dirname(input.genecall[0]),
    shadow: "minimal"
    run:
        write_fungal_proteins(input.parsed[0], params.indir, params.outdir)

###################
## READ COUNTING ##
###################
rule featurecount:
    """
    Count reads mapping to features using featureCounts
    """
    output:
        cnt = "results/annotation/{assembler}/{sample_id}/featureCounts/fc.tsv",
        summary = "results/annotation/{assembler}/{sample_id}/featureCounts/fc.tsv.summary"
    input:
        gff = "results/annotation/{assembler}/{sample_id}/transdecoder/final.fa.transdecoder.gff3",
        bam = "results/map/{assembler}/{sample_id}/{sample_id}.bam"
    log:
        "results/annotation/{assembler}/{sample_id}/featureCounts/fc.log"
    conda: "../../envs/featurecount.yaml"
    params:
        setting = config["fc_params"]
    container: "docker://quay.io/biocontainers/subread:2.0.8--h577a1d6_0"
    threads: 4
    shell:
        """
        featureCounts -T {threads} {params.setting} -a {input.gff} -o {output.cnt} {input.bam} > {log} 2>&1
        """

rule parse_featurecounts:
    """
    Parses featureCounts output
    """
    output:
        "results/annotation/{assembler}/{sample_id}/featureCounts/fc.raw.tsv"
    input:
        fc="results/annotation/{assembler}/{sample_id}/featureCounts/fc.tsv",
    run:
        df = pd.read_csv(input.fc, sep="\t", index_col=0, header=0, comment="#", usecols=[0,6], names=["gene_id", wildcards.sample_id])
        df.to_csv(output[0], sep="\t")

###################
## EGGNOG-MAPPER ##
###################
rule emapper_search:
    """
    Run eggNOG-mapper search on fungal proteins
    """
    output:
        "results/annotation/{assembler}/{sample_id}/eggNOG/annotation_results.emapper.seed_orthologs",
    input:
        faa="results/annotation/{assembler}/{sample_id}/genecall/fungal.faa",
        db=f"{config['emapper_db_dir']}/eggnog.db",
        mmseqs_db=f"{config['emapper_db_dir']}/mmseqs/mmseqs.db"
    params:
        data_dir = lambda wc, input: os.path.dirname(input.db),
        out = "annotation_results",
        tmpdir = "$TMPDIR/eggnog/{assembler}/{sample_id}",
        tmp_out = "$TMPDIR/eggnog/{assembler}/{sample_id}/annotation_results",
        outdir = lambda wc, output: os.path.dirname(output[0]),
    log: "results/annotation/{assembler}/{sample_id}/eggNOG/emapper.log"
    threads: 10
    #shadow: "shallow"
    conda: "../../envs/emapper.yaml"
    container: "docker://quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"
    shell:
        """
        mkdir -p {params.tmpdir}
        mkdir -p {params.outdir}
        emapper.py -m mmseqs --mmseqs_db {input.mmseqs_db} --no_annot --no_file_comments --itype proteins --cpu {threads} \
            -i {input[0]} -o {params.out} --temp_dir {params.tmpdir} \
            --output_dir {params.tmpdir} --data_dir {params.data_dir} 2>{log}
        mv {params.tmp_out}.emapper.seed_orthologs {output[0]}
        """

rule emapper_annotate_hits:
    """
    Annotate hits from eggNOG-mapper search
    """
    output:
        "results/annotation/{assembler}/{sample_id}/eggNOG/annotation_results.emapper.annotations"
    input:
        seed_orthologs=rules.emapper_search.output,
        db=f"{config['emapper_db_dir']}/eggnog.db",
        mmseqs_db=f"{config['emapper_db_dir']}/mmseqs/mmseqs.db"
    params:
        data_dir = lambda wc, input: os.path.dirname(input.db),
        tmpdir = "$TMPDIR/eggnog/{assembler}/{sample_id}",
        out = "results/annotation/{assembler}/{sample_id}/eggNOG/annotation_results",
    log: "results/annotation/{assembler}/{sample_id}/eggNOG/emapper.annotate.log"
    threads: 10
    conda: "../../envs/emapper.yaml"
    container: "docker://quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"
    resources:
        mem_mb = 50000
    shell:
        """
        emapper.py -m no_search --cpu {threads} --annotate_hits_table {input.seed_orthologs} -o {params.out} \
            --data_dir {params.data_dir} --dbmem 2>{log}
        """

##################
## PARSE EGGNOG ##
##################
rule parse_eggnog:
    """
    Parse eggNOG-mapper results and output tables for different annotation types
    """
    output:
        expand("results/annotation/{{assembler}}/{{sample_id}}/eggNOG/{db}.parsed.tsv",
            db = ["enzymes","kos","modules","pathways","tc","cazy"])
    input:
        f = rules.emapper_annotate_hits.output,
        db = expand("resources/kegg/{f}",
            f = ["kegg_ec2pathways.tsv","kegg_ko2ec.tsv",
                 "kegg_ko2pathways.tsv","kegg_kos.tsv","kegg_modules.tsv","kegg_pathways.tsv"])
    log: "results/annotation/{assembler}/{sample_id}/eggNOG/parser.log"
    params:
        src = workflow.source_path("../utils/eggnog-parser.py"),
        dldir = "resources/kegg",
        outdir = lambda wc, output: os.path.dirname(output[0])
    shell:
        """
        python {params.src} parse {params.dldir} {input.f} {params.outdir} 2>{log}
        """

rule quantify_eggnog:
    """
    Sum abundances of annotation types
    """
    output:
        "results/annotation/{assembler}/{sample_id}/eggNOG/{db}.raw.tsv"
    input:
        abundance = "results/annotation/{assembler}/{sample_id}/featureCounts/fc.raw.tsv",
        parsed = "results/annotation/{assembler}/{sample_id}/eggNOG/{db}.parsed.tsv"
    params:
        src = workflow.source_path("../utils/eggnog-parser.py"),
    shell:
        """
        python {params.src} quantify {input.abundance} {input.parsed} {output[0]}
        """

rule quantify_eggnog_normalized:
    """
    Normalize pathways further by the number of KEGG orthologs belonging to each pathway
    """
    output:
        "results/annotation/{assembler}/{sample_id}/eggNOG/{db}.norm.raw.tsv"
    input:
        abundance = "results/annotation/{assembler}/{sample_id}/featureCounts/fc.raw.tsv",
        parsed = "results/annotation/{assembler}/{sample_id}/eggNOG/{db}.parsed.tsv",
        norm = "resources/kegg/kegg_ko2{db}.tsv"
    params:
        src = workflow.source_path("../utils/eggnog-parser.py"),
    shell:
        """
        python {params.src} quantify --normalize {input.norm} {input.abundance} {input.parsed} {output[0]}
        """

rule eggnog_merge_and_sum:
    """
    Take count tables for an annotation type in each sample and merge them
    """
    output:
        "results/collated/{assembler}/eggNOG/{db}.raw.tsv"
    input:
        expand("results/annotation/{{assembler}}/{sample_id}/eggNOG/{{db}}.raw.tsv",
            sample_id = samples.keys())
    params:
        src = workflow.source_path("../utils/eggnog-parser.py"),
    shell:
        """
        python {params.src} merge --sum {input} {output[0]}
        """

##########################
## TAXONOMIC ANNOTATION ##
##########################

rule sum_taxonomy:
    """
    Sum taxonomic counts for fungal proteins
    """
    output:
        raw = "results/annotation/{assembler}/{sample_id}/taxonomy/taxonomy.raw.tsv"
    input:
        gene_tax = "results/annotation/{assembler}/{sample_id}/genecall/fungal.taxonomy.tsv",
        raw = "results/annotation/{assembler}/{sample_id}/featureCounts/fc.raw.tsv"
    run:
        ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
        gene_tax = pd.read_csv(input.gene_tax, header=0, sep="\t", index_col=0)
        raw = pd.read_csv(input.raw, header=0, sep="\t", index_col=0)
        gene_tax_raw = pd.merge(gene_tax, raw, left_index = True, right_index = True)
        species_raw = gene_tax_raw.groupby(ranks).sum().reset_index()
        # Write results
        species_raw.to_csv(output.raw, sep="\t", index=False)


def merge_tax(files):
    """
    Merge taxonomic count tables

    Parameters
    ----------
    files : list
        List of file paths

    Returns
    -------
    pandas.DataFrame
    """
    df = pd.DataFrame()
    for f in files:
        _df = pd.read_csv(f, header=0, sep="\t")
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
    """
    Collate taxonomic counts for all samples
    """
    output:
        raw = "results/collated/{assembler}/taxonomy/taxonomy.raw.tsv"
    input:
        raw = expand("results/annotation/{assembler}/{sample_id}/taxonomy/taxonomy.raw.tsv",
            assembler = config["assembler"], sample_id = samples.keys())
    run:
        raw = merge_tax(input.raw)
        raw.to_csv(output.raw, sep="\t", index=False)
