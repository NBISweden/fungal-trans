localrules:
    mmseqs_convertali,
    parse_mmseqs_first,
    firstpass_fungal_proteins,
    secondpass_fungal_proteins,
    collate_taxonomy,
    eggnog_merge_and_sum,
    #parse_featurecounts,
    parse_eggnog,
    #quantify_eggnog,
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
        expand(
            "results/annotation/{{assembler}}/{{sample_id}}/transdecoder/final.fa.transdecoder_dir/longest_orfs.{suff}",
            suff=["pep", "gff3", "cds"],
        ),
    input:
        fa="results/assembly/{assembler}/{sample_id}/final.fa",
    log:
        "results/annotation/{assembler}/{sample_id}/transdecoder/longest_orfs.log",
    container:
        "docker://trinityrnaseq/transdecoder:5.7.1"
    params:
        output_dir=lambda wildcards: f"results/annotation/{wildcards.assembler}/{wildcards.sample_id}/transdecoder",
    shell:
        """
        TransDecoder.LongOrfs -t {input.fa} -O {params.output_dir} > {log} 2>&1
        """


rule mmseqs_createquerydb_longorfs:
    """
    Create mmseqs database from longest ORFs
    """
    output:
        "results/annotation/{assembler}/{sample_id}/transdecoder/longorfs-queryDB",
    input:
        query=rules.transdecoder_longorfs.output[0],
    log:
        "results/annotation/{assembler}/{sample_id}/transdecoder/longorfs-queryDB.log",
    container:
        "docker://quay.io/biocontainers/mmseqs2:17.b804f--hd6d6fdc_0"
    conda:
        "../../envs/mmseqs.yaml"
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
        tax=expand(
            "results/annotation/{{assembler}}/{{sample_id}}/transdecoder/{{td_db}}-taxaDB.{suff}",
            suff=["index", "dbtype"],
        ),
        aln=expand(
            "results/annotation/{{assembler}}/{{sample_id}}/transdecoder/{{td_db}}-taxaDB_aln/first.{suff}",
            suff=["index", "dbtype"],
        ),
    input:
        query=rules.mmseqs_createquerydb_longorfs.output,
        target=os.path.join(config["mmseqs_db_dir"], "{td_db}"),
    log:
        "results/annotation/{assembler}/{sample_id}/transdecoder/mmseqs_firstpass_taxonomy.{td_db}.log",
    threads: 4
    container: "docker://quay.io/biocontainers/mmseqs2:17.b804f--hd6d6fdc_0"
    conda: "../../envs/mmseqs.yaml"
    shadow: "shallow"
    params:
        tmp=lambda wildcards: f"mmseqs_firstpass_taxonomy.{wildcards.assembler}.{wildcards.sample_id}.{wildcards.td_db}",
        split_memory_limit=lambda wildcards, resources: int(resources.mem_mb * 0.8),
        output=lambda wildcards, output: f"{os.path.dirname(output[0])}/{wildcards.td_db}-taxaDB",
        ranks="superkingdom,kingdom,phylum,class,order,family,genus,species",
        aln_dir=lambda wildcards, output: os.path.dirname(output.aln[0]),
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
        m8="results/annotation/{assembler}/{sample_id}/transdecoder/{td_db}.m8",
    input:
        alignment=rules.mmseqs_firstpass_taxonomy.output.aln,
        target=os.path.join(config["mmseqs_db_dir"], "{td_db}"),
        query=rules.mmseqs_createquerydb_longorfs.output,
    log:
        "results/annotation/{assembler}/{sample_id}/transdecoder/{td_db}.convertali.log",
    container:
        "docker://quay.io/biocontainers/mmseqs2:17.b804f--hd6d6fdc_0"
    conda:
        "../../envs/mmseqs.yaml"
    params:
        alignment=lambda wildcards, input: os.path.splitext(input.alignment[0])[0],
    shell:
        """
        mmseqs convertalis {input.query} {input.target} {params.alignment} {output.m8} --threads 1 --format-mode 0 > {log} 2>&1
        """


rule transdecoder_predict:
    """
    Run transdecoder predict using homology search results for ORF retention criteria
    """
    output:
        expand(
            "results/annotation/{{assembler}}/{{sample_id}}/transdecoder/final.fa.transdecoder.{suff}",
            suff=["pep", "gff3", "cds", "bed"],
        ),
    input:
        fa="results/assembly/{assembler}/{sample_id}/final.fa",
        m8=expand(
            "results/annotation/{{assembler}}/{{sample_id}}/transdecoder/{td_db}.m8",
            td_db=config["transdecoder_homology_db"],
        ),
    log:
        "results/annotation/{assembler}/{sample_id}/transdecoder/predict.log",
    container:
        "docker://trinityrnaseq/transdecoder:5.7.1"
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),
    shell:
        """
        TransDecoder.Predict -t {input.fa} --retain_blastp_hits {input.m8} -O {params.output_dir} > {log} 2>&1
        sed -i 's/*$//g' {output[0]}
        """


rule mmseqs_createtsv_first:
    """
    Create tsv file from mmseqs taxonomy output
    """
    output:
        tsv="results/annotation/{assembler}/{sample_id}/taxonomy/{td_db}/firstpass.tsv",
    input:
        query=rules.mmseqs_createquerydb_longorfs.output,
        result=rules.mmseqs_firstpass_taxonomy.output.tax,
    log:
        "results/annotation/{assembler}/{sample_id}/taxonomy/{td_db}/mmseqs_createtsv_first.log",
    container:
        "docker://quay.io/biocontainers/mmseqs2:17.b804f--hd6d6fdc_0"
    conda:
        "../../envs/mmseqs.yaml"
    params:
        result=lambda wildcards, input: f"{os.path.dirname(input.result[0])}/{wildcards.td_db}-taxaDB",
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
        tsv="results/annotation/{assembler}/{sample_id}/taxonomy/{td_db}/firstpass.parsed.tsv",
    input:
        tsv=rules.mmseqs_createtsv_first.output,
    log:
        "results/annotation/{assembler}/{sample_id}/taxonomy/{td_db}/firstpass.parsed.log",
    params:
        script=workflow.source_path("../utils/parse_mmseqs.py"),
        ranks=[
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
    shell:
        """
        python {params.script} -i {input.tsv} -o {output.tsv} -r {params.ranks} > {log} 2>&1
        """


rule firstpass_fungal_proteins:
    """
    Extract proteins from first taxonomy pass that are classified as fungi
    """
    output:
        "results/annotation/{assembler}/{sample_id}/taxonomy/{td_db}/firstpass.fungal.faa",
    input:
        parsed=rules.parse_mmseqs_first.output.tsv,
        fa=rules.transdecoder_predict.output[0],
    params:
        idlist=lambda wildcards, input: f"{os.path.dirname(input.parsed)}/fungal.ids",
    run:
        import pandas as pd

        df = pd.read_csv(input.parsed, sep="\t", header=0, index_col=0)
        fungi = df.loc[df["kingdom"] == "Fungi"].index
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
        expand(
            "results/annotation/{{assembler}}/{{sample_id}}/taxonomy/{mmseqs_db}/secondpass-taxresult_{suff}",
            suff=["lca.tsv", "report", "tophit_aln", "tophit_report"],
            mmseqs_db=config["mmseqs_db"],
        ),
    input:
        query=expand(
            "results/annotation/{{assembler}}/{{sample_id}}/taxonomy/{td_db}/firstpass.fungal.faa",
            td_db=config["transdecoder_homology_db"],
        ),
        target=get_mmseq_taxdb,
    log:
        expand(
            "results/annotation/{{assembler}}/{{sample_id}}/taxonomy/{mmseqs_db}/mmseqs_secondpass_taxonomy.log",
            mmseqs_db=config["mmseqs_db"],
        ),
    params:
        output=lambda wildcards, output: f"{os.path.dirname(output[0])}/secondpass-taxresult",
        tmp=lambda wildcards: f"{os.environ.get('TMPDIR', 'scratch')}/mmseqs_secondpass_taxonomy.{wildcards.assembler}.{wildcards.sample_id}",
        ranks="superkingdom,kingdom,phylum,class,order,family,genus,species",
        split_memory_limit=lambda wildcards, resources: int(resources.mem_mb * 0.8),
        target=lambda wildcards, input: (input.target).replace("_taxonomy", ""),
    container:
        "docker://quay.io/biocontainers/mmseqs2:17.b804f--hd6d6fdc_0"
    conda:
        "../../envs/mmseqs.yaml"
    threads: 10
    shell:
        """
        mmseqs easy-taxonomy {input.query} {params.target} {params.output} {params.tmp} \
            --lca-ranks {params.ranks} --lca-mode 3 --tax-lineage 1 --split-memory-limit {params.split_memory_limit}M \
            --threads {threads} > {log} 2>&1
        """


rule parse_mmseqs_second:
    """
    Parse mmseqs taxonomy output from second run
    """
    output:
        tsv="results/annotation/{assembler}/{sample_id}/taxonomy/{mmseqs_db}/secondpass.parsed.tsv",
    input:
        taxresult=rules.mmseqs_secondpass_taxonomy.output,
    log:
        "results/annotation/{assembler}/{sample_id}/taxonomy/{mmseqs_db}/secondpass.parsed.log",
    params:
        script=workflow.source_path("../utils/parse_mmseqs.py"),
        tsv=lambda wildcards, input: os.path.join(
            os.path.dirname(input.taxresult[0]), "secondpass-taxresult_lca.tsv"
        ),
        ranks=[
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
    shell:
        """
        python {params.script} -i {params.tsv} -o {output.tsv} -r {params.ranks} > {log} 2>&1
        """


rule secondpass_fungal_proteins:
    output:
        expand(
            "results/annotation/{{assembler}}/{{sample_id}}/genecall/fungal.{suff}",
            suff=["faa", "gff3", "cds", "bed", "taxonomy.tsv"],
        ),
    input:
        parsed=expand(
            "results/annotation/{{assembler}}/{{sample_id}}/taxonomy/{mmseqs_db}/secondpass.parsed.tsv",
            mmseqs_db=config["mmseqs_db"],
        ),
        genecall=rules.transdecoder_predict.output,
    log:
        "results/annotation/{assembler}/{sample_id}/genecall/secondpass_fungal_proteins.log",
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        indir = lambda wildcards, input: os.path.dirname(input.genecall[0]),
        src=workflow.source_path("../utils/write_fungal_proteins.py"),
    shell:
        """
        python {params.src} -i {input.parsed} \
            -g {params.indir}/final.fa.transdecoder.gff3 \
            -b {params.indir}/final.fa.transdecoder.bed \
            -c {params.indir}/final.fa.transdecoder.cds \
            -p {params.indir}/final.fa.transdecoder.pep \
            -o {params.outdir} \
            -a Parent --prefix fungal 2>{log}
        """

###################
## READ COUNTING ##
###################
rule featurecount:
    """
    Count reads mapping to features using featureCounts
    """
    output:
        cnt="results/annotation/{assembler}/{sample_id}/featureCounts/fc.tsv",
        summary="results/annotation/{assembler}/{sample_id}/featureCounts/fc.tsv.summary",
    input:
        gff="results/annotation/{assembler}/{sample_id}/transdecoder/final.fa.transdecoder.gff3",
        bam="results/map/{assembler}/{sample_id}/RSEM/bowtie2.bam",
    log:
        "results/annotation/{assembler}/{sample_id}/featureCounts/fc.log",
    conda:
        "../../envs/featurecount.yaml"
    params:
        setting=config["fc_params"],
    container:
        "docker://quay.io/biocontainers/subread:2.0.8--h577a1d6_0"
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
        "results/annotation/{assembler}/{sample_id}/featureCounts/fc.raw.tsv",
    input:
        fc="results/annotation/{assembler}/{sample_id}/featureCounts/fc.tsv",
    run:
        df = pd.read_csv(
            input.fc,
            sep="\t",
            index_col=0,
            header=0,
            comment="#",
            usecols=[0, 6],
            names=["gene_id", wildcards.sample_id],
        )
        df.to_csv(output[0], sep="\t")


rule interproscan:
    """
    Runs the interproscan nextflow pipeline to annotate proteins.
    """
    output:
        tsv="results/annotation/{assembler}/{sample_id}/interproscan/interproscan.tsv"
    input:
        input=rules.transdecoder_predict.output[0],
        datadir=rules.download_interproscan_data.output.data
    log:
        "results/annotation/{assembler}/{sample_id}/interproscan/interproscan.log",
    params:
        pipeline = "ebi-pf-team/interproscan6",
        revision = "6.0.0-beta",
        extra = lambda wildcards, resources: set_interproscan_extra(wildcards, resources),
        profile = config["interproscan_profiles"],
        outdir = lambda wildcards, output: os.path.dirname(output.tsv),
        outprefix = "interproscan"
    handover: True
    wrapper:
        "v7.2.0/utils/nextflow"

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
        #faa="results/annotation/{assembler}/{sample_id}/genecall/fungal.faa",
        faa=rules.transdecoder_predict.output[0],
        db=f"{config['emapper_db_dir']}/eggnog.db",
        mmseqs_db=f"{config['emapper_db_dir']}/mmseqs/mmseqs.db",
    params:
        data_dir=lambda wc, input: os.path.dirname(input.db),
        out="annotation_results",
        tmpdir="{assembler}/{sample_id}",
        tmp_out="{assembler}/{sample_id}/annotation_results",
        outdir=lambda wc, output: os.path.dirname(output[0]),
    log:
        "results/annotation/{assembler}/{sample_id}/eggNOG/emapper.log",
    threads: 10
    shadow: "shallow"
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
        "results/annotation/{assembler}/{sample_id}/eggNOG/annotation_results.emapper.annotations",
    input:
        seed_orthologs=rules.emapper_search.output,
        db=f"{config['emapper_db_dir']}/eggnog.db",
        mmseqs_db=f"{config['emapper_db_dir']}/mmseqs/mmseqs.db",
    params:
        data_dir=lambda wc, input: os.path.dirname(input.db),
        out="results/annotation/{assembler}/{sample_id}/eggNOG/annotation_results",
    log:
        "results/annotation/{assembler}/{sample_id}/eggNOG/emapper.annotate.log",
    threads: 10
    conda:
        "../../envs/emapper.yaml"
    container:
        "docker://quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"
    resources:
        mem_mb=50000,
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
    input:
        f=rules.emapper_annotate_hits.output,
        db=expand(
            "resources/kegg/{f}",
            f=[
                "kegg_ec2pathways.tsv",
                "kegg_ko2ec.tsv",
                "kegg_ko2pathways.tsv",
                "kegg_kos.tsv",
                "kegg_modules.tsv",
                "kegg_pathways.tsv",
            ],
        ),
    output:
        parsed = "results/annotation/{assembler}/{sample_id}/eggNOG/{db}.parsed.tsv",
    log:
        "results/annotation/{assembler}/{sample_id}/eggNOG/{db}.parsed.log",
    params:
        src=workflow.source_path("../utils/eggnog-parser.py"),
        dldir="resources/kegg",
        cmd=lambda wildcards: get_eggnog_parser_extra_cmd
    shell:
        """
        python {params.src} parse {input.f} {output.parsed} {params.cmd} > {log} 2>&1
        """

def get_quant_table(wildcards):
    if wildcards.quant_type in ["TPM","FPKM","expected_count"]:
        return f"results/map/{wildcards.assembler}/{wildcards.sample_id}/RSEM/RSEM.isoforms.results"
    elif wildcards.quant_type in ["tpm","est_counts"]:
        return f"results/map/{wildcards.assembler}/{wildcards.sample_id}/kallisto/abundance.tsv"
    
def get_protein_quant_table(wildcards):
    if wildcards.quant_type in ["TPM","FPKM","expected_count"]:
        return f"results/annotation/{wildcards.assembler}/{wildcards.sample_id}/abundance/RSEM/proteins.{wildcards.quant_type}.tsv"
    elif wildcards.quant_type in ["tpm","est_counts"]:
        return f"results/annotation/{wildcards.assembler}/{wildcards.sample_id}/abundance/kallisto/proteins.{wildcards.quant_type}.tsv"
    elif wildcards.quant_type=="raw":
        return f"results/annotation/{wildcards.assembler}/{wildcards.sample_id}/featureCounts/fc.raw.tsv"

rule sum_proteins:
    """
    Sum quantification to protein level
    """
    input:
        tsv=get_quant_table,
        gff="results/annotation/{assembler}/{sample_id}/transdecoder/final.fa.transdecoder.gff3"
    output:
        tsv="results/annotation/{assembler}/{sample_id}/abundance/{tool}/proteins.{quant_type}.tsv"
    run:
        import polars as pl
        gff_df = pl.scan_csv(input.gff, separator="\t", has_header=False, new_columns=["transcript_id","source","feature_type","start","end","score", "strand","frame","attr"])
        # generate transcript id to protein table
        transcript2protein = (
            gff_df.filter(pl.col("feature_type")=="CDS")
            .select(["transcript_id","attr"])
            .with_columns(protein=pl.col("attr").str.split(";").list.first().str.replace("ID=cds.",""))
            .drop("attr")
        ).collect(engine="streaming")
        tsv_df = pl.read_csv(input.tsv, separator="\t", ignore_errors=True)
        tsv_df = tsv_df.select([tsv_df.columns[0], wildcards.quant_type]).rename({tsv_df.columns[0]: "transcript_id"})
        joined = transcript2protein.join(tsv_df, on="transcript_id")
        protein_sum = joined.group_by("protein").agg(pl.sum(tsv_df.columns[1:]))
        protein_sum.write_csv(output.tsv, separator="\t")


rule quantify_eggnog:
    """
    Sum abundances of annotation types
    """
    input:
        abundance=get_protein_quant_table,
        parsed="results/annotation/{assembler}/{sample_id}/eggNOG/{db}.parsed.tsv",
    output:
        "results/annotation/{assembler}/{sample_id}/eggNOG/{db}.{quant_type}.tsv",
    run:
        abundance_df = pd.read_csv(input.abundance, sep="\t", index_col=0)
        abundance_df.columns = [wildcards.sample_id]
        parsed_df = pd.read_csv(input.parsed, sep="\t", index_col=0)
        df = pd.merge(
            parsed_df, abundance_df, left_index=True, right_index=True, how="right"
        )
        df.fillna("Unclassified", inplace=True)
        feature_cols = list(parsed_df.columns)
        df_sum = df.groupby(feature_cols).sum().reset_index()
        df_sum.set_index(feature_cols[0], inplace=True)
        df_sum.to_csv(output[0], sep="\t")

def merge_files(input):
    import polars as pl
    df = pl.DataFrame()
    for i, f in enumerate(sorted(input)):
        _df = pl.read_csv(f, separator="\t")
        if i==0:
            on = _df.select(pl.col(pl.String)).columns
            df = _df
            continue
        df = df.join(_df, on=on, how="full", coalesce=True).fill_null(0)
    return df

rule eggnog_merge_and_sum:
    """
    Take count tables for an annotation type in each sample and merge them
    """
    output:
        "results/collated/{assembler}/eggNOG/{db}.{quant_type}.tsv",
    input:
        expand(
            "results/annotation/{{assembler}}/{sample_id}/eggNOG/{{db}}.{{quant_type}}.tsv",
            sample_id=samples.keys(),
        ),
    run:
        df = merge_files(input)
        df.write_csv(output[0], separator="\t")

##########################
## TAXONOMIC ANNOTATION ##
##########################
rule sum_taxonomy:
    """
    Sum taxonomic counts for fungal proteins
    """
    output:
        quant="results/annotation/{assembler}/{sample_id}/taxonomy/taxonomy.{quant_type}.tsv",
    input:
        gene_tax="results/annotation/{assembler}/{sample_id}/genecall/fungal.taxonomy.tsv",
        quant=get_protein_quant_table
    run:
        ranks = [
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]
        gene_tax = pd.read_csv(input.gene_tax, header=0, sep="\t", index_col=0)
        quant = pd.read_csv(input.quant, header=0, sep="\t", index_col=0)
        quant.columns = [wildcards.sample_id]
        gene_tax_quant = pd.merge(gene_tax, quant, left_index=True, right_index=True)
        species_quant = gene_tax_quant.groupby(ranks).sum().reset_index()
        # Write results
        species_quant.to_csv(output.quant, sep="\t", index=False)

rule collate_taxonomy:
    """
    Collate taxonomic counts for all samples
    """
    output:
        "results/collated/{assembler}/taxonomy/taxonomy.{quant_type}.tsv",
    input:
        expand(
            "results/annotation/{assembler}/{sample_id}/taxonomy/taxonomy.{{quant_type}}.tsv",
            assembler=config["assembler"],
            sample_id=samples.keys(),
        ),
    run:
        df = merge_files(input)
        df.write_csv(output[0], separator="\t")
