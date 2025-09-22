localrules:
    #mmseqs_convertali_co,
    parse_mmseqs_first_co,
    parse_mmseqs_second_co,
    firstpass_fungal_proteins_co,
    #secondpass_fungal_proteins_co,
    collate_featurecount_co,
    #parse_featurecounts_co,
    parse_eggnog_co,
    #quantify_eggnog_co,
    sum_taxonomy_co,
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
    output:
        expand(
            "results/annotation/co-assembly/{{assembler}}/{{assembly}}/transdecoder/final.fa.transdecoder_dir/longest_orfs.{suff}",
            suff=["pep", "gff3", "cds"],
        ),
    input:
        fa="results/co-assembly/{assembler}/{assembly}/final.fa",
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/longest_orfs.log",
    container:
        "docker://trinityrnaseq/transdecoder:5.7.1"
    params:
        output_dir=lambda wildcards: f"results/annotation/co-assembly/{wildcards.assembler}/{wildcards.assembly}/transdecoder",
    resources:
        runtime=240,
    shell:
        """
        TransDecoder.LongOrfs -t {input.fa} -O {params.output_dir} > {log} 2>&1
        """


rule mmseqs_createquerydb_longorfs_co:
    """
    Create mmseqs database from longest ORFs
    """
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/longorfs-queryDB",
    input:
        query=rules.transdecoder_longorfs_co.output[0],
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/longorfs-queryDB.log",
    container:
        "docker://quay.io/biocontainers/mmseqs2:17.b804f--hd6d6fdc_0"
    conda:
        "../../envs/mmseqs.yaml"
    shell:
        """
        mmseqs createdb {input} {output} > {log} 2>&1
        """


rule mmseqs_firstpass_taxonomy_co:
    """
    Run mmseqs taxonomy on longest ORFs, also store the alignment results
    for use with transdecoder
    """
    output:
        tax=expand(
            "results/annotation/co-assembly/{{assembler}}/{{assembly}}/transdecoder/{{td_db}}-taxaDB.{suff}",
            suff=["index", "dbtype"],
        ),
        aln=expand(
            "results/annotation/co-assembly/{{assembler}}/{{assembly}}/transdecoder/{{td_db}}-taxaDB_aln/first.{suff}",
            suff=["index", "dbtype"],
        ),
    input:
        query=rules.mmseqs_createquerydb_longorfs_co.output,
        target=os.path.join(config["mmseqs_db_dir"], "{td_db}"),
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/mmseqs_firstpass_taxonomy_co.{td_db}.log",
    threads: 10
    container:
        "docker://quay.io/biocontainers/mmseqs2:17.b804f--hd6d6fdc_0"
    conda:
        "../../envs/mmseqs.yaml"
    params:
        tmp=lambda wildcards: f"{os.environ.get('TMPDIR', 'scratch')}/mmseqs_firstpass_taxonomy_co.{wildcards.assembler}.{wildcards.assembly}.{wildcards.td_db}",
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


rule mmseqs_convertali_co:
    """
    Convert mmseqs alignment to m8 format for transdecoder
    """
    output:
        m8="results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/{td_db}.m8",
    input:
        alignment=rules.mmseqs_firstpass_taxonomy_co.output.aln,
        target=os.path.join(config["mmseqs_db_dir"], "{td_db}"),
        query=rules.mmseqs_createquerydb_longorfs_co.output,
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/{td_db}.convertali.log",
    container:
        "docker://quay.io/biocontainers/mmseqs2:17.b804f--hd6d6fdc_0"
    conda:
        "../../envs/mmseqs.yaml"
    params:
        alignment=lambda wildcards, input: os.path.splitext(input.alignment[0])[0],
    shell:
        """
        mmseqs convertalis {input.query} {input.target} {params.alignment} {output.m8} --threads {threads} --format-mode 0 > {log} 2>&1
        """


rule transdecoder_predict_co:
    """
    Run transdecoder predict using homology search results for ORF retention criteria
    """
    output:
        expand(
            "results/annotation/co-assembly/{{assembler}}/{{assembly}}/transdecoder/final.fa.transdecoder.{suff}",
            suff=["pep", "gff3", "cds", "bed"],
        ),
    input:
        fa="results/co-assembly/{assembler}/{assembly}/final.fa",
        m8=expand(
            "results/annotation/co-assembly/{{assembler}}/{{assembly}}/transdecoder/{td_db}.m8",
            td_db=config["transdecoder_homology_db"],
        ),
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/predict.log",
    container:
        "docker://trinityrnaseq/transdecoder:5.7.1"
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),
    shell:
        """
        TransDecoder.Predict -t {input.fa} --retain_blastp_hits {input.m8} -O {params.output_dir} > {log} 2>&1
        sed -i 's/*$//g' {output[0]}
        """


rule mmseqs_createtsv_first_co:
    """
    Create tsv file from mmseqs taxonomy output
    """
    output:
        tsv="results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{td_db}/firstpass.tsv",
    input:
        query=rules.mmseqs_createquerydb_longorfs_co.output,
        result=rules.mmseqs_firstpass_taxonomy_co.output.tax,
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{td_db}/mmseqs_createtsv_first_co.log",
    container:
        "docker://quay.io/biocontainers/mmseqs2:17.b804f--hd6d6fdc_0"
    conda:
        "../../envs/mmseqs.yaml"
    params:
        result=lambda wildcards, input: f"{os.path.dirname(input.result[0])}/{wildcards.td_db}-taxaDB",
    threads: 1
    resources:
        mem_mb=1000,
    shell:
        """
        mmseqs createtsv {input.query} {params.result} {output.tsv} --threads {threads} > {log} 2>&1
        """


rule parse_mmseqs_first_co:
    """
    Parse mmseqs taxonomy output, adds columns with ranks
    """
    output:
        tsv="results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{td_db}/firstpass.parsed.tsv",
    input:
        tsv=rules.mmseqs_createtsv_first_co.output,
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{td_db}/firstpass.parsed.log",
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


rule firstpass_fungal_proteins_co:
    """
    Extract proteins from first taxonomy pass that are classified as fungi
    """
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{td_db}/firstpass.fungal.faa",
    input:
        parsed=rules.parse_mmseqs_first_co.output.tsv,
        fa=rules.transdecoder_predict_co.output[0],
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


rule mmseqs_secondpass_taxonomy_co:
    """
    Run mmseqs taxonomy on fungal proteins

    If the extra_genomes parameter is set in the config file, then the target in
    the input will be the custom mmseq database which is a concatenation of
    proteins from the extra genomes + fungal proteins in the database defined by
    the mmseqs_db parameter.
    """
    output:
        expand(
            "results/annotation/co-assembly/{{assembler}}/{{assembly}}/taxonomy/{mmseqs_db}/secondpass-taxresult_{suff}",
            suff=["lca.tsv", "report", "tophit_aln", "tophit_report"],
            mmseqs_db=config["mmseqs_db"],
        ),
    input:
        query=expand(
            "results/annotation/co-assembly/{{assembler}}/{{assembly}}/taxonomy/{td_db}/firstpass.fungal.faa",
            td_db=config["transdecoder_homology_db"],
        ),
        target=get_mmseq_taxdb,
    log:
        expand(
            "results/annotation/co-assembly/{{assembler}}/{{assembly}}/taxonomy/{mmseqs_db}/mmseqs_secondpass_taxonomy_co.log",
            mmseqs_db=config["mmseqs_db"],
        ),
    params:
        output=lambda wildcards, output: f"{os.path.dirname(output[0])}/secondpass-taxresult",
        tmp=lambda wildcards: f"{os.environ.get('TMPDIR', 'scratch')}/mmseqs_secondpass_taxonomy_co.{wildcards.assembler}.{wildcards.assembly}",
        ranks="superkingdom,kingdom,phylum,class,order,family,genus,species",
        split_memory_limit=lambda wildcards, resources: int(resources.mem_mb * 0.8),
        target=lambda wildcards, input: (input.target).replace("_taxonomy", ""),
    container:
        "docker://quay.io/biocontainers/mmseqs2:17.b804f--hd6d6fdc_0"
    conda:
        "../../envs/mmseqs.yaml"
    threads: 10
    resources:
        mem_mb=8000,
    shell:
        """
        mkdir -p {params.tmp}
        mmseqs easy-taxonomy {input.query} {params.target} {params.output} {params.tmp} \
            --lca-ranks {params.ranks} --lca-mode 3 --tax-lineage 1 --split-memory-limit {params.split_memory_limit}M \
            --threads {threads} > {log} 2>&1
        rm -r {params.tmp}
        """


rule parse_mmseqs_second_co:
    """
    Parse mmseqs taxonomy output from second run
    """
    output:
        tsv="results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{mmseqs_db}/secondpass.parsed.tsv",
    input:
        taxresult=rules.mmseqs_secondpass_taxonomy_co.output,
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/taxonomy/{mmseqs_db}/secondpass.parsed.log",
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


rule secondpass_fungal_proteins_co:
    """
    Output coding sequences, proteins, gff3, bed and taxonomy for fungal proteins
    """
    input:
        parsed=expand(
            "results/annotation/co-assembly/{{assembler}}/{{assembly}}/taxonomy/{mmseqs_db}/secondpass.parsed.tsv",
            mmseqs_db=config["mmseqs_db"],
        ),
        genecall=rules.transdecoder_predict_co.output,
    output:
        expand(
            "results/annotation/co-assembly/{{assembler}}/{{assembly}}/genecall/fungal.{suff}",
            suff=["faa", "gff3", "cds", "bed", "taxonomy.tsv"],
        ),
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/genecall/secondpass_fungal_proteins_co.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        indir=lambda wildcards, input: os.path.dirname(input.genecall[0]),
        src=workflow.source_path("../utils/write_fungal_proteins.py"),
    shadow:
        "shallow"
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

#################################
## READ COUNTING CO-ASSEMBLIES ##
#################################
rule featurecount_co:
    """
    Count read assignments to co-assemblies

    This rule counts assignments to all predicted genes, not only fungal ones.
    """
    input:
        gff="results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/final.fa.transdecoder.gff3",
        bam="results/map/co-assembly/{assembler}/{assembly}/{sample_id}.bam",
    output:
        cnt="results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.fc.tsv",
        summary="results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.fc.tsv.summary",
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.fc.log",
    resources:
        runtime=30,
        mem_mb=1000,
    conda:
        "../../envs/featurecount.yaml"
    container:
        "docker://quay.io/biocontainers/subread:2.0.8--h577a1d6_0"
    params:
        setting=config["fc_params"],
    threads: 4
    shell:
        """
        featureCounts -T {threads} {params.setting} -a {input.gff} -o {output.cnt} {input.bam} > {log} 2>&1
        """


rule parse_featurecounts_co:
    """
    Parses featureCounts output
    """
    input:
        fc="results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.fc.tsv",
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/featureCounts/{sample_id}.raw.tsv",
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


rule collate_featurecount_co:
    """
    Collate featureCounts output for all samples
    """
    input:
        expand(
            "results/annotation/co-assembly/{{assembler}}/{{assembly}}/featureCounts/{sample_id}.raw.tsv",
            sample_id=samples.keys(),
        ),
    output:
        "results/collated/co-assembly/{assembler}/{assembly}/abundance/featureCounts/raw.tsv",
    run:
        df = pd.DataFrame()
        for f in input:
            _df = pd.read_csv(f, sep="\t", header=0, index_col=0)
            df = pd.merge(df, _df, right_index=True, left_index=True, how="outer")
        df.to_csv(output[0], sep="\t", index=True, header=True)

################################
## INTERPROSCAN CO-ASSEMBLIES ##
################################
def set_interproscan_extra(wildcards, resources):
    extra = "-c config/interproscan.config"
    if "kebnekaise" in config["interproscan_profiles"] or "dardel" in config["interproscan_profiles"]:
        extra += f" --project {resources.slurm_account}"
    return extra

rule interproscan_co:
    """
    Runs the interproscan nextflow pipeline to annotate proteins.
    """
    output:
        tsv="results/annotation/co-assembly/{assembler}/{assembly}/interproscan/interproscan.tsv"
    input:
        input=rules.transdecoder_predict_co.output[0],
        datadir=rules.download_interproscan_data.output.data
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/interproscan/interproscan.log"
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

#################################
## EGGNOG-MAPPER CO-ASSEMBLIES ##
#################################
rule emapper_search_co:
    """
    Run eggNOG-mapper on all predicted proteins from co-assembly
    """
    input:
        faa=rules.transdecoder_predict_co.output[0],
        db=f"{config['emapper_db_dir']}/eggnog.db",
        mmseqs_db=f"{config['emapper_db_dir']}/mmseqs/mmseqs.db",
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/annotation_results.emapper.seed_orthologs",
    params:
        data_dir=lambda wc, input: os.path.dirname(input.db),
        out="annotation_results",
        tmpdir="{assembly}",
        tmp_out="{assembly}/annotation_results",
        outdir=lambda wc, output: os.path.dirname(output[0]),
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/emapper.log",
    threads: 10
    conda: "../../envs/emapper.yaml"
    container: "docker://quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"
    shadow: "shallow"
    shell:
        """
        mkdir -p {params.tmpdir}
        mkdir -p {params.outdir}
        emapper.py -m mmseqs --mmseqs_db {input.mmseqs_db} --no_annot --no_file_comments --itype proteins --cpu {threads} \
            -i {input[0]} -o {params.out} --temp_dir {params.tmpdir} \
            --output_dir {params.outdir} --data_dir {params.data_dir} 2>{log}
        """


rule emapper_annotate_hits_co:
    """
    Annotate hits from eggNOG-mapper
    """
    input:
        seed_orthologs=rules.emapper_search_co.output,
        db=f"{config['emapper_db_dir']}/eggnog.db",
        mmseqs_db=f"{config['emapper_db_dir']}/mmseqs/mmseqs.db",
    output:
        "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/annotation_results.emapper.annotations",
    params:
        data_dir=lambda wc, input: os.path.dirname(input.db),
        tmpdir="$TMPDIR/eggnog/{assembly}/",
        out="results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/annotation_results",
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/emapper.annotate.log",
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


################################
## PARSE EGGNOG CO-ASSEMBLIES ##
################################
rule parse_eggnog_co:
    """
    Parses the eggnog-mapper output and outputs tab-separated files for pathways, kegg orthologs, enzymes and modules.

    Each tab-separated file contains the ORF id in the first column. ORFs annotated to multiple features are duplicated in the output.

    For example, an ORF annotated like so in the emapper output:
    #query ... KEGG_ko              KEGG_Pathway                                                       ...
    ORF1   ... ko:K07901,ko:K18158 ko04144,ko04152,ko04530,ko04972,map04144,map04152,map04530,map04972 ...

    Will be output as:
    # kos.parsed.tsv
    ORF1     K18158  K18158  NCA2; nuclear control of ATPase protein 2
    ORF1     K07901  K07901  RAB8A, MEL; Ras-related protein Rab-8A
    # pathways.parsed.tsv
    ORF1     map04144        Endocytosis [PATH:ko04144]      09140 Cellular Processes        09141 Transport and catabolism
    ORF1     map04152        AMPK signaling pathway [PATH:ko04152]   09130 Environmental Information Processing      09132 Signal transduction
    ORF1     map04530        Tight junction [PATH:ko04530]   09140 Cellular Processes        09144 Cellular community - eukaryotes
    ORF1     map04972        Pancreatic secretion [PATH:ko04972]     09150 Organismal Systems        09154 Digestive system
    """
    input:
        f=rules.emapper_annotate_hits_co.output,
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
        expand(
            "results/annotation/co-assembly/{{assembler}}/{{assembly}}/eggNOG/{db}.parsed.tsv",
            db=["enzymes", "kos", "modules", "pathways", "tc", "cazy"],
        ),
    log:
        "results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/parser.log",
    params:
        src=workflow.source_path("../utils/eggnog-parser.py"),
        dldir="resources/kegg",
        outdir=lambda wc, output: os.path.dirname(output[0]),
    shell:
        """
        python {params.src} parse {params.dldir} {input.f} {params.outdir} 2>{log}
        """

def get_quant_table(wildcards):
    if wildcards.quant_type in ["TPM","FPKM","expected_count"]:
        return f"results/collated/co-assembly/{wildcards.assembler}/{wildcards.assembly}/abundance/RSEM/isoforms.{wildcards.quant_type}.tsv"
    elif wildcards.quant_type in ["tpm","est_counts"]:
        return f"results/collated/co-assembly/{wildcards.assembler}/{wildcards.assembly}/abundance/kallisto/{wildcards.quant_type}.tsv"
    
def get_protein_quant_table(wildcards):
    if wildcards.quant_type in ["TPM","FPKM","expected_count"]:
        return f"results/collated/co-assembly/{wildcards.assembler}/{wildcards.assembly}/abundance/RSEM/proteins.{wildcards.quant_type}.tsv"
    elif wildcards.quant_type in ["tpm","est_counts"]:
        return f"results/collated/co-assembly/{wildcards.assembler}/{wildcards.assembly}/abundance/kallisto/proteins.{wildcards.quant_type}.tsv"
    elif wildcards.quant_type=="raw":
            return f"results/collated/co-assembly/{wildcards.assembler}/{wildcards.assembly}/abundance/featureCounts/{wildcards.quant_type}.tsv"

rule sum_proteins_co:
    """
    Sum quantification to protein level
    """
    input:
        tsv=get_quant_table,
        gff="results/annotation/co-assembly/{assembler}/{assembly}/transdecoder/final.fa.transdecoder.gff3"
    output:
        tsv="results/collated/co-assembly/{assembler}/{assembly}/abundance/{tool}/proteins.{quant_type}.tsv"
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
        tsv_df = pl.read_csv(input.tsv, separator="\t")
        joined = transcript2protein.join(tsv_df, on="transcript_id")
        protein_sum = joined.group_by("protein").agg(pl.sum(tsv_df.columns[1:]))
        protein_sum.write_csv(output.tsv, separator="\t")


# TODO: Fix quantifying
rule quantify_eggnog_co:
    """
    Sums up quantification results for each feature in the eggNOG database
    Note that since ORFs can be annotated to multiple features, the same ORF can be counted multiple times.
    """
    input:
        abundance=get_protein_quant_table,
        parsed="results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/{db}.parsed.tsv",
    output:
        "results/collated/co-assembly/{assembler}/{assembly}/eggNOG/{db}.{quant_type}.tsv",
    run:
        abundance_df = pd.read_csv(input.abundance, sep="\t", index_col=0)
        parsed_df = pd.read_csv(input.parsed, sep="\t", index_col=0)
        df = pd.merge(
            parsed_df, abundance_df, left_index=True, right_index=True, how="right"
        )
        df.fillna("Unclassified", inplace=True)
        feature_cols = list(parsed_df.columns)
        df_sum = df.groupby(feature_cols).sum().reset_index()
        df_sum.set_index(feature_cols[0], inplace=True)
        df_sum.to_csv(output[0], sep="\t")

rule eggnog_tax_annotations:
    """Collate gene annotations and abundances for a certain rank:taxon combination"""
    input:
        gene_tax="results/annotation/co-assembly/{assembler}/{assembly}/genecall/fungal.taxonomy.tsv",
        parsed="results/annotation/co-assembly/{assembler}/{assembly}/eggNOG/{db}.parsed.tsv",
        abundance=get_protein_quant_table,
    output:
        expand(
            "results/collated/co-assembly/{{assembler}}/{{assembly}}/eggNOG_taxonomy/{tax_rank}.{tax_name}.{{db}}.{{quant_type}}.tsv",
            tax_rank=config["tax_rank"],
            tax_name=config["tax_name"],
        ),
    run:
        import pandas as pd
        abundance = pd.read_csv(input.abundance, header=0, sep="\t", index_col=0)
        parsed = pd.read_csv(input.parsed, header=0, sep="\t", index_col=0)
        gene_tax = pd.read_csv(input.gene_tax, header=0, sep="\t", index_col=0)
        features = len(set(parsed.loc[:, parsed.columns[0]]))
        orig_orfs = abundance.shape[0]
        try:
            orfs = list(
                set(
                    gene_tax.loc[
                        gene_tax[config["tax_rank"]] == config["tax_name"]
                    ].index
                )
            )
        except KeyError:
            orfs = []
        if len(orfs) == 0:
            shell("touch {output}")
        else:
            abundance = abundance.loc[
                list(set(orfs).intersection(set(abundance.index)))
            ]
            parsed = parsed.loc[list(set(orfs).intersection(set(parsed.index)))]
            tax_features = len(set(parsed.loc[:, parsed.columns[0]]))
            print(
                "{}/{} features, {}/{} orfs matched".format(
                    tax_features, features, len(orfs), orig_orfs
                )
            )
            df = pd.merge(parsed, abundance, left_index=True, right_index=True)
            df_sum = df.groupby(list(parsed.columns)).sum().reset_index()
            df_sum.to_csv(output[0], sep="\t", index=False, header=True)


##########################
## TAXONOMIC ANNOTATION ##
##########################

rule sum_taxonomy_co:
    input:
        gene_tax="results/annotation/co-assembly/{assembler}/{assembly}/genecall/fungal.taxonomy.tsv",
        abundance=get_protein_quant_table,
    output:
        "results/collated/co-assembly/{assembler}/{assembly}/taxonomy/taxonomy.{quant_type}.tsv"
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
        abundance = pd.read_csv(input.abundance, header=0, sep="\t", index_col=0)
        gene_tax_abundance = pd.merge(
            gene_tax, abundance, left_index=True, right_index=True
        )
        species_abundance = gene_tax_abundance.groupby(ranks).sum().reset_index()
        # Write results
        species_abundance.to_csv(output[0], sep="\t", index=False)
