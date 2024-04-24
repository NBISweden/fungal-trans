localrules:
    extract_fungi_reads,
    get_all_mapped_fungal_refs,
    alignment_report,
    count_reads,
    link_unfiltered

################
## Unfiltered ##
################
rule link_unfiltered:
    input:
        "results/preprocess/{sample_id}_R{i}.cut.trim.fastq.gz"
    output:
        "results/unfiltered/{sample_id}/{sample_id}_R{i}.cut.trim.fastq.gz"
    shell:
        """
        ln -s $(pwd)/{input} $(pwd)/{output}
        """

##########################
## Bowtie2/STAR mapping ##
##########################
include: "paired_strategy.smk"

rule bowtie_build_fungi:
    input:
        rules.concat_transcripts.output
    output:
        expand("resources/fungi/fungi_transcripts.fasta.{index}.bt2l",
               index=range(1,5))
    log:
        "resources/fungi/bowtie2.log"
    params:
        fasta = "$TMPDIR/fungi_transcripts/fungi_transcripts.fasta",
        tmpdir = "$TMPDIR/fungi_transcripts",
        outdir = "resources/fungi/"
    threads: 20
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*120
    conda: "../../envs/bowtie2.yaml"
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input} > {params.fasta}
        bowtie2-build \
            --threads {threads} \
            --large-index {input} \
            {params.fasta} >{log} 2>&1
        mv {params.fasta}*.bt2l {params.outdir}/
        rm -r {params.tmpdir}
        """

if not config["host_fna_url"] and not config["host_fna"] and config["host_aligner"]=="star":
    sys.exit("ERROR: No host genome fasta file given for STAR mapping")

def star_build_input(wildcards):
    input = []
    # if host_fna_url is give, put this download under resources/host/host.fna
    if config["host_fna"]:
        input.append(config["host_fna"])
    elif config["host_fna_url"]:
        input.append("resources/host/host.fna")
        config["host_fna"] = "resources/host/host.fna"
    if config["host_gff"]:
        input.append(config["host_gff"])
    elif config["host_gff_url"]:
        input.append("resources/host/host.gtf")
        config["host_gff"] = "resources/host/host.gff"
    return input

rule star_build_host:
    input:
        star_build_input,
    output:
        expand("resources/host/{f}",
            f = ["Genome", "SA", "SAindex", "chrLength.txt", "chrName.txt",
                 "chrNameLength.txt", "chrStart.txt", "genomeParameters.txt"])
    log:
        "resources/host/Log.out"
    params:
        overhang = "--sjdbOverhang "+str(config["star_overhang"]) if config["host_gff"] != "" else "",
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        gtfstring = lambda wildcards, input: f"--sjdbGTFfile {config['host_gff']}" if config['host_gff']!="" else ""
    threads: 20
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 10
    conda: "../../envs/star.yaml"
    shell:
        """
        exec &>{log}
         STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.outdir} \
            --genomeFastaFiles {input[0]} {params.overhang} {params.gtfstring}
        """


rule bowtie_build_host:
    input:
        "resources/host/host.fna"
    output:
        expand("resources/host/host.fna.{index}.bt2l", index=range(1,5))
    log:
        "resources/host/bowtie2.log"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    conda: "../../envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            --large-index {input} \
            {input} >{log} 2>&1
        """

rule bowtie_map_fungi:
    """
    Maps preprocessed reads against fungal transcripts.
    
    Reads that do or do not map concordantly in this step are saved directly
    in separate fastq files. If the 'paired_strategy' is set to 'concordant' 
    these files will be used to separate host and fungal reads.
    """
    input:
        R1="results/preprocess/{sample_id}_R1.cut.trim.fastq.gz",
        R2="results/preprocess/{sample_id}_R2.cut.trim.fastq.gz",
        db=expand("resources/fungi/fungi_transcripts.fasta.{index}.bt2l", index=range(1,5))
    output:
        bam="results/bowtie2/{sample_id}/{sample_id}.fungi.bam",
        R1="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.conc.fastq.gz",
        R2="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.conc.fastq.gz",
        R1u="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.noconc.fastq.gz",
        R2u="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.noconc.fastq.gz",
    params:
        prefix = "resources/fungi/fungi_transcripts.fasta",
        al_conc_path = "$TMPDIR/{sample_id}/{sample_id}_R%.fungi.conc.fastq.gz",
        un_conc_path = "$TMPDIR/{sample_id}/{sample_id}_R%.fungi.noconc.fastq.gz",
        R1 = "$TMPDIR/{sample_id}/{sample_id}_R1.fungi.conc.fastq.gz",
        R2 = "$TMPDIR/{sample_id}/{sample_id}_R2.fungi.conc.fastq.gz",
        R1u = "$TMPDIR/{sample_id}/{sample_id}_R1.fungi.noconc.fastq.gz",
        R2u = "$TMPDIR/{sample_id}/{sample_id}_R2.fungi.noconc.fastq.gz",
        tmpdir = "$TMPDIR/{sample_id}",
        temp_bam = "$TMPDIR/{sample_id}/{sample_id}.fungi.bam",
        setting = config["bowtie2_params"]
    log:
        bt2 = "results/bowtie2/{sample_id}/{sample_id}.bowtie2.fungi.log",
        st_sort = "results/bowtie2/{sample_id}/{sample_id}.samtools_sort.fungi.log",
    threads: 10
    resources:
        runtime = lambda wildcards, attempt: attempt**2*240
    shell:
        """
        mkdir -p {params.tmpdir}
        bowtie2 {params.setting} -p {threads} -x {params.prefix} -1 {input.R1} \
            -2 {input.R2} --al-conc-gz {params.al_conc_path} \
            --un-conc-gz {params.un_conc_path} 2> {log.bt2} | samtools sort -n -O BAM - >{params.temp_bam} 2>{log.st_sort} 
        mv {params.temp_bam} {output.bam}
        mv {params.R1} {output.R1} 
        mv {params.R2} {output.R2}
        mv {params.R1u} {output.R1u}
        mv {params.R2u} {output.R2u}
        """

rule get_mapped_fungal_refs:
    input:
        "results/bowtie2/{sample_id}/{sample_id}.fungi.bam"
    output:
        "results/bowtie2/{sample_id}/{sample_id}.fungi.refs"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*30
    shell:
        """
        samtools view -F 4 -f 64 {input[0]} | cut -f3 | sort -u > {output[0]}
        """

rule get_all_mapped_fungal_refs:
    input:
        expand("results/bowtie2/{sample_id}/{sample_id}.fungi.refs",
               sample_id = samples.keys())
    output:
        "results/collated/bowtie2/fungi.refs.ids"
    shell:
        """
        cat {input} | sort -u > {output[0]}
        """

rule star_map_host:
    """
    Maps putative fungal reads against host with STAR
    """
    input:
        db=expand("resources/host/{f}",
            f=["Genome", "SA", "SAindex", "chrLength.txt", "chrName.txt",
               "chrNameLength.txt", "chrStart.txt", "genomeParameters.txt"]),
        R1="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.fastq.gz",
        R2="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.fastq.gz"
    output:
        "results/star/{sample_id}/{sample_id}.host.bam"
    log:
        star="results/star/{sample_id}/{sample_id}.host.log",
        all="results/star/{sample_id}/map.log",
        stat="results/star/{sample_id}/{sample_id}.Log.final.out"
    params:
        prefix = "$TMPDIR/{sample_id}/{sample_id}.",
        genomedir = lambda wildcards, input: os.path.dirname(input.db[0]),
        temp_bam = "$TMPDIR/{sample_id}/{sample_id}.Aligned.out.bam",
        temp_log = "$TMPDIR/{sample_id}/{sample_id}.Log.out",
        temp_stat = "$TMPDIR/{sample_id}/{sample_id}.Log.final.out",
        setting = config["star_params"]
    conda: "../../envs/star.yaml"
    threads: 20
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 10
    shell:
        """
        exec &>{log.all}
        STAR --outFileNamePrefix {params.prefix} --runThreadN {threads} \
            --genomeDir {params.genomedir} --readFilesIn {input.R1} {input.R2} \
            --readFilesCommand 'gunzip -c' --outSAMtype BAM Unsorted \
            --outSAMunmapped Within KeepPairs {params.setting}
        mv {params.temp_bam} {output[0]}
        mv {params.temp_log} {log.star}
        mv {params.temp_stat} {log.stat}
        """

rule bowtie_map_host:
    """
    Maps putative fungal reads against host sequences with bowtie2 
    """
    input:
        db = expand("resources/host/host.fna.{index}.bt2l", index=range(1,5)),
        R1="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.fastq.gz",
        R2="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.fastq.gz"
    output:
        bam="results/bowtie2/{sample_id}/{sample_id}.host.bam",
        R1f="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.host.noconc.fastq.gz",
        R2f="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.host.noconc.fastq.gz",
        R1h="results/bowtie2/{sample_id}/{sample_id}_R1.fungi.host.conc.fastq.gz",
        R2h="results/bowtie2/{sample_id}/{sample_id}_R2.fungi.host.conc.fastq.gz"
    params:
        temp_bam = "$TMPDIR/{sample_id}/{sample_id}.host.bam",
        no_al_path = "$TMPDIR/{sample_id}/{sample_id}_R%.fungi.host.noconc.fastq.gz",
        al_path = "$TMPDIR/{sample_id}/{sample_id}_R%.fungi.host.conc.fastq.gz",
        R1h = "$TMPDIR/{sample_id}/{sample_id}_R1.fungi.host.conc.fastq.gz",
        R2h ="$TMPDIR/{sample_id}/{sample_id}_R2.fungi.host.conc.fastq.gz",
        R1f = "$TMPDIR/{sample_id}/{sample_id}_R1.fungi.host.noconc.fastq.gz",
        R2f = "$TMPDIR/{sample_id}/{sample_id}_R2.fungi.host.noconc.fastq.gz",
        tmpdir = "$TMPDIR/{sample_id}",
        prefix = "resources/host/host.fna",
        setting = config["bowtie2_params"]
    log:
        bt2 = "results/bowtie2/{sample_id}/{sample_id}.host.log",
        samtools = "results/bowtie2/{sample_id}/{sample_id}.samtools.host.log"
    threads: 10
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    shell:
        """
        exec &> {log.samtools}
        mkdir -p {params.tmpdir}
        bowtie2 {params.setting} -p {threads} -x {params.prefix} -1 {input.R1} \
            -2 {input.R2} --al-conc-gz {params.al_path} --un-conc-gz {params.no_al_path} \
            2>{log.bt2} | samtools sort -n -O BAM - > {params.temp_bam} 2>{log.samtools}
        mv {params.temp_bam} {output.bam}
        mv {params.R1f} {output.R1f}
        mv {params.R2f} {output.R2f}
        mv {params.R1h} {output.R1h}
        mv {params.R2h} {output.R2h}
        """

rule host_reads:
    input:
        R1="results/preprocess/{sample_id}_R1.cut.trim.fastq.gz",
        R2="results/preprocess/{sample_id}_R2.cut.trim.fastq.gz",
        R1_1 = "results/bowtie2/{sample_id}/{sample_id}_R1.nonfungi.fastq.gz",
        R2_1 = "results/bowtie2/{sample_id}/{sample_id}_R2.nonfungi.fastq.gz",
        R1_2 = "results/"+config["host_aligner"]+"/{sample_id}/{sample_id}_R1.fungi.putative-host.fastq.gz",
        R2_2 = "results/"+config["host_aligner"]+"/{sample_id}/{sample_id}_R2.fungi.putative-host.fastq.gz"
    output:
        R1 = "results/host/{sample_id}_R1.host.fastq.gz",
        R2 = "results/host/{sample_id}_R2.host.fastq.gz"
    params:
        tmpids="$TMPDIR/{sample_id}.tmp",
        ids = "$TMPDIR/{sample_id}.ids",
        R1 = "$TMPDIR/{sample_id}_R1.host.fastq.gz",
        R2 = "$TMPDIR/{sample_id}_R2.host.fastq.gz"
    log:
        "results/host/{sample_id}.log"
    shell:
        """
        set +o pipefail;
        exec &>{log}
        touch {params.tmpids}
        for f in {input.R1_1} {input.R1_2};
        do
            if [ -s $f ]; then
                gunzip -c $f | egrep "^@" | cut -f1 -d ' ' | sed 's/@//g' >> {params.tmpids}
            fi
        done
        for f in {input.R2_1} {input.R2_2};
        do
            if [ -s $f ]; then
                gunzip -c $f | egrep "^@" | cut -f1 -d ' ' | sed 's/@//g' >> {params.tmpids}
            fi
        done
        cat {params.tmpids} | sort -u > {params.ids}
        seqtk subseq {input.R1} {params.ids} | gzip -c > {params.R1}
        seqtk subseq {input.R2} {params.ids} | gzip -c > {params.R2}
        mv {params.R1} {output.R1}
        mv {params.R2} {output.R2}
        rm {params.ids}
        rm {params.tmpids}
        """

def get_host_logs(config, samples):
    if config["host_aligner"] == "star":
        return expand("results/star/{sample_id}/{sample_id}.host.log", sample_id = samples.keys())
    elif config["host_aligner"] == "bowtie2":
        return expand("results/bowtie2/{sample_id}/{sample_id}.host.log", sample_id = samples.keys())

rule alignment_report:
    input:
        hostlogs = get_host_logs(config, samples),
        fungallogs = expand("results/bowtie2/{sample_id}/{sample_id}.bowtie2.fungi.log",
            sample_id = samples.keys())
    output:
        "results/report/filtering/filter_report.html"
    params:
        tmpdir = "multiqc_filter",
        config = "config/multiqc_filter_config.yaml"
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.hostlogs} {params.tmpdir}
        cp {input.fungallogs} {params.tmpdir}
        multiqc \
            -f -c {params.config} \
            -n filter_report \
            -o results/report/filtering {params.tmpdir}
        rm -r {params.tmpdir}
        """

########################
## Taxmapper searches ##
########################

rule taxmapper_search:
    input:
        R1="results/preprocess/{sample_id}_R1.cut.trim.fastq.gz",
        R2="results/preprocess/{sample_id}_R2.cut.trim.fastq.gz",
        db="resources/taxmapper/databases/taxonomy/meta_database.db"
    output:
        temp(expand("results/taxmapper/{{sample_id}}/hits_{i}.aln", i = [1,2]))
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*6
    params:
        tmp_dir = "$TMPDIR/{sample_id}",
        out_dir = "results/taxmapper/{sample_id}"
    threads: 20
    conda: "../../envs/taxmapper.yaml"
    shell:
        """
        mkdir -p {params.tmp_dir}
        gunzip -c {input.R1} > {params.tmp_dir}/R1.fastq
        gunzip -c {input.R2} > {params.tmp_dir}/R2.fastq
        taxmapper search \
            -f {params.tmp_dir}/R1.fastq \
            -r {params.tmp_dir}/R2.fastq \
            -d {input.db} -t {threads} -o {params.tmp_dir}/hits
         mv {params.tmp_dir}/*.aln {params.out_dir}
         rm -r {params.tmp_dir}
        """

rule taxmapper_map:
    input:
        R1 = "results/taxmapper/{sample_id}/hits_1.aln",
        R2 = "results/taxmapper/{sample_id}/hits_2.aln",
    output:
        "results/taxmapper/{sample_id}/taxa.tsv.gz",
        "results/taxmapper/{sample_id}/taxa_identities.tsv"
    resources:
        runtime = lambda wildcards, attempt: attempt*60*3
    conda: "../../envs/taxmapper.yaml"
    threads: 4
    params:
        out_dir = "results/taxmapper/{sample_id}",
        unzipped = "results/taxmapper/{sample_id}/taxa.tsv"
    shell:
        """
        taxmapper map \
            -f {input.R1} -r {input.R2} \
            -o {params.out_dir}/taxa.tsv -t {threads} -m 100
        gzip {params.unzipped}
        """

rule taxmapper_filter_reads:
    input:
        "results/taxmapper/{sample_id}/taxa.tsv.gz"
    output:
        "results/taxmapper/{sample_id}/taxa_filtered.tsv.gz"
    resources:
        runtime = lambda wildcards, attempt: attempt*10
    conda: "../../envs/taxmapper.yaml"
    threads: 1
    params:
        tmp_in = "$TMPDIR/taxa.tsv",
        tmp_out = "$TMPDIR/taxa_filtered.tsv"
    shell:
        """
        gunzip -c {input[0]} > {params.tmp_in}
        taxmapper filter -i {params.tmp_in} -o {params.tmp_out} -a 0.25
        gzip {params.tmp_out}
        mv {params.tmp_out}.gz {output[0]}
        """

rule taxmapper_count:
    input:
        "results/taxmapper/{sample_id}/taxa_filtered.tsv.gz"
    output:
        lvl1 = "results/taxmapper/{sample_id}/taxa_counts_level1.tsv",
        lvl2 = "results/taxmapper/{sample_id}/taxa_counts_level2.tsv"
    params:
        tmp_in = "$TMPDIR/taxa_filtered.tsv"
    conda: "../../envs/taxmapper.yaml"
    resources:
        runtime = lambda wildcards, attempt: attempt*attempt*10
    shell:
        """
        gunzip -c {input[0]} > {params.tmp_in}
        taxmapper count \
            -i {params.tmp_in} --out1 {output.lvl1} --out2 {output.lvl2}
        """

rule extract_fungi_reads:
    input:
        tsv = "results/taxmapper/{sample_id}/taxa_filtered.tsv.gz",
        R1="results/preprocess/{sample_id}_R1.cut.trim.fastq.gz",
        R2="results/preprocess/{sample_id}_R2.cut.trim.fastq.gz"
    output:
        R1 = "results/taxmapper/{sample_id}/{sample_id}_R1.cut.trim.filtered.fastq.gz",
        R2 = "results/taxmapper/{sample_id}/{sample_id}_R2.cut.trim.filtered.fastq.gz"
    params:
        R1_ids = "results/taxmapper/{sample_id}/fungi.ids1",
        R2_ids = "results/taxmapper/{sample_id}/fungi.ids2"
    shell:
        """
        gunzip -c {input.tsv} | grep -w "Fungi" | cut -f1 | sed "s/\/2/\/1/g" > {params.R1_ids}
        gunzip -c {input.tsv} | grep -w "Fungi" | cut -f1 | sed "s/\/1/\/2/g" > {params.R2_ids}
        seqtk subseq {input.R1} {params.R1_ids} | gzip -c > {output.R1}
        seqtk subseq {input.R2} {params.R2_ids} | gzip -c > {output.R2}
        rm {params.R1_ids} {params.R2_ids}
        """

####################################
## Count reads at different steps ##
####################################
def get_ids(f):
    l = []
    for line in shell("seqtk comp {f} | cut -f1", iterable = True):
        l.append(line.rstrip())
    return set(l)

rule count_reads:
    input:
        R1tm = expand("results/taxmapper/{sample_id}/{sample_id}_R1.cut.trim.filtered.fastq.gz",
            sample_id = samples.keys()) if config["read_source"] in ["taxmapper", "filtered"] else [],
        R1bt = expand("results/bowtie2/{sample_id}/{sample_id}_R1.fungi.fastq.gz",
            sample_id = samples.keys()),
        R1btf = expand("results/{aligner}/{sample_id}/{sample_id}_R1.fungi.nohost.fastq.gz",
            sample_id = samples.keys(), aligner = config["host_aligner"]),
        R1bts = expand("results/{aligner}/{sample_id}/{sample_id}_R1.fungi.putative-host.fastq.gz",
            sample_id = samples.keys(), aligner = config["host_aligner"])
    output:
        "results/report/filtering/filtered_read_counts.tsv"
    run:
        counts = {}
        for i,sample in enumerate(samples.keys(), start=1):
            print("{} ({}/{})".format(sample, i, len(samples.keys())))
            counts[sample] = {"taxmapper_tot": 0,"bowtie_tot": 0,"bowtie_nonhost_tot": 0,"bowtie_host_tot": 0,
                "taxmapper_bowtie": 0, "taxmapper_bowtie_nonhost": 0, "taxmapper_bowtie_host": 0,
                "union": 0}
            # Taxmapper read ids
            tmfile = "results/taxmapper/{sample_id}/{sample_id}_R1.cut.trim.filtered.fastq.gz".format(sample_id=sample)
            tmids = get_ids(tmfile)
            # Bowtie read ids
            btfile = "results/bowtie2/{sample_id}/{sample_id}_R1.fungi.fastq.gz".format(sample_id=sample)
            btids = get_ids(btfile)
            # Bowtie read ids non host
            btffile = "results/{aligner}/{sample_id}/{sample_id}_R1.fungi.nohost.fastq.gz".format(sample_id=sample)
            btfids = get_ids(btffile)
            # Bowtie read ids host
            btsfile = "results/bowtie2/{sample_id}/{sample_id}_R1.fungi.host.fastq.gz".format(sample_id=sample)
            btsids = get_ids(btsfile)
            # Count reads common and unique
            counts[sample]["taxmapper_tot"] = len(tmids)
            counts[sample]["bowtie_tot"] = len(btids)
            counts[sample]["bowtie_nonhost_tot"] = len(btfids)
            counts[sample]["bowtie_host_tot"] = len(btsids)
            counts[sample]["taxmapper_bowtie"] = len(tmids.intersection(btids))
            counts[sample]["taxmapper_bowtie_nonhost"] = len(tmids.intersection(btfids))
            counts[sample]["taxmapper_bowtie_host"] = len(tmids.intersection(btsids))
            counts[sample]["union"] = len(tmids.union(btfids))
        df = pd.DataFrame(counts).T
        df.index.name="Sample"
        df = df[["taxmapper_tot","bowtie_nonhost_tot","union","bowtie_tot","bowtie_host_tot","taxmapper_bowtie","taxmapper_bowtie_nonhost","taxmapper_bowtie_host"]]
        df.to_csv(output[0], sep="\t")

###########################
## Make union reads file ##
###########################
def make_idfile(in1, in2, out, tmpdir):
    shell("seqtk comp {in1} | cut -f1 > {tmpdir}/tmp")
    shell("seqtk comp {in2} | cut -f1 >> {tmpdir}/tmp")
    shell("sort -u {tmpdir}/tmp > {out}")
    shell("rm {tmpdir}/tmp")

def find_reads(infile, fhout, wanted):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    with open(os.path.expandvars(infile)) as fhin:
        for title, seq, qual in FastqGeneralIterator(fhin):
            seqid = title.split(None, 1)[0]
            if seqid in wanted:
                fhout.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))
                wanted.remove(seqid)
    return wanted

def extract_union_reads(in1, in2, out, idfile):
    with open(os.path.expandvars(idfile)) as id_handle:
        # Taking first word on each line as an identifer
        wanted = set(line.rstrip("\n").split(None,1)[0] for line in id_handle)
    with open(os.path.expandvars(out), 'w') as fhout:
        wanted = find_reads(in1, fhout, wanted)
        wanted = find_reads(in2, fhout, wanted)
    return len(wanted)

rule union_filtered_reads:
    input:
        R1_taxmapper = "results/taxmapper/{sample_id}/{sample_id}_R1.cut.trim.filtered.fastq.gz",
        R2_taxmapper = "results/taxmapper/{sample_id}/{sample_id}_R2.cut.trim.filtered.fastq.gz",
        R1_bowtie = "results/bowtie2/{sample_id}/{sample_id}_R1.fungi.nohost.fastq.gz",
        R2_bowtie = "results/bowtie2/{sample_id}/{sample_id}_R2.fungi.nohost.fastq.gz"
    output:
        R1 = "results/filtered/{sample_id}/{sample_id}_R1.filtered.union.fastq.gz",
        R2 = "results/filtered/{sample_id}/{sample_id}_R2.filtered.union.fastq.gz"
    params:
        R1_ids = "$TMPDIR/{sample_id}/R1.ids",
        R2_ids = "$TMPDIR/{sample_id}/R2.ids",
        R1_fastq = "$TMPDIR/{sample_id}/R1.fastq",
        R2_fastq = "$TMPDIR/{sample_id}/R2.fastq",
        tm_tmp_fastq = "$TMPDIR/{sample_id}/tm.fastq",
        bt_tmp_fastq = "$TMPDIR/{sample_id}/bt.fastq",
        tmpdir = "$TMPDIR/{sample_id}"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    run:
        shell("mkdir -p {params.tmpdir}")
        print("Storing read ids")
        make_idfile(input.R1_taxmapper, input.R1_bowtie, params.R1_ids, params.tmpdir)
        make_idfile(input.R2_taxmapper, input.R2_bowtie, params.R2_ids, params.tmpdir)
        # Handle R1
        shell("gunzip -c {input.R1_taxmapper} > {params.tm_tmp_fastq}")
        shell("gunzip -c {input.R1_bowtie} > {params.bt_tmp_fastq}")
        print("Extracting reads from {} and {}".format(input.R1_taxmapper, input.R1_bowtie))
        remaining = extract_union_reads(params.tm_tmp_fastq, params.bt_tmp_fastq, params.R1_fastq, params.R1_ids)
        if remaining != 0:
            sys.exit("WARNING: Could not get all reads for {output.R1}\n")
        else:
            print("All reads extracted to {}".format(params.R1_fastq))
        # Handle R2
        shell("gunzip -c {input.R2_taxmapper} > {params.tm_tmp_fastq}")
        shell("gunzip -c {input.R2_bowtie} > {params.bt_tmp_fastq}")
        print("Extracting reads from {} and {}".format(input.R2_taxmapper, input.R2_bowtie))
        remaining = extract_union_reads(params.tm_tmp_fastq, params.bt_tmp_fastq, params.R2_fastq, params.R2_ids)
        if remaining != 0:
            sys.exit("WARNING: Could not get all reads for {output.R2}\n")
        else:
            print("All reads extracted to {}".format(params.R2_fastq))
        # Zip to outfiles
        print("Compressing output to {}".format(output.R1))
        shell("gzip -c {params.R1_fastq} > {params.R1_fastq}.gz; mv {params.R1_fastq}.gz {output.R1}")
        print("Compressing output to {}".format(output.R2))
        shell("gzip -c {params.R2_fastq} > {params.R2_fastq}.gz; mv {params.R2_fastq}.gz {output.R2}")
        # Cleanup
        shell("rm -rf {params.tmpdir}/*")
