rule fastqc:
    input:
        "{path}/{sample}/fastq/{cell}.{pair}.fastq.gz",
    output:
        html="{path}/{sample}/fastqc/{cell}_{pair}_fastqc.html",
        zip="{path}/{sample}/fastqc/{cell}_{pair}_fastqc.zip",
    params:
        "--quiet",
    log:
        "{path}/log/fastqc/{sample}/{cell}_{pair}.log",
    threads: 1
    resources:
        mem_mb=get_mem_mb,
        time="10:00:00",
    wrapper:
        "v1.7.0/bio/fastqc"


rule bwa_index:
    input:
        config["reference"],
    output:
        idx=multiext(config["reference"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "{}.log".format(config["reference"]),
    params:
        algorithm="bwtsw",
    resources:
        mem_mb=get_mem_mb,
        time="10:00:00",
    wrapper:
        "v1.7.0/bio/bwa/index"


rule bwa_strandseq_to_reference_alignment:
    input:
        mate1="{path}/{sample}/fastq/{cell}.1.fastq.gz",
        mate2="{path}/{sample}/fastq/{cell}.2.fastq.gz",
        ref_index=config["reference"],
    output:
        bam=temp("{path}/{sample}/all/{cell}.bam"),
    log:
        bwa="{path}/{sample}/log/{cell}.bwa.log",
        samtools="{path}/{sample}/log/{cell}.samtools.log",
    threads: 6
    params:
        idx_prefix=lambda wildcards, input: input.ref_index.rsplit(".", 1)[0],
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "bwa mem -t {threads}"
        ' -R "@RG\\tID:{wildcards.cell}\\tPL:Illumina\\tSM:{wildcards.sample}"'
        " -v 2 {input.ref_index} {input.mate1} {input.mate2} 2> {log.bwa} | "
        " samtools view -b /dev/stdin > {output.bam} 2> {log.samtools}"


rule samtools_sort_bam:
    input:
        "{path}/{sample}/all/{cell}.bam",
    output:
        temp("{path}/{sample}/all/{cell}.sort.bam"),
    log:
        "{path}/{sample}/log/samtools_sort/{cell}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "samtools sort -O BAM -o {output} {input} 2>&1 > {log}"

rule mark_duplicates:
    input:
        bam="{path}/{sample}/all/{cell}.sort.bam",
    output:
        "{path}/{sample}/all/{cell}.sort.mdup.bam",
    log:
        "{path}/{sample}/log/markdup/{cell}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
        time="10:00:00",
    shell:
        "sambamba markdup {input.bam} {output} 2>&1 > {log}"

if config["mosaicatcher_pipeline"] is False:
    rule samtools_index:
        input:
            "{path}/{sample}/all/{cell}.sort.mdup.bam",
        output:
            "{path}/{sample}/all/{cell}.sort.mdup.bam.bai",
        log:
            "{path}/{sample}/log/samtools_index/{cell}.log",
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        shell:
            "samtools index {input} 2>&1 > {log}"

### ASHLEYS AUTOMATED ANALYSIS
if config["hand_selection"] is False:

    rule install_ashleys:
        output:
            touch("{path}/config/ashleys_install_success.txt"),
        log:
            "{path}/log/config/install_ashleys.log",
        conda:
            "../envs/ashleys.yaml"
        shell:
            """
            FOLDER='ashleys-qc'
            if [ ! -d "$FOLDER" ] ; then
                git clone https://github.com/friendsofstrandseq/ashleys-qc.git
                cd "$FOLDER"
                python setup.py install 2>&1 > {log}
            fi
            """



    rule generate_features:
        input:
            # ashleys="{path}/config/ashleys_install_success.txt",
            bam=expand(
                "{path}/{sample}/all/{cell}.sort.mdup.bam",
                zip,
                path=input_bam_location_expand,
                sample=samples_expand,
                cell=cell_expand,
            ),
        output:
            "{path}/{sample}/predictions/ashleys_features.tsv",
        log:
            "{path}/log/ashleys/{sample}/features.log",
        conda:
            "../envs/ashleys.yaml" 
        params:
            windows="5000000 2000000 1000000 800000 600000 400000 200000",
            jobs=23,
            extension=".sort.mdup.bam",
            path=lambda wildcards, input: "{}all".format(input.bam[0].split("all")[0]),
            chroms=config["chromosomes"]
        shell:
            "ashleys -j {params.jobs} features -c {params.chroms} -f {params.path} -w {params.windows} -o {output} --recursive_collect -e {params.extension}"

    rule predict:
        input:
            path="{path}/{sample}/predictions/ashleys_features.tsv",
        output:
            "{path}/{sample}/predictions/predictions_raw.tsv",
        log:
            "{path}/log/ashleys/{sample}/prediction_ashleys.log",
        conda:
            "../envs/ashleys.yaml"
        params:
            model_default="workflow/models/svc_default.pkl",
            model_stringent="workflow/models/svc_stringent.pkl",
        resources:
            mem_mb=get_mem_mb,
            time="10:00:00",
        shell:
            "ashleys predict -p {input.path} -o {output} -m {params.model_default}"


### HAND SELECTION VIA JUPYTER NB
elif config["hand_selection"] is True:
    rule generate_exclude_file_for_mosaic_count:
        """
        rule fct: 
        input:
        output:
        """
        input:
            ancient("{path}/config/config_df_ashleys.tsv".format(path=config["input_bam_location"])),
            # ancient("config/samples.tsv"),
            bam=config["input_bam_location"],
        output:
            "{path}/config/exclude_file",
        log:
            "{path}/log/config/exclude_file.log",
        params:
            chroms=config["chromosomes"],
        # conda:
        #     "../envs/mc_base.yaml"
        script:
            "../scripts/utils/generate_exclude_file.py"
            
    rule mosaic_count:
        """
        rule fct: Call mosaic count C++ function to count reads in each BAM file according defined window
        input: For the moment, individual BAM file in the selected input_bam_location of the associated sample
        output: counts: read counts for the BAM file according defined window ; info file : summary statistics 
        """
        input:
            bam=lambda wc: expand(
                "{path}/{sample}/all/{cell}.sort.mdup.bam",
                path=config["input_bam_location"],
                sample=samples,
                cell=cell_per_sample[str(wc.sample)]
                if wc.sample in cell_per_sample
                else "FOOBAR",
            ),
            bai=lambda wc: expand(
                "{path}/{sample}/all/{cell}.sort.mdup.bam.bai",
                path=config["input_bam_location"],
                sample=samples,
                cell=cell_per_sample[str(wc.sample)],
            )
            if wc.sample in cell_per_sample
            else "FOOBAR",
            excl="{path}/config/exclude_file",
        output:
            counts="{path}/{sample}/ashleys_counts/{sample}.all.txt.fixme.gz",
            info="{path}/{sample}/ashleys_counts/{sample}.all.info",
        log:
            "{path}/log/counts/{sample}/mosaic_count.log",
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        params:
            window=config["window"],
        resources:
            mem_mb=get_mem_mb,
        shell:
            """
            mosaicatcher count \
                --verbose \
                --do-not-blacklist-hmm \
                -o {output.counts} \
                -i {output.info} \
                -x {input.excl} \
                -w {params.window} \
                {input.bam} \
            > {log} 2>&1
            """


    rule order_mosaic_count_output:
        input:
            "{path}/{sample}/ashleys_counts/{sample}.all.txt.fixme.gz",
        output:
            "{path}/{sample}/ashleys_counts/{sample}.all.txt.gz",
        log:
            "{path}/log/ashleys_counts/{sample}/{sample}.log",
        run:
            df = pd.read_csv(input[0], compression="gzip", sep="\t")
            df = df.sort_values(by=["sample", "cell", "chrom", "start"])
            df.to_csv(output[0], index=False, compression="gzip", sep="\t")


    rule plot_mosaic_counts:
        """
        rule fct: Plot function of read counts for each bam file
        input: mosaic count outputs (counts & info)
        output: Generate figure based on couting results
        """
        input:
            counts="{path}/{sample}/ashleys_counts/{sample}.all.txt.gz",
            info="{path}/{sample}/ashleys_counts/{sample}.all.info",
        output:
            "{path}/{sample}/plots/ashleys_counts/CountComplete_{sample}.pdf",
        log:
            "{path}/log/plot_mosaic_counts/{sample}.log",
        conda:
            "../envs/rtools.yaml"
        resources:
            mem_mb=get_mem_mb,
        shell:
            """
            LC_CTYPE=C Rscript workflow/scripts/plotting/qc.R {input.counts} {input.info} {output} > {log} 2>&1
            """
    
    # PDF must be in the jupyter directory
    # symlink not possible due to jupyter errors (too many symlink)
    rule cp_pdf_for_jupyter:
        input:  
            pdf = expand("{path}/{sample}/plots/ashleys_counts/CountComplete_{sample}.pdf", path=config["input_bam_location"], sample=samples),
        output:
            ".snakemake/scripts/CountComplete_{sample}.pdf"
        shell:
            "ln -s {input.pdf} {output}"

    rule notebook_hand_selection:
        input:
            # pdf_raw = "{path}/plots/{sample}/counts/CountComplete_{sample}.pdf",
            pdf_symlink = ".snakemake/scripts/CountComplete_{sample}.pdf",
        output:
            path = "{path}/{sample}/predictions/predictions_raw.tsv",
        log:
            "{path}/log/hand_selection/{sample}/prediction_probabilities.log"
        params:
            cell_per_sample = cell_per_sample
        conda:
            "../envs/notebook.yaml"
        notebook:
            "../notebooks/hand_selection.py.ipynb"


########################################################
#                      DEV PART
########################################################

if config["use_light_data"] is False:
    rule cp_predictions:
        input:
            path = "{path}/{sample}/predictions/predictions_raw.tsv",
        output:
            path = "{path}/{sample}/predictions/predictions.tsv",
        shell:
            "cp {input.path} {output.path}"

elif config["use_light_data"] is True:
    rule dev_all_cells_correct:
        input:
            path = "{path}/{sample}/predictions/predictions_raw.tsv",
        output:
            path = "{path}/{sample}/predictions/predictions.tsv",
        run:
            df = pd.read_csv(input.path, sep="\t")
            df["prediction"] = 1
            df["probability"] = 1
            df.to_csv(output.path, sep="\t", index=False)
        




rule symlink_bam_bai:
    input:
        bam = "{path}/{sample}/all/{cell}.sort.mdup.bam",
        bai = "{path}/{sample}/all/{cell}.sort.mdup.bam.bai",
    output: 
        bam = "{path}/{sample}/selected/{cell}.sort.mdup.bam",
        bai = "{path}/{sample}/selected/{cell}.sort.mdup.bam.bai",
    log:
        "{path}/log/symlink/{sample}/{cell}.log"
    conda:
        "../envs/ashleys.yaml"
    shell:
        """
        ln -s {input.bam} {output.bam} && touch -h {output.bam} > {log} 2>&1
        ln -s {input.bai} {output.bai} && touch -h {output.bai} > {log} 2>&1
        """

rule remove_unselected:
    input:
        predictions = "{path}/{sample}/predictions/predictions.tsv",
        bam = lambda wc: expand(
                "{path}/{sample}/selected/{cell}.sort.mdup.bam",
                path=config["input_bam_location"],
                sample=samples,
                cell=cell_per_sample[str(wc.sample)]
                if wc.sample in cell_per_sample
                else "FOOBAR",
            ),
        bai = lambda wc: expand(
                "{path}/{sample}/selected/{cell}.sort.mdup.bam.bai",
                path=config["input_bam_location"],
                sample=samples,
                cell=cell_per_sample[str(wc.sample)]
                if wc.sample in cell_per_sample
                else "FOOBAR",
            ),
    output:
        "{path}/config/{sample}_selected_cells.ok"
    log:
        "{path}/log/remove_unselected/{sample}.log"
    params:
        path = config["input_bam_location"]
    conda:
        "../envs/ashleys.yaml"
    script:
        "../scripts/utils/rm_unselected_cells.py"