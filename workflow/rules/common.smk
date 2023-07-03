import pandas as pd
import os, sys
import collections
import yaml
import subprocess

if config["mosaicatcher_pipeline"] == False:

    exclude = [
        "._.DS_Store",
        ".DS_Store",
        "all",
        "ashleys_counts",
        "bam",
        "cell_selection",
        "config",
        "counts",
        "fastq",
        "fastqc",
        "haplotag",
        "log",
        "merged_bam",
        "mosaiclassifier",
        "normalizations",
        "ploidy",
        "plots",
        "predictions",
        "segmentation",
        "snv_calls",
        "stats",
        "strandphaser",
    ]

    if config["chromosomes_to_exclude"]:
        chroms_init = config["chromosomes"]
        chroms = [e for e in chroms_init if e not in config["chromosomes_to_exclude"]]
        config["chromosomes"] = chroms

    if config["reference"] == "mm10":
        config["chromosomes"] = [
            "chr" + str(e) for e in list(range(1, 20)) + ["X", "Y"]
        ]

    from scripts.utils import make_log_useful_ashleys, pipeline_aesthetic_start_ashleys

    if config["list_commands"] is True:
        pipeline_aesthetic_start_ashleys.argparse_help(config)

    onstart:
        pipeline_aesthetic_start_ashleys.pipeline_aesthetic_start(config)
        subprocess.Popen(
            "mkdir -p {folder_path}/config".format(
                    folder_path=config["data_location"]
                ),
                shell=True,
                stdout=subprocess.PIPE,
            )
        subprocess.Popen(
            "rsync --ignore-existing -avzh config/config.yaml {folder_path}/config".format(
                    folder_path=config["data_location"]
                ),
                shell=True,
                stdout=subprocess.PIPE,
            )

            # for sample in [e for e in os.listdir(config["data_location"]) if e not in exclude]:
            #     print(sample)
            #     if len(sample.split("_")) == 4:
            #         assert len(sample.split("_")) != 4, "Your sample name is using 4 times the '_' character, which is currently not supported by ashleys-qc"


    def onsuccess_fct(log):

        make_log_useful_ashleys.make_log_useful(log, "SUCCESS", config)
        shell(
            'mail -s "[Snakemake] smk-wf-catalog/ashleys-qc-pipeline v{} - Run on {} - SUCCESS" {} < {{log}}'.format(
                config["version"], config["data_location"], config["email"]
            )
        )

    def onerror_fct(log):
        make_log_useful_ashleys.make_log_useful(log, "ERROR", config)
        shell(
            'mail -s "[Snakemake] smk-wf-catalog/ashleys-qc-pipeline v{} - Run on {} - ERRROR" {} < {{log}}'.format(
                config["version"], config["data_location"], config["email"]
            )
        )


# Simple class to retrieve automatically files in the fastq/bam folder and create a config dataframe
class HandleInput:
    def __init__(
        self,
        input_path,
        output_path,
        check_sm_tag=False,
        bam=True,
        genecore=False,
        genecore_path=str,
    ):
        if genecore is False:
            df_config_files = self.handle_input_data(thisdir=input_path, bam=bam)
        elif genecore is True:
            df_config_files, d_master = self.handle_input_data_genecore(
                thisdir=genecore_path
            )
            self.d_master = d_master

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        df_config_files.to_csv(output_path, sep="\t", index=False)
        self.df_config_files = df_config_files

    @staticmethod
    def handle_input_data_genecore(thisdir):
        """_summary_

        Args:
            thisdir (_type_): _description_
            exclude_list (_type_, optional): _description_. Defaults to list.

        Returns:
            _type_: _description_
        """
        complete_df_list = list()

        # List of folders/files to not consider (restrict to samples only)
        l = sorted(
            [
                e
                for e in os.listdir(
                    "{genecore_prefix}/{date_folder}".format(
                        genecore_prefix=config["genecore_prefix"],
                        date_folder=config["genecore_date_folder"],
                    )
                )
                if e.endswith(".txt.gz")
            ]
        )

        # Create a list of  files to process for each sample
        d_master = collections.defaultdict(dict)
        sub_l = list()
        for j, e in enumerate(l):
            sub_l.append(e)
            if (j + 1) % 192 == 0:
                common_element = findstem(sub_l)
                l_elems = common_element.split("lane1")
                # print(sub_l)
                # print(common_element)
                # print(l_elems)
                # print(l_elems[1].split("{regex_element}".format(regex_element=config["genecore_regex_element"]))
                prefix = l_elems[0]
                # technician_name = l_elems[0].split("_")[-2]
                sample = l_elems[1].split(
                    "{regex_element}".format(
                        regex_element=config["genecore_regex_element"]
                    )
                )[0]
                index = l_elems[1].split(
                    "{regex_element}".format(
                        regex_element=config["genecore_regex_element"]
                    )
                )[1]
                # pe_index = common_element[-1]
                sub_l = list()

                d_master[sample]["prefix"] = prefix
                # d_master[sample]["technician_name"] = technician_name
                d_master[sample]["index"] = index
                d_master[sample]["common_element"] = common_element
        # from pprint import pprint
        # pprint(d_master)
        # exit()
        samples_to_process = (
            config["samples_to_process"]
            if len(config["samples_to_process"]) > 0
            else list(d_master.keys())
        )

        config["data_location"] = "{data_location}/{genecore_date_folder}".format(
            data_location=config["data_location"],
            genecore_date_folder=config["genecore_date_folder"],
        )

        genecore_list = [
            expand(
                "{data_location}/{sample}/fastq/{sample}{regex_element}{index}{cell_nb}.{pair}.fastq.gz",
                data_location=config["data_location"],
                sample=sample,
                regex_element=config["genecore_regex_element"],
                index=d_master[sample]["index"],
                cell_nb=[str(e).zfill(2) for e in list(range(1, 97))],
                pair=["1", "2"],
            )
            for sample in d_master
            if sample in samples_to_process
        ]
        genecore_list = [sub_e for e in genecore_list for sub_e in e]
        # pprint(genecore_list)
        complete_df_list = list()

        for sample in d_master:
            df = pd.DataFrame(
                [
                    {"File": os.path.basename(f), "Folder": os.path.dirname(f)}
                    for f in genecore_list
                    if sample in f
                ]
            )
            if df.shape[0] > 0:
                df["File"] = df["File"].str.replace(".fastq.gz", "", regex=True)
                df["Sample"] = sample
                df["Pair"] = df["File"].apply(lambda r: r.split(".")[1])
                df["Cell"] = df["File"].apply(lambda r: r.split(".")[0])
                df["Full_path"] = df[["Folder", "File"]].apply(
                    lambda r: f"{r['Folder']}/{r['File']}.fastq.gz", axis=1
                )
                df["Genecore_path"] = df["File"].apply(
                    lambda r: f"{config['genecore_prefix']}/{config['genecore_date_folder']}/{d_master[sample]['prefix']}lane1{r.replace('.', '_')}_sequence.txt.gz"
                )
                df["Genecore_file"] = df["File"].apply(
                    lambda r: f"{d_master[sample]['prefix']}lane1{r.replace('.', '_')}"
                )
                df["Genecore_file"] = df["Genecore_file"].apply(
                    lambda r: "_".join(r.split("_")[:-1])
                )

                # Concat dataframes for each sample & output
                complete_df_list.append(df)

        complete_df = pd.concat(complete_df_list)

        complete_df = complete_df.sort_values(by=["Cell", "File"]).reset_index(
            drop=True
        )
        pd.options.display.max_colwidth = 200
        # print(complete_df)
        return complete_df, d_master

    @staticmethod
    def handle_input_data(thisdir, exclude_list=list, bam=bool):
        """_summary_

        Args:
            thisdir (_type_): _description_
            exclude_list (_type_, optional): _description_. Defaults to list.

        Returns:
            _type_: _description_
        """
        # Extension & folder based on bam boolean input
        ext = ".bam" if bam is True else ".fastq.gz"
        folder = "bam" if bam is True else "fastq"
        complete_df_list = list()
        # List of folders/files to not consider (restrict to samples only)
        exclude = [
            "._.DS_Store",
            ".DS_Store",
            "all",
            "ashleys_counts",
            "bam",
            "cell_selection",
            "config",
            "counts",
            "fastq",
            "fastqc",
            "haplotag",
            "log",
            "merged_bam",
            "mosaiclassifier",
            "normalizations",
            "ploidy",
            "plots",
            "predictions",
            "segmentation",
            "snv_calls",
            "stats",
            "strandphaser",
        ]

        for sample in [e for e in os.listdir(thisdir) if e not in exclude]:
            # Create a list of  files to process for each sample
            l_files_all = [
                f
                for f in os.listdir(
                    "{thisdir}/{sample}/{folder}/".format(
                        thisdir=thisdir, sample=sample, folder=folder
                    )
                )
                if f.endswith(ext)
            ]

            # print(l_files_all)
            # Dataframe creation
            df = pd.DataFrame([{"File": f} for f in l_files_all])
            df["File"] = df["File"].str.replace(ext, "", regex=True)
            df["Folder"] = thisdir
            df["Sample"] = sample
            df["Cell"] = df["File"].apply(lambda r: r.split(".")[0])
            df["Full_path"] = "{thisdir}/{sample}/{folder}/".format(
                thisdir=thisdir, sample=sample, folder=folder
            )
            df["Full_path"] = df["Full_path"] + df["File"] + ext

            complete_df_list.append(df)

        # Concat dataframes for each sample & output
        complete_df = pd.concat(complete_df_list)
        complete_df = complete_df.sort_values(by=["Cell", "File"]).reset_index(
            drop=True
        )
        return complete_df


# GENECORE


def findstem(arr):

    # Determine size of the array
    n = len(arr)

    # Take first word from array
    # as reference
    s = arr[0]
    l = len(s)

    res = ""

    for i in range(l):
        for j in range(i + 1, l + 1):

            # generating all possible substrings
            # of our reference string arr[0] i.e s
            stem = s[i:j]
            k = 1
            for k in range(1, n):

                # Check if the generated stem is
                # common to all words
                if stem not in arr[k]:
                    break

            # If current substring is present in
            # all strings and its length is greater
            # than current result
            if k + 1 == n and len(res) < len(stem):
                res = stem

    return res


# Create configuration file with samples

c = HandleInput(
    input_path=config["data_location"],
    genecore_path="{genecore_prefix}/{genecore_date_folder}".format(
        genecore_prefix=config["genecore_prefix"],
        genecore_date_folder=config["genecore_date_folder"],
    ),
    output_path="{data_location}/config/config_df_ashleys.tsv".format(
        data_location=config["data_location"]
    ),
    check_sm_tag=False,
    bam=False,
    genecore=config["genecore"],
)
df_config_files = c.df_config_files
if config["genecore"] is True:
    d_master = c.d_master


samples = list(sorted(list(df_config_files.Sample.unique().tolist())))

# genecore_mapping = df_config_files.groupby("Genecore_file")["Cell"].unique().apply(lambda r: r[0]).to_dict()

# Dictionnary of libraries available per sample
cell_per_sample = (
    df_config_files.groupby("Sample")["Cell"].unique().apply(list).to_dict()
)

# Plottype options for QC count plot
plottype_counts = (
    config["plottype_counts"]
    if config["multistep_normalisation"] is True
    else config["plottype_counts"][0]
)
# print(plottype_counts)


def get_final_output():
    """
    Function called by snakemake rule all to run the pipeline
    """
    final_list = list()

    # MultiQC outputs

    if config["MultiQC"] is True:
        final_list.extend(
            expand(
                "{path}/{sample}/multiqc/multiqc_report/multiqc_report.html",
                path=config["data_location"],
                sample=samples,
            ),
        )

    if config["mosaicatcher_pipeline"] is False:

        final_list.extend(
            expand(
                "{path}/{sample}/cell_selection/labels.tsv",
                path=config["data_location"],
                sample=samples,
            )
        )

        # QC count plots (classic only or classic + corrected based on config multistep_normalisation option)

        final_list.extend(
            expand(
                "{output_folder}/{sample}/plots/counts/CountComplete.{plottype_counts}.pdf",
                output_folder=config["data_location"],
                sample=samples,
                plottype_counts=plottype_counts,
            ),
        )

    # Plate plots

    for sample in samples:

        if len(cell_per_sample[sample]) == 96:

            final_list.extend(
                [
                    sub_e
                    for e in [
                        expand(
                            "{path}/{sample}/plots/plate/ashleys_plate_{plate_plot}.pdf",
                            path=config["data_location"],
                            sample=sample,
                            plate_plot=["predictions", "probabilities"],
                        )
                    ]
                    for sub_e in e
                ]
            )

    if config["publishdir"] != "":
        final_list.extend(
            expand(
                "{folder}/{sample}/config/publishdir_outputs.ok",
                folder=config["data_location"],
                sample=samples,
            )
        )

    # print(final_list)
    return final_list


def publishdir_fct():
    """
    Restricted for ASHLEYS at the moment
    Backup files on a secondary location
    """
    list_files_to_copy = [
        "{folder}/{sample}/cell_selection/labels_raw.tsv",
        "{folder}/{sample}/cell_selection/labels.tsv",
        "{folder}/{sample}/counts/{sample}.info_raw",
        "{folder}/{sample}/counts/{sample}.txt.raw.gz",
        "config/config.yaml",
    ]
    final_list = [
        expand(e, folder=config["data_location"], sample=samples)
        for e in list_files_to_copy
    ]
    final_list.extend(
        expand(
            "{folder}/{sample}/plots/counts/CountComplete.{plottype_counts}.pdf",
            folder=config["data_location"],
            sample=samples,
            plottype_counts=plottype_counts,
        )
    )

    if config["use_light_data"] is False:

        final_list.extend(
            expand(
                "{folder}/{sample}/plots/plate/ashleys_plate_{plate_plot}.pdf",
                folder=config["data_location"],
                sample=samples,
                plate_plot=["predictions", "probabilities"],
            )
        )
        final_list.extend(
            expand(
                "{folder}/{sample}/cell_selection/labels_positive_control_corrected.tsv",
                folder=config["data_location"],
                sample=samples,
            )
        )
        final_list.extend(
            expand(
                "{folder}/{sample}/config/bypass_cell.txt",
                folder=config["data_location"],
                sample=samples,
            )
        )

    return final_list


def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [2, 4, 8, 16, 32]
    return mem_avail[attempt - 1] * 1000


def get_mem_mb_heavy(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [16, 32, 64, 128, 256]
    return mem_avail[attempt - 1] * 1000
