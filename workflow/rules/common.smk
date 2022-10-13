import pandas as pd
import os, sys


class HandleInput:
    def __init__(self, input_path, output_path, check_sm_tag=False, bam=True):
        df_config_files = self.handle_input_data(thisdir=input_path, bam=bam)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        df_config_files.to_csv(output_path, sep="\t", index=False)
        self.df_config_files = df_config_files

    @staticmethod
    def handle_input_data(thisdir, exclude_list=list, bam=bool):
        """_summary_

        Args:
            thisdir (_type_): _description_
            exclude_list (_type_, optional): _description_. Defaults to list.

        Returns:
            _type_: _description_
        """
        ext = ".bam" if bam is True else ".fastq.gz"
        folder = "bam" if bam is True else "fastq"
        complete_df_list = list()
        exclude = ["._.DS_Store", ".DS_Store", "all", "ashleys_counts", "bam", "cell_selection", "config", "counts", "fastq", "fastqc", "haplotag", "log", "merged_bam", "mosaiclassifier", "normalizations", "ploidy", "plots", "predictions", "segmentation", "snv_calls", "stats", "strandphaser" ]
        for sample in [
            e
            for e in os.listdir(thisdir)
            if e not in exclude
        ]:
            l_files_all = [
                f
                for f in os.listdir(
                    "{thisdir}/{sample}/{folder}/".format(
                        thisdir=thisdir, sample=sample, folder=folder
                    )
                )
                if f.endswith(ext)
            ]
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

        complete_df = pd.concat(complete_df_list)
        complete_df = complete_df.sort_values(by=["Cell", "File"]).reset_index(
            drop=True
        )
        return complete_df


# Create configuration file with samples

c = HandleInput(
    input_path=config["data_location"],
    output_path="{data_location}/config/config_df_ashleys.tsv".format(
        data_location=config["data_location"]
    ),
    check_sm_tag=False,
    bam=False,
)
df_config_files = c.df_config_files


samples = list(sorted(list(df_config_files.Sample.unique().tolist())))

# Dictionnary of libraries available per sample
cell_per_sample = (
    df_config_files.groupby("Sample")["Cell"].unique().apply(list).to_dict()
)

plottype_counts = (
    config["plottype_counts"]
    if config["GC_analysis"] is True
    else config["plottype_counts"][0]
)

plottype_counts = ["classic"]

if config["GC_analysis"] is True:
 
    import string
    import collections
    import numpy as np

    d = collections.defaultdict(dict)
    orientation = (8,12) if config["plate_orientation"] == "landscape" else (12,8)
    for sample in samples:
        if len(cell_per_sample[sample]) == 96:
            for j, e in enumerate(np.reshape(np.array(sorted(cell_per_sample[sample])), orientation)):
                d[sample][list(string.ascii_uppercase)[j]] = e



def get_final_output():
    """
    Function called by snakemake rule all to run the pipeline
    """
    final_list = list()
    final_list.extend(
        expand(
            "{path}/{sample}/cell_selection/labels.tsv",
            path=config["data_location"],
            sample=samples,
        )
    )

    final_list.extend(
        (
            [
                sub_e
                for e in [
                    expand(
                        "{path}/{sample}/fastqc/{cell}_{pair}_fastqc.html",
                        path=config["data_location"],
                        sample=sample,
                        cell=cell_per_sample[sample],
                        pair=[1, 2],
                    )
                    for sample in samples
                ]
                for sub_e in e
            ]
        )
    )




    if config["GC_analysis"] is True:


        final_list.extend(
            (
                [
                    sub_e
                    for e in [
                        expand(
                            "{path}/{sample}/plots/alfred/{bam}_gc_{alfred_plot}.png",
                            path=config["data_location"],
                            sample=sample,
                            bam=cell_per_sample[sample],
                            alfred_plot=config["alfred_plots"]
                        )
                        for sample in samples
                    ]
                    for sub_e in e
                ]
            )
        )
        final_list.extend(
            (
                [
                    sub_e
                    for e in [
                        expand(
                            "{path}/{sample}/plots/alfred/MERGE/merged_bam_gc_{alfred_plot}.merge.png",
                            path=config["data_location"],
                            sample=sample,
                            alfred_plot=config["alfred_plots"]
                        )
                        for sample in samples
                    ]
                    for sub_e in e
                ]
            )
        )



        final_list.extend(
            (
                [
                    sub_e
                    for e in [
                        expand(
                            "{path}/{sample}/plots/alfred/PLATE_ROW/{row}_gc_{alfred_plot}.row.png",
                            path=config["data_location"],
                            sample=sample,
                            row=list(string.ascii_uppercase)[:orientation[0]],
                            alfred_plot=config["alfred_plots"]
                        )
                        for sample in samples if len(cell_per_sample[sample]) == 96
                    ]
                    for sub_e in e
                ]
            )
        )
        # final_list.extend(
        #     (
        #         [
        #             sub_e
        #             for e in [
        #                 expand(
        #                     "{path}/{sample}/plots/alfred/{row}_gc_{alfred_plot}.png",
        #                     path=config["data_location"],
        #                     sample=sample,
        #                     row=row,
        #                     alfred_plot=config["alfred_plots"]
        #                 )
        #                 for sample in samples for row in d[sample]
        #             ]
        #             for sub_e in e
        #         ]
        #     )
        # )
        # final_list.extend(
        #     expand(
        #         "{output_folder}/{sample}/plots/{sample}/alfred/gc_devi.png",
        #         output_folder=config["data_location"],
        #         sample=samples,
        #     ),
        # )
        # final_list.extend(
        #     expand(
        #         "{output_folder}/{sample}/plots/{sample}/alfred/gc_dist.png",
        #         output_folder=config["data_location"],
        #         sample=samples,
        #     ),
        # )

        final_list.extend(
            expand(
                "{output_folder}/{sample}/plots/counts/CountComplete.{plottype_counts}.pdf",
                output_folder=config["data_location"],
                sample=samples,
                plottype_counts=plottype_counts,
            ),
        )

    return final_list


def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [2, 4, 8, 16, 64, 128, 256]
    return mem_avail[attempt - 1] * 1000


def get_mem_mb_heavy(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [16, 64, 128, 256]
    return mem_avail[attempt - 1] * 1000
