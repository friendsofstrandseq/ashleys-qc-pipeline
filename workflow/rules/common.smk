import pandas as pd
import os, sys


# MOVED OUTSIDE A SCRIPT TO PREVENT PATH ISSUES
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
        folder = "all" if bam is True else "fastq"
        complete_df_list = list()
        # print(thisdir)
        for sample in [
            e
            for e in os.listdir(thisdir)
            if e not in ["config", "log", ".DS_Store", "._.DS_Store"]
        ]:
            # print("{thisdir}/{sample}/{folder}/".format(thisdir=thisdir, sample=sample, folder=folder))
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
            if bam is True:
                l_files_selected = [
                    f
                    for f in os.listdir(
                        "{thisdir}/{sample}/selected/".format(
                            thisdir=thisdir, sample=sample
                        )
                    )
                    if f.endswith(".bam")
                ]
                print(l_files_selected)

                join = list(set(l_files_all).intersection(set(l_files_selected)))
                df["Selected"] = False
                df.loc[df["File"].isin(join), "Selected"] = True

            complete_df_list.append(df)

        complete_df = pd.concat(complete_df_list)
        complete_df = complete_df.sort_values(by=["Cell", "File"]).reset_index(
            drop=True
        )
        # complete_df = complete_df.loc[~complete_df["Cell"].isin(exclude_list)]
        return complete_df


# Create configuration file with samples

c = HandleInput(
    input_path=config["input_bam_location"],
    output_path="{input_bam_location}/config/config_df_ashleys.tsv".format(
        input_bam_location=config["input_bam_location"]
    ),
    check_sm_tag=False,
    bam=False,
)
df_config_files = c.df_config_files

# Read config file previously produced
# df_config_files = pd.read_csv("{input_bam_location}/config/config_df_ashleys.tsv".format(input_bam_location=config["input_bam_location"]), sep="\t")
# print(df_config_files)
# exit()
# List of available samples
samples = list(sorted(list(df_config_files.Sample.unique().tolist())))

# Dictionnary of libraries available per sample
cell_per_sample = (
    df_config_files.groupby("Sample")["Cell"].unique().apply(list).to_dict()
)


# FIXME: this is a workaround used to use *zip* function inside expand statements
samples_expand = [[k] * len(cell_per_sample[k]) for k in cell_per_sample.keys()]
samples_expand = [sub_e for e in samples_expand for sub_e in e]

cell_expand = [sub_e for e in list(cell_per_sample.values()) for sub_e in e]

input_bam_location_expand = [
    [config["input_bam_location"]] * len(cell_per_sample[k])
    for k in cell_per_sample.keys()
]
input_bam_location_expand = [sub_e for e in input_bam_location_expand for sub_e in e]


def get_final_output():
    """
    Function called by snakemake rule all to run the pipeline
    """
    final_list = list()
    final_list.extend(
        expand(
            "{path}/{sample}/cell_selection/labels.tsv",
            path=config["input_bam_location"],
            sample=samples,
        )
    )
    # final_list.extend(expand("{path}/config/{sample}_selected_cells.ok", path=config["input_bam_location"], sample=samples,))
    final_list.extend(
        (
            [
                sub_e
                for e in [
                    expand(
                        "{path}/{sample}/fastqc/{cell}_{pair}_fastqc.html",
                        path=config["input_bam_location"],
                        sample=samples,
                        cell=cell_per_sample[sample],
                        pair=[1, 2],
                    )
                    for sample in samples
                ]
                for sub_e in e
            ]
        )
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
