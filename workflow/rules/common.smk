from scripts.utils import handle_input
import pandas as pd

# Create configuration file with samples

c = handle_input.HandleInput(
    input_path=config["input_bam_location"],
    output_path="{input_bam_location}/config/config_df_ashleys.tsv".format(input_bam_location=config["input_bam_location"]),
    check_sm_tag=False,
    bam=False
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
    df_config_files
    .groupby("Sample")["Cell"]
    .unique().apply(list)
    .to_dict()
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
    final_list.extend(expand("{path}/{sample}/cell_selection/labels.tsv", path=config["input_bam_location"], sample=samples,))
    # final_list.extend(expand("{path}/config/{sample}_selected_cells.ok", path=config["input_bam_location"], sample=samples,))
    final_list.extend(([sub_e for e in [expand("{path}/{sample}/fastqc/{cell}_{pair}_fastqc.html", path=config["input_bam_location"], sample=samples, cell=cell_per_sample[sample], pair=[1,2]) for sample in samples] for sub_e in e]))
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
