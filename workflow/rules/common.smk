from workflow.scripts.utils import handle_input
import pandas as pd

# Create configuration file with samples
c = handle_input.HandleInput(
    input_path=config["folder"],
    output_path="{folder}/config/config_df_ashleys.tsv".format(folder=config["folder"]),
    check_sm_tag=False,
    bam=False
    )

# Read config file previously produced
df_config_files = pd.read_csv("{folder}/config/config_df_ashleys.tsv".format(folder=config["folder"]), sep="\t")

# List of available samples
samples = list(sorted(list(df_config_files.Sample.unique()())))

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

folder_expand = [
    [config["folder"]] * len(cell_per_sample[k])
    for k in cell_per_sample.keys()
]
folder_expand = [sub_e for e in folder_expand for sub_e in e]


def get_final_output():
    """
    Function called by snakemake rule all to run the pipeline
    """
    final_list = list()
    final_list.extend(expand("{path}/config/{sample}_selected_cells.ok", path=config["folder"], sample=samples,))
    final_list.extend(([sub_e for e in [expand("{path}/{sample}/fastqc/{cell}_{pair}_fastqc.html", path=config["folder"], sample=samples, cell=cell_per_sample[sample], pair=[1,2]) for sample in samples] for sub_e in e]))
    return final_list

def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [1, 2, 4, 8, 16, 64, 128, 256]
    return mem_avail[attempt - 1] * 1000
