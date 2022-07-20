import pandas as pd
import subprocess
import os
df = pd.read_csv(snakemake.input.predictions, sep="\t")
cells_unselected = df.loc[df["prediction"] == 0, "cell"].tolist()

# ADDING NEW COLUMN TO CONFIG FILE
df_config = pd.read_csv("{folder}/config/config_df_ashleys.tsv".format(folder=snakemake.config["folder"]), sep='\t')
df_config["Selected"] = True
df_config.loc[df_config["Cell"].isin([e.split('.')[0] for e in cells_unselected]), "Selected"] = False
df_config.to_csv("{folder}/config/config_df_ashleys.tsv".format(folder=snakemake.config["folder"]), sep='\t', index=False)

with open(snakemake.output[0], 'w') as out:
    out.write("Folder processed: {folder}\n".format(folder=snakemake.params.path))
    out.write("Removed following cells:\n")
    for cell in cells_unselected:
        # print("rm {path}/{sample}/selected/{cell}".format(path=snakemake.params.path, sample=snakemake.wildcards.sample, cell=cell))
        subprocess.call("rm {path}/{sample}/selected/{cell}".format(path=snakemake.params.path, sample=snakemake.wildcards.sample, cell=cell), shell=True)
        subprocess.call("rm {path}/{sample}/selected/{cell}.bai".format(path=snakemake.params.path, sample=snakemake.wildcards.sample, cell=cell), shell=True)
        out.write("- {cell}\n".format(cell=cell))