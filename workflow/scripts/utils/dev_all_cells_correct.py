import pandas as pd

df = pd.read_csv(snakemake.input.folder, sep="\t")
df["prediction"] = 1
df["probability"] = 1
df.loc[df["cell"].str.contains("05"), "prediction"] = 0
df.loc[df["cell"].str.contains("12"), "prediction"] = 0
df.to_csv(snakemake.output.path, sep="\t", index=False)