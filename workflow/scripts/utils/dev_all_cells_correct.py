import pandas as pd

df = pd.read_csv(snakemake.input.folder, sep="\t")
df["prediction"] = 1
df["probability"] = 1
df.loc[df["cell"].str.contains("BM510x04_PE20305"), "prediction"] = 0
df.loc[df["cell"].str.contains("BM510x04_PE20312"), "prediction"] = 0
df.to_csv(snakemake.output.folder, sep="\t", index=False)