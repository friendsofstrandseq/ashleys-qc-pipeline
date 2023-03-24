import os, sys
import pandas as pd
import scipy

# LOAD MOSAIC COUNTS INFO
info = snakemake.input.info
info_df = pd.read_csv(info, sep="\t", skiprows=13)
info_df["cell"] = info_df["cell"] + ".sort.mdup.bam"

# Load ashleys predictions
labels_path = snakemake.input.labels
labels = pd.read_csv(labels_path, sep="\t").sort_values(by="cell")
labels["sample"] = snakemake.wildcards.sample

# Retrieve correct prediction
labels_corrected = labels.loc[labels["prediction"] == 1]

# Merge counts & labels
labels_corrected = pd.merge(info_df[["cell", "good"]], labels_corrected, on="cell")
# Compute z-score on reads nb
labels_corrected["z_score"] = scipy.stats.zscore(labels_corrected["good"])

print(labels_corrected )

# Output, Correct outliers predictions & proba
z_score_cutoff = 5
labels_corrected.loc[labels_corrected["z_score"] >= z_score_cutoff, "cell"].to_csv(snakemake.output.bypass_cell, index=False, sep="\t")
labels_corrected.loc[labels_corrected["z_score"] >= z_score_cutoff, "new_prediction"] = 0
labels_corrected.loc[labels_corrected["z_score"] >= z_score_cutoff, "new_probability"] = 0
labels_corrected.loc[labels_corrected["z_score"] < z_score_cutoff, "new_probability"] = labels_corrected.loc[
    labels_corrected["z_score"] < z_score_cutoff, "probability"
]
labels_corrected["new_prediction"] = labels_corrected["new_prediction"].fillna(1)
labels_corrected["new_prediction"] = labels_corrected["new_prediction"].astype(int)

# Back to full dataframe
labels.loc[labels["cell"].isin(labels_corrected.cell.values.tolist()), "prediction"] = labels_corrected.new_prediction.values.tolist()
labels.loc[labels["cell"].isin(labels_corrected.cell.values.tolist()), "probability"] = labels_corrected.new_probability.values.tolist()

# Output
labels.to_csv(snakemake.output.labels_corrected, index=False, sep="\t")
