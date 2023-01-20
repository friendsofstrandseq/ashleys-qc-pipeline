import os, sys
import pandas as pd
import scipy


# LOAD MOSAIC COUNTS INFO
# counts = snakemake.input.counts
# counts_df = pd.read_csv(counts, sep="\t", compression="gzip")
# counts_df["cell"] = counts_df["cell"] + ".sort.mdup.bam"

# # Groupby cell & sum reads
# counts_gb_df = counts_df.groupby("cell")[["c", "w"]].sum()
# counts_gb_df["nb_reads"] = counts_gb_df["c"] + counts_gb_df["w"]
info = snakemake.input.info
info_df = pd.read_csv(info, sep="\t", skiprows=13)
info_df["cell"] = info_df["cell"] + ".sort.mdup.bam"

# Load ashleys predictions
# labels = pd.read_csv(directory + "/cell_selection/labels_raw.tsv".format(sample=sample), sep="\t").sort_values(by="cell")
labels_path = snakemake.input.labels
labels = pd.read_csv(labels_path, sep="\t").sort_values(by="cell")
# labels["sample"] = sample
labels["sample"] = snakemake.wildcards.sample

# Retrieve correct prediction
labels_corrected = labels.loc[labels["prediction"] == 1]

# Merge counts & labels
# labels_corrected = pd.merge(info_df.reset_index()[["cell", "good"]], labels_corrected, on="cell")
labels_corrected = pd.merge(info_df[["cell", "good"]], labels_corrected, on="cell")
# Compute z-score on reads nb
labels_corrected["z_score"] = scipy.stats.zscore(labels_corrected["good"])

# Output, Correct outliers predictions & proba
labels_corrected.loc[labels_corrected["z_score"] >= 5, "cell"].to_csv(snakemake.output.bypass_cell, index=False, sep="\t")
labels_corrected.loc[labels_corrected["z_score"] >= 5, "new_prediction"] = 0
labels_corrected.loc[labels_corrected["z_score"] >= 5, "new_probability"] = 0
labels_corrected.loc[labels_corrected["z_score"] < 5, "new_probability"] = labels_corrected.loc[
    labels_corrected["z_score"] < 5, "probability"
]
labels_corrected["new_prediction"] = labels_corrected["new_prediction"].fillna(1)
labels_corrected["new_prediction"] = labels_corrected["new_prediction"].astype(int)

# Back to full dataframe
labels.loc[labels["cell"].isin(labels_corrected.cell.values.tolist()), "prediction"] = labels_corrected.new_prediction.values.tolist()
labels.loc[labels["cell"].isin(labels_corrected.cell.values.tolist()), "probability"] = labels_corrected.new_probability.values.tolist()

# Output
labels.to_csv(snakemake.output.labels_corrected, index=False, sep="\t")
