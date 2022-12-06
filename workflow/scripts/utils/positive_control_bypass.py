import os, sys
import pandas as pd
import scipy


# LOAD MOSAIC COUNTS
# counts_df = pd.read_csv(counts, sep="\t", compression="gzip")
counts = snakemake.input.counts
counts_df["cell"] = counts_df["cell"] + ".sort.mdup.bam"

# Groupby cell & sum reads
counts_gb_df = counts_df.groupby("cell")[["c","w"]].sum()
counts_gb_df["nb_reads"] = counts_gb_df["c"] + counts_gb_df["w"]

# Load ashleys predictions
# labels = pd.read_csv(directory + "/cell_selection/labels_raw.tsv".format(sample=sample), sep="\t").sort_values(by="cell")
labels_path = snakemake.input.labels
labels = pd.read_csv(labels_path, sep="\t").sort_values(by="cell")
labels["sample"] = sample

# Retrieve correct prediction
labels_corrected = labels.loc[labels["prediction"] == 1]

# Merge counts & labels
labels_corrected = pd.merge(counts_gb_df.reset_index()[["cell", "nb_reads"]], labels_corrected, on="cell")
# Compute z-score on reads nb
labels_corrected["z_score"] = scipy.stats.zscore(labels_corrected["nb_reads"])

# Correct outliers predictions & proba
labels_corrected.loc[labels_corrected["z_score"] >= 5, "new_prediction"] = 0
labels_corrected.loc[labels_corrected["z_score"] >= 5, "new_probability"] = 0
labels_corrected.loc[labels_corrected["z_score"] < 5, "new_probability"] = labels_corrected.loc[labels_corrected["z_score"] < 5, "probability"]
labels_corrected["new_prediction"] = labels_corrected["new_prediction"].fillna(1)
labels_corrected["new_prediction"] = labels_corrected["new_prediction"].astype(int)

# Back to full dataframe
labels.loc[labels["cell"].isin(labels_corrected.cell.values.tolist()), "prediction"] = labels_corrected.new_prediction.values.tolist()
labels.loc[labels["cell"].isin(labels_corrected.cell.values.tolist()), "probability"] = labels_corrected.new_probability.values.tolist()

# Output
labels.to_csv(snakemake.output.labels_corrected)
