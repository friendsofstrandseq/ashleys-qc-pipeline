import pandas as pd

# Read 200kb bins file

binbed = pd.read_csv(
    snakemake.input.bin_bed[0],
    # "../../../../mosaicatcher-update/workflow/data/bin_200kb_all.bed",
    sep="\t",
    # names=["chrom", "start", "end", "scalar", "class"],
)
binbed["bin_id"] = (
    binbed["chrom"].astype(str)
    + "_"
    + binbed["start"].astype(str)
    + "_"
    + binbed["end"].astype(str)
)

binbed = binbed.drop(columns=["scalar"])


# Turn chrom into categorical
max_chrom = 23 if snakemake.config["reference"] != "mm10" else 20
binbed["chrom"] = pd.Categorical(
    binbed["chrom"],
    categories=["chr{}".format(e) for e in range(1, 23)] + ["chrX", "chrY"],
    ordered=True,
)

# Sort & filter out chrY #TMP / can be changed
binbed = binbed.sort_values(by=["chrom", "start", "end"]).reset_index(drop=True)
binbed["w"], binbed["c"] = 0, 0

# Read SV file
# df = pd.read_csv("../../../../mosaicatcher-update/.tests/data_CHR17/RPE-BM510/counts/RPE-BM510.txt.raw.gz", sep="\t")

# sep = "," if "/multistep_normalisation/" in snakemake.input.counts else "\t"
sep = "\t"
df = pd.read_csv(snakemake.input.counts, sep=sep, compression="gzip")
df["bin_id"] = (
    df["chrom"].astype(str)
    + "_"
    + df["start"].astype(str)
    + "_"
    + df["end"].astype(str)
)
df["w"] = df["w"].round(0).astype(int)
df["c"] = df["c"].round(0).astype(int)
df["class"] = df["class"].astype(str)
if sep == ",":
    df["tot_count"] = df["tot_count"].round(0).astype(int)


# Update 'class' in df using norm_df
df = df.merge(binbed[["bin_id", "class"]], on="bin_id", how="left")
df["class"] = df.apply(
    lambda x: "None" if x["class_y"] != "good" else x["class_x"], axis=1
)

# set w & c to 0 for None class if None in class_y
df.loc[df["class_y"] == "None", "w"] = 0
df.loc[df["class_y"] == "None", "c"] = 0


df = df.drop(columns=["class_x", "class_y"])
# df.loc[df["class"] == "None", "w"] = 0
# df.loc[df["class"] == "None", "c"] = 0

## Populate counts df for each cell in order to have all bins represented
l = list()

# Loop over cells
for cell in df.cell.unique().tolist()[1:6]:
    # Outer join to retrieve both real count values from specified chromosome and empty bins
    tmp_df = pd.concat(
        [
            binbed.loc[
                ~binbed["bin_id"].isin(
                    df.loc[df["cell"] == cell].bin_id.values.tolist()
                )
            ],
            df.loc[df["cell"] == cell],
        ]
    )

    # Filla cell & sample columns
    tmp_df["cell"] = cell
    tmp_df["sample"] = df.loc[df["cell"] == cell, "sample"].values.tolist()[0]
    l.append(tmp_df)

# Concat list of DF and output
populated_df = pd.concat(l).sort_values(by=["cell", "chrom", "start"])
print(populated_df)
populated_df["start"] = populated_df["start"].astype(int)
populated_df["end"] = populated_df["end"].astype(int)
populated_df["w"] = populated_df["w"].astype(int)
populated_df["c"] = populated_df["c"].astype(int)
# populated_df.to_csv("test.txt.gz", compression="gzip", sep="\t", index=False)
populated_df.to_csv(
    snakemake.output.populated_counts, compression="gzip", sep="\t", index=False
)
