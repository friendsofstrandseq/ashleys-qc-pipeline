import pandas as pd

# Read 200kb bins file

binbed = pd.read_csv(
    snakemake.input.bin_bed,
    # "../../../../mosaicatcher-update/workflow/data/bin_200kb_all.bed",
    sep="\t",
    names=["chrom", "start", "end", "bin_id"],
)
binbed["ID"] = (
    binbed["chrom"].astype(str)
    + "_"
    + binbed["start"].astype(str)
    + "_"
    + binbed["end"].astype(str)
)

# Turn chrom into categorical
binbed["chrom"] = pd.Categorical(
    binbed["chrom"],
    categories=["chr{}".format(e) for e in range(1, 23)] + ["chrX", "chrY"],
    ordered=True,
)

# Sort & filter out chrY #TMP / can be changed
binbed = binbed.sort_values(by=["chrom", "start", "end"]).reset_index(drop=True)
binbed["w"], binbed["c"], binbed["class"] = 0, 0, None


# Determine if using mm10 data and load normalization file if so
reference = snakemake.config["reference"]
if reference == "mm10":
    print(reference)
    norm_file = "workflow/data/normalization/mm10/HGSVC.200000.txt"
    norm_df = pd.read_csv(norm_file, sep="\t")
    norm_df["ID"] = (
        norm_df["chrom"].astype(str)
        + "_"
        + norm_df["start"].astype(str)
        + "_"
        + norm_df["end"].astype(str)
    )
    norm_df["class"] = norm_df["class"].fillna("None")
    norm_df["class"] = norm_df["class"].astype(str)
else:
    norm_df = pd.DataFrame()

# Read SV file
# df = pd.read_csv("../../../../mosaicatcher-update/.tests/data_CHR17/RPE-BM510/counts/RPE-BM510.txt.raw.gz", sep="\t")

# sep = "," if "/multistep_normalisation/" in snakemake.input.counts else "\t"
sep = "\t"
df = pd.read_csv(snakemake.input.counts, sep=sep, compression="gzip")
df["ID"] = (
    df["chrom"].astype(str)
    + "_"
    + df["start"].astype(str)
    + "_"
    + df["end"].astype(str)
)
df["w"] = df["w"].round(0).astype(int)
df["c"] = df["c"].round(0).astype(int)
if sep == ",":
    df["tot_count"] = df["tot_count"].round(0).astype(int)

# Adjust counts based on normalization data if applicable
if not norm_df.empty:
    df = df.merge(norm_df[["ID", "class"]], on="ID", how="left")
    df.loc[df["class_y"] == "None", ["w", "c"]] = 0
    df.drop(columns=["class_y"], inplace=True)
    df = df.rename(columns={"class_x": "class"})
    df["class"] = df["class"].fillna("None")
    print(df)


## Populate counts df for each cell in order to have all bins represented
l = list()

# Loop over cells
for cell in df.cell.unique().tolist():
    # Outer join to retrieve both real count values from specified chromosome and empty bins
    tmp_df = pd.concat(
        [
            binbed.loc[
                ~binbed["ID"].isin(df.loc[df["cell"] == cell].ID.values.tolist())
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
# populated_df.to_csv("test.txt.gz", compression="gzip", sep="\t", index=False)
populated_df.to_csv(
    snakemake.output.populated_counts, compression="gzip", sep="\t", index=False
)
