#!/usr/bin/env python3

import pandas as pd


############################################################
# Input / output
############################################################

SRA_FILE = "SRA_characteristics.tsv"
ANNOTATION_FILE = "biosample_metadata_infection_v8.tsv"

OUTPUT_FILE = "SRA_metadata_annotated_v2.tsv"


############################################################
# Load
############################################################

print("Loading SRA metadata...")

sra = pd.read_csv(
    SRA_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)

print("SRA rows:", len(sra))
print("Unique runs:", sra["Run_accession"].nunique())
print("Unique BioSamples:", sra["BioSample"].nunique())


print("\nLoading BioSample annotations...")

anno = pd.read_csv(
    ANNOTATION_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)

print("Annotation rows:", len(anno))
print("Unique BioSamples:", anno["BioSample"].nunique())


############################################################
# Check BioSample uniqueness
############################################################

duplicated = anno["BioSample"].duplicated(keep=False)

if duplicated.any():

    print("\nWARNING: duplicated BioSamples in annotation file!")

    print(
        anno.loc[duplicated, "BioSample"]
        .value_counts()
        .head(20)
    )

    raise ValueError(
        "BioSample is not unique in annotation file."
    )


############################################################
# Select annotation columns
#
# Do not duplicate columns already present in SRA file.
############################################################

annotation_columns = [
    col
    for col in anno.columns
    if col == "BioSample"
    or col not in sra.columns
]

anno = anno[annotation_columns]


############################################################
# Merge
############################################################

print("\nMerging...")

merged = sra.merge(
    anno,
    on="BioSample",
    how="left",
    validate="many_to_one",
    indicator=True
)


############################################################
# QC
############################################################

print("\n=== Merge status ===")

print(
    merged["_merge"].value_counts()
)

print("\nInput SRA rows :", len(sra))
print("Output rows    :", len(merged))

if len(sra) != len(merged):
    raise ValueError(
        "Number of rows changed during merge!"
    )


############################################################
# Show unmatched BioSamples
############################################################

unmatched = merged[
    merged["_merge"] == "left_only"
]

print("\nUnmatched SRA rows:", len(unmatched))
print(
    "Unmatched unique BioSamples:",
    unmatched["BioSample"].nunique()
)

if len(unmatched) > 0:

    unmatched[
        [
            "Run_accession",
            "Sample_accession",
            "BioSample"
        ]
    ].drop_duplicates().to_csv(
        "SRA_unmatched_BioSamples.tsv",
        sep="\t",
        index=False
    )

    print(
        "Saved unmatched BioSamples: "
        "SRA_unmatched_BioSamples.tsv"
    )


############################################################
# Remove merge indicator
############################################################

merged = merged.drop(
    columns="_merge"
)


############################################################
# Save
############################################################

merged.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)

print(
    "\nSaved:",
    OUTPUT_FILE
)
