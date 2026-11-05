#!/usr/bin/env python3

import re
import pandas as pd


# ============================================================
# Settings
# ============================================================

INPUT_FILE = "biosample_characteristics.tsv"
OUTPUT_FILE = "biosample_metadata.tsv"


# ============================================================
# Attribute keys used for Tissue extraction
# ============================================================

TISSUE_KEYS = [
    "tissue",
    "tissue_type",
    "organ",
    "organism_part",
    "body_site",
    "anatomical_site",
    "source_tissue",
    # added based on uncharacterized samples
    #"source_name",
    #"isolation_source",
    "cell/tissue_type",
    "tissue_or_dev_stage",
]

# ============================================================
# Attribute keys used for Cell type extraction
# ============================================================

CELL_KEYS = [
    "cell_type",
    "celltype",
    "cell_line",
    "cell_line_name",
    # added based on uncharacterized samples
    "cell_or_tissue_type",
]

# ============================================================
# Normalize attribute-key names
#
# Examples:
#
# "organism part" -> "organism_part"
# "Organism Part" -> "organism_part"
# "cell-type"     -> "cell_type"
#
# ============================================================

def normalize_key(key):

    key = str(key).strip().lower()

    key = re.sub(
        r"[\s\-]+",
        "_",
        key
    )

    key = re.sub(
        r"_+",
        "_",
        key
    )

    return key


# ============================================================
# Pre-normalize key lists
# ============================================================

TISSUE_KEYS_NORMALIZED = {
    normalize_key(x)
    for x in TISSUE_KEYS
}

CELL_KEYS_NORMALIZED = {
    normalize_key(x)
    for x in CELL_KEYS
}


# ============================================================
# Parse All_characteristics
#
# Input example:
#
# organism=Gallus gallus |
# organism part=bone marrow |
# description=Bone marrow sample
#
# Output:
#
# [
#   ("organism", "Gallus gallus"),
#   ("organism part", "bone marrow"),
#   ("description", "Bone marrow sample")
# ]
#
# ============================================================

def parse_characteristics(text):

    if pd.isna(text):
        return []

    text = str(text).strip()

    if not text:
        return []

    attributes = []

    for item in text.split(" | "):

        item = item.strip()

        if not item:
            continue

        if "=" not in item:
            continue

        # Split only at the first "="
        #
        # Important because values themselves may contain "=".
        key, value = item.split(
            "=",
            1
        )

        key = key.strip()
        value = value.strip()

        if not key:
            continue

        if not value:
            continue

        attributes.append(
            (
                key,
                value
            )
        )

    return attributes


# ============================================================
# Extract attributes matching a predefined set of keys
#
# Returns:
#
# value_string
# source_string
#
# Example:
#
# organism part=bone marrow
# tissue=bone marrow
#
# ->
#
# value:
# bone marrow
#
# source:
# organism part;tissue
#
# Duplicate values are removed while preserving order.
# ============================================================

def extract_attributes(
    attributes,
    target_keys
):

    values = []
    sources = []

    for key, value in attributes:

        key_norm = normalize_key(
            key
        )

        if key_norm not in target_keys:
            continue

        # ----------------------------------------------------
        # A single BioSample characteristic can contain
        # multiple values joined by "; ".
        #
        # Preserve these as separate candidate values.
        # ----------------------------------------------------

        split_values = [
            x.strip()
            for x in value.split(";")
            if x.strip()
        ]

        for x in split_values:

            if x not in values:
                values.append(x)

        if key not in sources:
            sources.append(key)


    if values:
        value_string = "; ".join(
            values
        )
    else:
        value_string = "Unknown"


    if sources:
        source_string = "; ".join(
            sources
        )
    else:
        source_string = ""


    return (
        value_string,
        source_string
    )


# get source_name, isolation_source independently
def extract_single_key(attributes, target_key):

    values = []

    target_key = normalize_key(target_key)

    for key, value in attributes:

        if normalize_key(key) == target_key:

            if value not in values:
                values.append(value)

    if values:
        return "; ".join(values)

    return ""

# ============================================================
# Process one BioSample
# ============================================================

def extract_metadata(row):

    all_characteristics = row.get(
        "All_characteristics",
        ""
    )

    attributes = parse_characteristics(
        all_characteristics
    )


    tissue, tissue_source = (
        extract_attributes(
            attributes,
            TISSUE_KEYS_NORMALIZED
        )
    )


    cell_type, cell_type_source = (
        extract_attributes(
            attributes,
            CELL_KEYS_NORMALIZED
        )
    )

    source_name = extract_single_key(
        attributes,
        "source_name"
    )
    
    isolation_source = extract_single_key(
        attributes,
        "isolation_source"
    )

    return pd.Series({
        "Tissue": tissue,
        "Cell_type": cell_type,
        "Tissue_source": tissue_source,
        "Cell_type_source": cell_type_source,
        "Source_name": source_name,
        "Isolation_source": isolation_source
    })


# ============================================================
# Main
# ============================================================

def main():

    print(
        f"Reading: {INPUT_FILE}"
    )


    df = pd.read_csv(
        INPUT_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False
    )


    # --------------------------------------------------------
    # Required columns
    # --------------------------------------------------------

    required_columns = {
        "BioSample",
        "All_characteristics",
        "BioSample_status"
    }


    missing_columns = (
        required_columns
        - set(df.columns)
    )


    if missing_columns:

        raise ValueError(
            "Missing required columns: "
            + ", ".join(
                sorted(
                    missing_columns
                )
            )
        )


    print(
        "Number of BioSamples:",
        len(df)
    )


    # --------------------------------------------------------
    # Extract Tissue and Cell type
    # --------------------------------------------------------

    print(
        "Extracting Tissue / Cell type..."
    )


    extracted = df.apply(
        extract_metadata,
        axis=1
    )


    # --------------------------------------------------------
    # Construct output table
    # --------------------------------------------------------

    out = pd.concat(
        [
            df[
                [
                    "BioSample"
                ]
            ],

            extracted,

            df[
                [
                    "All_characteristics",
                    "BioSample_status"
                ]
            ]
        ],
        axis=1
    )


    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    out.to_csv(
        OUTPUT_FILE,
        sep="\t",
        index=False
    )


    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print(
        "=== Tissue extraction ==="
    )


    tissue_known = (
        out["Tissue"]
        != "Unknown"
    ).sum()


    print(
        "Known:",
        tissue_known
    )

    print(
        "Unknown:",
        len(out)
        - tissue_known
    )


    print()
    print(
        "=== Tissue sources ==="
    )


    print(
        out.loc[
            out["Tissue_source"]
            != "",
            "Tissue_source"
        ]
        .value_counts()
        .head(30)
        .to_string()
    )


    print()
    print(
        "=== Cell type extraction ==="
    )


    cell_known = (
        out["Cell_type"]
        != "Unknown"
    ).sum()


    print(
        "Known:",
        cell_known
    )

    print(
        "Unknown:",
        len(out)
        - cell_known
    )


    print()
    print(
        "=== Cell type sources ==="
    )


    print(
        out.loc[
            out["Cell_type_source"]
            != "",
            "Cell_type_source"
        ]
        .value_counts()
        .head(30)
        .to_string()
    )


    print()
    print(
        "Saved:"
    )

    print(
        OUTPUT_FILE
    )


if __name__ == "__main__":
    main()
