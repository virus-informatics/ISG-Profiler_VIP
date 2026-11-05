#!/usr/bin/env python3

import re
import pandas as pd


############################################################
# 1. Input / output
############################################################

# Tissue normalization済みファイルを入力にする
INPUT_FILE = "biosample_metadata_tissue_normalized_v2.tsv"

OUTPUT_FILE = "biosample_metadata_celltype_normalized_v3.tsv"

UNCLASSIFIED_FILE = "celltype_unclassified_counts.tsv"
MULTIPLE_FILE = "celltype_multiple_counts.tsv"
CELLLINE_FILE = "cellline_counts.tsv"
RESCUED_FILE = "celltype_rescued.tsv"

############################################################
# 2. Cell type normalization dictionary
############################################################

CELL_TYPE_PATTERNS = {

    "Fibroblast": [
        r"\bfibroblasts?\b",
        r"\bfibroblastic\b",
        r"\bfibroblast[- ]like\b",
    ],

    "PBMC": [
        r"\bpbmcs?\b",
        r"\bperipheral blood mononuclear cells?\b",
    ],

    "Mononuclear cell": [
        r"\bmononuclear cells?\b",
    ],

    "T cell": [
        r"(?<![-\w])t[ _-]?cells?\b",
        r"\bcd3\+?[ _-]?t[ _-]?cells?\b",
        r"\bcd4\+?[ _-]?t[ _-]?cells?\b",
        r"\bcd8\+?[ _-]?t[ _-]?cells?\b",
        r"\bcd4\+\b",
        r"\bcd8\+\b",
    ],

    "B cell": [
        r"\bb[ _-]?cells?\b",
        r"\bb[ _-]?lymphocytes?\b",
        r"\bcd19\+?[ _-]?cells?\b",
        r"\bcd20\+?[ _-]?cells?\b",
    ],

    "NK cell": [
        r"\bnatural killer cells?\b",
        r"\bnk[ _-]?cells?\b",
    ],

    "Macrophage": [
        r"\bmacrophages?\b",
        r"\balveolar macrophages?\b",
        r"\bbone marrow[- ]derived macrophages?\b",
        r"\bbmdms?\b",
    ],

    "Monocyte": [
        r"\bmonocytes?\b",
    ],

    "Dendritic cell": [
        r"\bdendritic cells?\b",
        r"\bdcs?\b",
    ],

    "Microglia": [
        r"\bmicroglia\b",
        r"\bmicroglial\b",
    ],

    "Leukocyte": [
        r"\bleukocytes?\b",
        r"\bwhite blood cells?\b",
        r"\bwhite blood cell\b",
    ],

    "Epithelial cell": [
        r"\bepithelial cells?\b",
        r"\bepithelial\b",
        r"\bepithelium\b",
    ],

    "Endothelial cell": [
        r"\bendothelial cells?\b",
        r"\bendothelium\b",
    ],

    "Cardiomyocyte": [
        r"\bcardiomyocytes?\b",
    ],

    "Myoblast": [
        r"\bmyoblasts?\b",
    ],

    "Myocyte": [
        r"\bmyocytes?\b",
    ],

    "Neuron": [
        r"\bneurons?\b",
        r"\bneuronal\b",
    ],

    "Astrocyte": [
        r"\bastrocytes?\b",
    ],

    "Oligodendrocyte": [
        r"\boligodendrocytes?\b",
    ],

    "Granulosa cell": [
        r"\bgranulosa(?:[ _-]?cells?)?\b",
    ],

    "Theca cell": [
        r"\btheca(?:[ _-]?cells?)?\b",
    ],

    "Hepatocyte": [
        r"\bhepatocytes?\b",
    ],

    "Thyrocyte": [
        r"\bthyrocytes?\b",
    ],

    "Adipocyte": [
        r"\badipocytes?\b",
    ],

    "Mesenchymal stem cell": [
        r"\bmesenchymal stem cells?\b",
        r"\bmesenchymal stromal cells?\b",
    ],

    "iPSC": [
        r"\binduced pluripotent stem cells?\b",
        r"\bipscs?\b",
        r"\biPS cells?\b",
    ],

    "Embryonic stem cell": [
        r"\bembryonic stem cells?\b",
        r"\bescs?\b",
    ],

    "Neural crest cell": [
        r"\bneural crest(?:[ _-]?cells?)?\b",
    ],

    "Oocyte": [
        r"\boocytes?\b",
    ],

    "Blastomere": [
        r"\bblastomeres?\b",
    ],

    "Inner cell mass": [
        r"\binner cell mass(?: cell)?\b",
        r"\bicm cells?\b",
    ],

    "Mononuclear phagocyte": [
        r"\bmononuclear phagocytes?\b",
    ],
    
    "Osteoblast": [
        r"\bosteoblasts?\b",
    ],
}


############################################################
# 3. Known cell lines
#
# Cell lines are kept separately from biological cell types.
############################################################

CELL_LINE_PATTERNS = {

    "Vero E6": [
        r"\bvero[ _-]?e6\b",
        r"\bveroe6\b",
    ],

    "Vero": [
        r"^vero$",
        r"\bvero cells?\b",
    ],

    "MDCK": [
        r"\bmdck\b",
        r"\bg[- ]?mdck\b",
        r"\bmadin[- ]darby canine kidney\b",
    ],

    "COS-7": [
        r"\bcos[- ]?7\b",
    ],

    "CHO-S": [
        r"\bcho[- ]?s\b",
    ],

    "CHO": [
        r"\bcho cells?\b",
    ],

    "IPEC-J2": [
        r"\bipec[- ]?j2\b",
    ],

    "GH4C1": [
        r"\bgh4c1\b",
    ],

    "HEK293": [
        r"\bhek[- ]?293\b",
        r"\b293t\b",
        r"\bhek293t\b",
    ],

    "HeLa": [
        r"\bhela\b",
    ],

    "A549": [
        r"\ba549\b",
    ],

    "RAW264.7": [
        r"\braw[ _-]?264\.?7\b",
    ],

    "DF-1": [
        r"\bdf[- ]?1\b",
    ],

    "MAC-T": [
        r"\bmac[- ]?t\b",
    ],    

    # additional lines identified by QC
    "PCCl3": [
        r"\bpccl3\b",
    ],

    "H9C2": [
        r"\bh9c[- ]?2\b",
        r"\bh9c2(?:\(2-1\))?\b",
        r"\bh9c2[-_][\w-]+\b",
    ],

    "PK-15": [
        r"\bpk[- ]?15\b",
        r"\bpk15[a-z0-9_-]*\b",
    ],

    "LMH": [
        r"\blmh\b",
    ],

    "MARC-145": [
        r"\bmarc[- ]?145\b",
    ],

    "BRL-3A": [
        r"\bbrl[- ]?3a\b",
    ],

    "CHO-K1": [
        r"\bcho[- ]?k1\b",
        r"\bcho k1\b",
    ],

    "MDBK": [
        r"\bmdbk\b",
        r"\bmadin[- ]darby bovine kidney\b",
    ],

    "BHK21": [
        r"\bbhk[- ]?21\b",
    ],

    "RLE-6TN": [
        r"\brle[- ]?6tn\b",
    ],

    "HSC-T6": [
        r"\bhsc[- ]?t6\b",
    ],

    "CC-531": [
        r"\bcc[- ]?531\b",
    ],

    "3D4/21": [
        r"\b3d4/21\b",
    ],

    "BV2": [
        r"\bbv2\b",
    ],

    "PC12": [
        r"\bpc[- ]?12\b",
        r"\bpc12\b",
    ],
}


############################################################
# 4. Values treated as missing / uninformative
############################################################

MISSING_VALUES = {
    "",
    "unknown",
    "missing",
    "not collected",
    "not applicable",
    "not provided",
    "not specified",
    "na",
    "n/a",
    "none",
    "null",
}


############################################################
# 5. Compile regex once
############################################################

COMPILED_CELL_TYPE_PATTERNS = {
    category: [
        re.compile(pattern, flags=re.IGNORECASE)
        for pattern in patterns
    ]
    for category, patterns in CELL_TYPE_PATTERNS.items()
}


COMPILED_CELL_LINE_PATTERNS = {
    cell_line: [
        re.compile(pattern, flags=re.IGNORECASE)
        for pattern in patterns
    ]
    for cell_line, patterns in CELL_LINE_PATTERNS.items()
}


############################################################
# 6. Helper
############################################################

def find_hits(text, pattern_dict):

    hits = []

    for category, patterns in pattern_dict.items():

        if any(pattern.search(text) for pattern in patterns):
            hits.append(category)

    return list(dict.fromkeys(hits))

def is_missing_raw(value):

    text = str(value).strip().lower()

    return text in MISSING_VALUES

def extract_attribute_values(text, target_keys):

    if not text:
        return []

    values = []

    for item in text.split(" | "):

        if "=" not in item:
            continue

        key, value = item.split("=", 1)

        key_norm = re.sub(
            r"[\s_-]+",
            " ",
            key.strip().lower()
        )

        if key_norm in target_keys:

            value = value.strip()

            if value and not is_missing_raw(value):
                values.append(value)

    return values

############################################################
# Cell type / Cell line normalization with rescue
#
# Priority:
#
# 1. Explicit Cell_type field
#
# 2. Original Tissue field
#    ONLY when Tissue_status == "No_match"
#
# 3. Source_name
#
# 4. Isolation_source
#
# Explicit Cell_type annotations always have priority.
############################################################

############################################################
# 7. Basic Cell type / Cell line classification
#
# This function classifies ONE text string.
############################################################

def normalize_celltype(value):

    ########################################################
    # Missing
    ########################################################

    if is_missing_raw(value):

        return (
            "Unknown",
            "",
            "",
            "Missing"
        )


    text = str(value).strip()


    ########################################################
    # Detect biological cell types
    ########################################################

    cell_hits = find_hits(
        text,
        COMPILED_CELL_TYPE_PATTERNS
    )


    ########################################################
    # Detect known cell lines
    ########################################################

    cell_line_hits = find_hits(
        text,
        COMPILED_CELL_LINE_PATTERNS
    )


    ########################################################
    # Resolve known overlaps
    ########################################################

    # PBMC is more specific than generic Mononuclear cell
    if (
        "PBMC" in cell_hits
        and "Mononuclear cell" in cell_hits
    ):
        cell_hits.remove(
            "Mononuclear cell"
        )


    # PBMC-derived macrophages
    if (
        "PBMC" in cell_hits
        and "Macrophage" in cell_hits
        and re.search(
            r"\bpbmc[- ]derived macrophages?\b",
            text,
            flags=re.IGNORECASE
        )
    ):
        cell_hits.remove(
            "PBMC"
        )

    if (
        "PBMC" in cell_hits
        and len(cell_hits) > 1
        and re.search(
            r"\b(?:from|derived from)\s+PBMCs?\b",
            text,
            flags=re.IGNORECASE
        )
    ):
        cell_hits.remove("PBMC")

    # Monocyte-derived macrophages
    if (
        "Monocyte" in cell_hits
        and "Macrophage" in cell_hits
        and re.search(
            r"\bmonocyte[- ]derived\b",
            text,
            flags=re.IGNORECASE
        )
    ):
        cell_hits.remove(
            "Monocyte"
        )

    ########################################################
    # Resolve differentiation / derivation relationships
    ########################################################
    
    # ------------------------------------------------------
    # X-derived iPSCs
    #
    # The final cell type is iPSC.
    #
    # Example:
    #   iPSCs derived from PBMCs
    #       -> iPSC
    # ------------------------------------------------------
    
    if (
        "iPSC" in cell_hits
        and re.search(
            r"\b(?:ipscs?|induced pluripotent stem cells?)"
            r"\s+derived\s+from\b",
            text,
            flags=re.IGNORECASE
        )
    ):
        cell_hits = ["iPSC"]
    
    
    # ------------------------------------------------------
    # iPSC-derived differentiated cells
    #
    # The differentiated cell type is retained,
    # not the parental iPSC.
    #
    # Examples:
    #   iPSC-derived endothelial cells
    #       -> Endothelial cell
    #
    #   induced pluripotent stem cell-derived
    #   mesenchymal stem cells
    #       -> Mesenchymal stem cell
    # ------------------------------------------------------
    
    elif (
        "iPSC" in cell_hits
        and len(cell_hits) > 1
        and re.search(
            r"\b(?:ipscs?|induced pluripotent stem cells?)"
            r"[- ]derived\b",
            text,
            flags=re.IGNORECASE
        )
    ):
        cell_hits.remove("iPSC")
    
    
    # ------------------------------------------------------
    # Fibroblast-derived iPSC
    #
    # The final cell type is iPSC.
    #
    # Example:
    #   fetal fibroblast-derived iPSC line
    #       -> iPSC
    # ------------------------------------------------------
    
    if (
        "Fibroblast" in cell_hits
        and "iPSC" in cell_hits
        and re.search(
            r"\bfibroblasts?[- ]derived\b.*\b"
            r"(?:ipscs?|induced pluripotent stem cells?)\b",
            text,
            flags=re.IGNORECASE
        )
    ):
        cell_hits.remove("Fibroblast")
    
    
    # ------------------------------------------------------
    # Fibroblast-derived induced neurons
    #
    # The final cell type is Neuron.
    #
    # Example:
    #   fibroblast-derived induced neurons
    #       -> Neuron
    # ------------------------------------------------------
    
    if (
        "Fibroblast" in cell_hits
        and "Neuron" in cell_hits
        and re.search(
            r"\bfibroblasts?[- ]derived\s+induced\s+neurons?\b",
            text,
            flags=re.IGNORECASE
        )
    ):
        cell_hits.remove("Fibroblast")

    ########################################################
    # Resolve negative expressions
    ########################################################

    # cardiomyocytes and non-myocytes
    # -> Cardiomyocyte
    if (
        "Myocyte" in cell_hits
        and re.search(
            r"\bnon[- ]myocytes?\b",
            text,
            flags=re.IGNORECASE
        )
    ):
        cell_hits.remove(
            "Myocyte"
        )

    if (
        "Neural crest cell" in cell_hits
        and re.search(
            r"\bnon[- ]neural[ -]crest\b"
            r"|\bneural[ -]crest[ -]negative\b",
            text,
            flags=re.IGNORECASE
        )
    ):
        cell_hits.remove("Neural crest cell")

    ########################################################
    # Cell line only
    ########################################################

    if (
        len(cell_hits) == 0
        and len(cell_line_hits) > 0
    ):

        return (
            "Unknown",
            "",
            ";".join(
                cell_line_hits
            ),
            "Cell_line_only"
        )


    ########################################################
    # No recognized cell type / cell line
    ########################################################

    if len(cell_hits) == 0:

        return (
            "Other",
            "",
            "",
            "No_match"
        )


    ########################################################
    # Exactly one biological cell type
    ########################################################

    if len(cell_hits) == 1:

        if cell_line_hits:

            return (
                cell_hits[0],
                cell_hits[0],
                ";".join(
                    cell_line_hits
                ),
                "Single_match_with_cell_line"
            )

        return (
            cell_hits[0],
            cell_hits[0],
            "",
            "Single_match"
        )


    ########################################################
    # Multiple biological cell types
    ########################################################

    return (
        "Mixed",
        ";".join(
            cell_hits
        ),
        ";".join(
            cell_line_hits
        ),
        "Multiple_cell_types"
    )


############################################################
# 8. Decide whether a classification is usable
############################################################

def has_cell_information(result):

    cell_type_normalized = result[0]
    cell_line_name = result[2]
    status = result[3]

    return (
        status not in {
            "Missing",
            "No_match"
        }
        or bool(
            cell_line_name
        )
    )


############################################################
# 9. Cell type / Cell line normalization with rescue
#
# Priority:
#
# 1. Explicit Cell_type
#
# 2. Tissue
#    ONLY when Tissue_status == No_match
#
# 3. Source_name
#
# 4. Isolation_source
#
# Explicit Cell_type always has priority.
############################################################

def normalize_celltype_with_rescue(row):

    ########################################################
    # 1. Explicit Cell_type
    ########################################################

    cell_type_raw = str(
        row.get(
            "Cell_type",
            ""
        )
    ).strip()


    primary = normalize_celltype(
        cell_type_raw
    )


    if has_cell_information(
        primary
    ):

        return (
            primary[0],
            primary[1],
            primary[2],
            primary[3],
            "Cell_type"
        )


    ########################################################
    # 2. Cross-rescue from Tissue
    #
    # Only Tissue values that FAILED anatomical tissue
    # normalization are tested here.
    ########################################################

    tissue_status = str(
        row.get(
            "Tissue_status",
            ""
        )
    ).strip()


    tissue_raw = str(
        row.get(
            "Tissue",
            ""
        )
    ).strip()


    if (
        tissue_status == "No_match"
        and not is_missing_raw(
            tissue_raw
        )
    ):

        tissue_result = normalize_celltype(
            tissue_raw
        )


        if has_cell_information(
            tissue_result
        ):

            return (
                tissue_result[0],
                tissue_result[1],
                tissue_result[2],
                tissue_result[3],
                "Tissue"
            )


    ########################################################
    # 3. Rescue from Source_name
    ########################################################

    source_name = str(
        row.get(
            "Source_name",
            ""
        )
    ).strip()


    if not is_missing_raw(
        source_name
    ):

        source_result = normalize_celltype(
            source_name
        )


        if has_cell_information(
            source_result
        ):

            return (
                source_result[0],
                source_result[1],
                source_result[2],
                source_result[3],
                "Source_name"
            )


    ########################################################
    # 4. Rescue from Isolation_source
    ########################################################

    isolation_source = str(
        row.get(
            "Isolation_source",
            ""
        )
    ).strip()


    if not is_missing_raw(
        isolation_source
    ):

        isolation_result = normalize_celltype(
            isolation_source
        )


        if has_cell_information(
            isolation_result
        ):

            return (
                isolation_result[0],
                isolation_result[1],
                isolation_result[2],
                isolation_result[3],
                "Isolation_source"
            )

    ############################################################
    # Rescue known cell line from BioSample cell-line attribute
    ############################################################
    
    all_char = str(
        row.get(
            "All_characteristics",
            ""
        )
    ).strip()
    
    cell_line_values = extract_attribute_values(
        all_char,
        {
            "cell line",
            "cell lines",
        }
    )
    
    for value in cell_line_values:
    
        cell_line_hits = find_hits(
            value,
            COMPILED_CELL_LINE_PATTERNS
        )
    
        if cell_line_hits:
    
            return (
                "Unknown",
                "",
                ";".join(cell_line_hits),
                "Cell_line_only",
                "Cell_line_attribute"
            )


    ########################################################
    # 5. Nothing rescued
    ########################################################

    return (
        primary[0],
        primary[1],
        primary[2],
        primary[3],
        "Cell_type"
    )


############################################################
# 10. Read metadata
############################################################

print(
    "Reading:",
    INPUT_FILE
)

df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)

print(
    "Number of BioSamples:",
    len(df)
)


############################################################
# 11. Check required columns
############################################################

REQUIRED_COLUMNS = {
    "Cell_type",
    "Tissue",
    "Tissue_status",
    "Source_name",
    "Isolation_source",
}

missing_columns = (
    REQUIRED_COLUMNS
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


############################################################
# 12. Normalize Cell type / Cell line
############################################################

print(
    "Normalizing Cell type / Cell line..."
)

results = df.apply(
    normalize_celltype_with_rescue,
    axis=1
)


df["Cell_type_normalized"] = [
    x[0]
    for x in results
]

df["Cell_type_match"] = [
    x[1]
    for x in results
]

df["Cell_line_name"] = [
    x[2]
    for x in results
]

df["Cell_type_status"] = [
    x[3]
    for x in results
]

df["Cell_type_normalization_source"] = [
    x[4]
    for x in results
]


############################################################
# 13. Save complete table
############################################################

df.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)


############################################################
# 14. Summary
############################################################

print()
print(
    "=== Cell type normalized ==="
)

print(
    df[
        "Cell_type_normalized"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)


print()
print(
    "=== Cell type status ==="
)

print(
    df[
        "Cell_type_status"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)


print()
print(
    "=== Cell type normalization source ==="
)

print(
    df[
        "Cell_type_normalization_source"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)


print()
print(
    "=== Cell lines ==="
)

print(
    df.loc[
        df["Cell_line_name"] != "",
        "Cell_line_name"
    ]
    .value_counts()
    .head(50)
    .to_string()
)


############################################################
# 15. Save unclassified values
#
# Important:
# Use the ACTUAL source text that ultimately failed.
#
# For now, retain the explicit Cell_type value for review.
############################################################

unclassified = (
    df.loc[
        df["Cell_type_status"]
        == "No_match",
        "Cell_type"
    ]
    .value_counts()
    .rename_axis(
        "Cell_type_original"
    )
    .reset_index(
        name="Count"
    )
)

unclassified.to_csv(
    UNCLASSIFIED_FILE,
    sep="\t",
    index=False
)


############################################################
# 16. Save multiple-cell-type values
############################################################

multiple = (
    df.loc[
        df["Cell_type_status"]
        == "Multiple_cell_types",
        [
            "Cell_type",
            "Tissue",
            "Source_name",
            "Isolation_source",
            "Cell_type_match",
            "Cell_line_name",
            "Cell_type_normalization_source",
        ]
    ]
    .value_counts()
    .reset_index(
        name="Count"
    )
)

multiple.to_csv(
    MULTIPLE_FILE,
    sep="\t",
    index=False
)


############################################################
# 17. Save cell-line counts
############################################################

cellline_counts = (
    df.loc[
        df["Cell_line_name"] != "",
        "Cell_line_name"
    ]
    .value_counts()
    .rename_axis(
        "Cell_line_name"
    )
    .reset_index(
        name="Count"
    )
)

cellline_counts.to_csv(
    CELLLINE_FILE,
    sep="\t",
    index=False
)


############################################################
# 18. Save rescued samples
############################################################

rescued = df.loc[
    df[
        "Cell_type_normalization_source"
    ].isin(
        [
            "Tissue",
            "Source_name",
            "Isolation_source"
        ]
    ),
    [
        "BioSample",
        "Tissue",
        "Cell_type",
        "Tissue_status",
        "Source_name",
        "Isolation_source",
        "Cell_type_normalized",
        "Cell_type_match",
        "Cell_line_name",
        "Cell_type_status",
        "Cell_type_normalization_source",
    ]
]


rescued.to_csv(
    RESCUED_FILE,
    sep="\t",
    index=False
)


############################################################
# 19. Finish
############################################################

print()
print("Saved:")
print(OUTPUT_FILE)
print(UNCLASSIFIED_FILE)
print(MULTIPLE_FILE)
print(CELLLINE_FILE)
print(RESCUED_FILE)

print("\n=== Cell type normalization source ===")

print(
    df["Cell_type_normalization_source"]
    .value_counts()
)


print("\n=== Cell type status ===")

print(
    df["Cell_type_status"]
    .value_counts()
)


print("\n=== Cell lines ===")

print(
    df.loc[
        df["Cell_line_name"] != "",
        "Cell_line_name"
    ]
    .value_counts()
    .head(50)
)

rescued = df[
    df["Cell_type_normalization_source"].isin(
        [
            "Tissue",
            "Source_name",
            "Isolation_source"
        ]
    )
]

rescued[
    [
        "Tissue",
        "Cell_type",
        "Source_name",
        "Isolation_source",
        "Cell_type_normalized",
        "Cell_type_match",
        "Cell_line_name",
        "Cell_type_status",
        "Cell_type_normalization_source"
    ]
].to_csv(
    "celltype_rescued.tsv",
    sep="\t",
    index=False
)
