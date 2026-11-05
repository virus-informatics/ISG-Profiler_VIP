#!/usr/bin/env python3

import re
import pandas as pd


############################################################
# 1. Input / output
############################################################

INPUT_FILE = "biosample_metadata_celltype_normalized_v3.tsv"
OUTPUT_FILE = "biosample_metadata_vivo_vitro_v3.tsv"


############################################################
# 2. Regex patterns
############################################################

IN_VITRO_PATTERNS = [
    r"\bin[ _-]?vitro\b",
    r"\bcell[ _-]?cultures?\b",
    r"\bcultured[ _-]?cells?\b",
    r"\bcells?[ _-]?cultured\b",
    r"\bcultured[ _-]?in\b",
    r"\bprimary[ _-]?cultures?\b",
    r"\bprimary[ _-]?cell[ _-]?cultures?\b",
    r"\bcell[ _-]?monolayers?\b",
    # Require the word "passage"; do not use a generic P0/P28 rule.
    r"\bpassage[ _-]?(?:number[ _-]?)?[=:]?[ _-]?[pP]?\d+\b",
    r"\bimmortali[sz]ed\b",
    r"\bipsc[- ]derived\b",
    r"\binduced pluripotent stem cell[- ]derived\b",
]

# Additional strong culture phrases identified during QC of samples
# that had originally been classified as in vivo.
# These are deliberately specific; generic "culture" or "primary cell"
# alone is NOT treated as sufficient evidence.
STRONG_CULTURE_PATTERNS = [
    r"\b\d+\s*(?:days?|hours?|hrs?)\s+in\s+culture(?:\s+(?:medium|media))?\b",
    r"\bcultured\s+for\s+\d+\s*(?:days?|hours?|hrs?)?\b",
    r"\bculture\s+time\s*[:=]?\s*\d+\b",
    r"\bcells?\s+(?:were\s+)?cultured\b",
]

EX_VIVO_PATTERNS = [
    r"\bex[ _-]vivo\b",
    r"\bexplants?\b",
    r"\btissue explants?\b",
]

# Stem-cell / culture-type cell annotations that strongly imply in vitro
IN_VITRO_CELL_TYPES = {
    "iPSC",
    "Embryonic stem cell",
}

# Sample contexts that strongly imply in vitro
IN_VITRO_CONTEXTS = {
    "Organoid",
}


############################################################
# 3. Compile regex once
############################################################

COMPILED_IN_VITRO = [
    re.compile(x, flags=re.IGNORECASE)
    for x in IN_VITRO_PATTERNS
]

COMPILED_EX_VIVO = [
    re.compile(x, flags=re.IGNORECASE)
    for x in EX_VIVO_PATTERNS
]

COMPILED_STRONG_CULTURE = [
    re.compile(x, flags=re.IGNORECASE)
    for x in STRONG_CULTURE_PATTERNS
]


############################################################
# 4. Helper
############################################################

def first_regex_match(text, patterns):

    for pattern in patterns:
        m = pattern.search(text)
        if m:
            return m.group(0)

    return ""

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

            if value:
                values.append(value)

    return values


def extract_all_attribute_values(text):
    """
    Return only BioSample attribute VALUES from All_characteristics.
    This avoids treating an attribute key such as "cell line=" as
    evidence by itself.
    """

    if not text:
        return []

    values = []

    for item in text.split(" | "):

        if "=" not in item:
            continue

        _, value = item.split("=", 1)
        value = value.strip()

        if value:
            values.append(value)

    return values


def build_metadata_value_text(row, all_char):
    """
    Build a search string from metadata VALUES only.

    All_characteristics contains key=value pairs.  We search only the
    values so that a malformed/irrelevant "cell line" attribute key
    does not automatically force an in vitro classification.
    """

    columns = [
        "Tissue",
        "Cell_type",
        "Tissue_source",
        "Cell_type_source",
        "Source_name",
        "Isolation_source",
    ]

    values = []

    for col in columns:
        value = row.get(col, "")

        if pd.isna(value):
            continue

        value = str(value).strip()

        if value:
            values.append(value)

    values.extend(
        extract_all_attribute_values(all_char)
    )

    return " | ".join(values)

############################################################
# 5. Classification function
############################################################

def classify_vivo_vitro(row):

    ########################################################
    # Prepare values
    ########################################################

    all_char = str(
        row.get("All_characteristics", "")
    ).strip()

    cell_line = str(
        row.get("Cell_line_name", "")
    ).strip()

    cell_type = str(
        row.get("Cell_type_normalized", "")
    ).strip()

    tissue = str(
        row.get("Tissue_normalized", "")
    ).strip()

    tissue_status = str(
        row.get("Tissue_status", "")
    ).strip()

    sample_context = str(
        row.get("Sample_context", "")
    ).strip()

    # Search metadata values broadly, but do not use All_characteristics
    # attribute KEYS themselves as evidence.
    metadata_value_text = build_metadata_value_text(
        row,
        all_char
    )


    ########################################################
    # 1. Explicit ex vivo evidence
    #
    # Check this before generic culture terminology.
    ########################################################

    ex_vivo_match = first_regex_match(
        metadata_value_text,
        COMPILED_EX_VIVO
    )

    if ex_vivo_match:
        return (
            "ex vivo",
            f"All_characteristics: {ex_vivo_match}",
            "Explicit_ex_vivo"
        )


    ########################################################
    # 2. Known cell line -> in vitro
    ########################################################

    if (
        cell_line
        and cell_line.lower()
        not in {"unknown", "other", "nan"}
    ):
        return (
            "in vitro",
            f"Cell_line_name: {cell_line}",
            "Cell_line"
        )

    cell_line_values = extract_attribute_values(
        all_char,
        {
            "cell line",
            "cell lines",
        }
    )
    
    for value in cell_line_values:
    
        if re.search(
            r"\bcell[ _-]?lines?\b",
            value,
            flags=re.IGNORECASE
        ):
            return (
                "in vitro",
                "Explicit_in_vitro",
                "cell_line_attribute_value"
            )

    ########################################################
    # 3. Organoid -> in vitro
    ########################################################

    if sample_context in IN_VITRO_CONTEXTS:
        return (
            "in vitro",
            f"Sample_context: {sample_context}",
            "Organoid"
        )


    ########################################################
    # 4. Strong cell-type evidence -> in vitro
    ########################################################

    if cell_type in IN_VITRO_CELL_TYPES:
        return (
            "in vitro",
            f"Cell_type_normalized: {cell_type}",
            "Cultured_cell_type"
        )


    ########################################################
    # 6. Explicit culture / in vitro terminology
    ########################################################

    in_vitro_match = first_regex_match(
        metadata_value_text,
        COMPILED_IN_VITRO
    )

    if in_vitro_match:
        return (
            "in vitro",
            f"Metadata: {in_vitro_match}",
            "Explicit_in_vitro"
        )


    ########################################################
    # 7. Additional strong culture phrases from QC
    #
    # Generic "culture", "culture medium", "primary cell", or
    # "monolayer" alone are deliberately NOT sufficient.
    ########################################################

    strong_culture_match = first_regex_match(
        metadata_value_text,
        COMPILED_STRONG_CULTURE
    )

    if strong_culture_match:
        return (
            "in vitro",
            f"Metadata: {strong_culture_match}",
            "Explicit_culture"
        )


    ########################################################
    # 8. Clear tissue and no culture evidence -> in vivo
    #
    # Only use samples with a single interpretable tissue.
    ########################################################

    if tissue_status in {
        "Single_match",
        "Single_match_with_context"
    }:

        if tissue not in {
            "",
            "Unknown",
            "Other",
            "Mixed"
        }:
            return (
                "in vivo",
                f"Tissue_normalized: {tissue}",
                "Tissue_without_invitro_evidence"
            )


    ########################################################
    # 9. Otherwise unknown
    ########################################################

    return (
        "Unknown",
        "",
        "Insufficient_evidence"
    )


############################################################
# 6. Read metadata
############################################################

print("Reading:", INPUT_FILE)

df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)

print("Number of BioSamples:", len(df))


############################################################
# 7. Classify
############################################################

print("Classifying in vivo / in vitro...")

results = df.apply(
    classify_vivo_vitro,
    axis=1
)

df["Vivo_vitro"] = [
    x[0] for x in results
]

df["Vivo_vitro_evidence"] = [
    x[1] for x in results
]

df["Vivo_vitro_status"] = [
    x[2] for x in results
]


############################################################
# 8. Save
############################################################

df.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)


############################################################
# 9. Summary
############################################################

print()
print("=== Vivo / vitro ===")

print(
    df["Vivo_vitro"]
    .value_counts(dropna=False)
    .to_string()
)

print()
print("=== Classification status ===")

print(
    df["Vivo_vitro_status"]
    .value_counts(dropna=False)
    .to_string()
)


############################################################
# 10. Save uncertain samples for inspection
############################################################

unknown = df.loc[
    df["Vivo_vitro"] == "Unknown",
    [
        "BioSample",
        "Tissue",
        "Cell_type",
        "Tissue_normalized",
        "Cell_type_normalized",
        "Cell_line_name",
        "Sample_context",
        "All_characteristics"
    ]
]

unknown.to_csv(
    "vivo_vitro_unknown.tsv",
    sep="\t",
    index=False
)


############################################################
# 11. Save classified evidence summary
############################################################

evidence_summary = (
    df[
        [
            "Vivo_vitro",
            "Vivo_vitro_status",
            "Vivo_vitro_evidence"
        ]
    ]
    .value_counts()
    .reset_index(name="Count")
)

evidence_summary.to_csv(
    "vivo_vitro_evidence_counts.tsv",
    sep="\t",
    index=False
)


print()
print("Saved:")
print(OUTPUT_FILE)
print("vivo_vitro_unknown.tsv")
print("vivo_vitro_evidence_counts.tsv")

