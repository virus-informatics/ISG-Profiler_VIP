#!/usr/bin/env python3

import re
import pandas as pd


############################################################
# 1. Input / output
############################################################

INPUT_FILE = "biosample_metadata.tsv"

OUTPUT_FILE = "biosample_metadata_tissue_normalized_v2.tsv"

UNCLASSIFIED_FILE = "tissue_unclassified_counts_v2.tsv"

MIXED_FILE = "tissue_mixed_counts_v2.tsv"


############################################################
# 2. Tissue normalization dictionary
############################################################

TISSUE_PATTERNS = {

    "Blood": [
        r"\bwhole blood\b",
        r"\bperipheral blood\b",
        r"\bblood\b",
        r"\bpbmcs?\b",
        r"\bperipheral blood mononuclear",
        r"\bwhite blood cells?\b",
        r"\bleukocytes?\b",
    ],

    "Liver": [
        r"\bliver\b",
        r"\bhepatic\b",
    ],

    # IMPORTANT:
    # Do NOT use generic "cortex"
    # because kidney cortex / adrenal cortex / ovarian cortex exist.
    "Brain/CNS": [
        r"\bbrain\b",
        r"\bcerebral\b",
        r"\bcerebrum\b",
        r"\bcerebral cortex\b",
        r"\bprefrontal cortex\b",
        r"\bfrontal cortex\b",
        r"\bmotor cortex\b",
        r"\bvisual cortex\b",
        r"\bauditory cortex\b",
        r"\banterior cingulate cortex\b",

        r"\bhippocamp",
        r"\bhipppocampus\b",  # typo in metadata
        r"\bbrian tissue\b",  # typo
        r"\bhypothamalus\b", #typo
        r"\bdentate gyrus\b",
        r"\bcornu ammonis\b",

        r"\bhypothalam",
        r"\bpre[ -]?optic area\b",
        r"\bpreoptic area\b",
        r"\bmedial preoptic area\b",

        r"\bcerebell",
        r"\bamygdal",

        # striat だと "striated muscle" にhitするので限定
        r"\bstriatum\b",
        r"\bstriatal\b",

        r"\bnucleus[ _-]?accumbens\b",
        r"\bcaudate nucleus\b",
        r"\bcaudate\b",
        r"\bputamen\b",
        r"\blateral septal nucleus\b",
        r"\bnucleus of solitary tract\b",

        r"\bthalamus\b",
        r"\bforebrain\b",
        r"\btelencephalon\b",
        r"\bmedulla oblongata\b",
        r"\barea postrema\b",
        r"\bsuperior temporal sulcus\b",
        r"\barea 46\b",

        r"\bspinal cord\b",
        r"\bolfactory bulb\b",

        #add manually
        r"\bspinal\b",
        r"\bmidbrain\b",
        r"\barcuate nucleus\b",
        r"\blateral geniculate nucleus\b",
        r"\bammon'?s horn\b",
        r"\bhyperpallium\b",
        r"\bmesopallium\b",
    ],

    "Lung/Respiratory": [
        r"\blung\b",
        r"\blungs\b",
        # r"\bpulmonary\b",
        r"\bbronchoalveolar[ _-]+lavage\b",
    ],

    "Kidney": [
        r"\bkidney\b",
        r"\bkidneys\b",
        r"\brenal\b",
        # typo
        r"\bkindey\b",
    ],

    "Heart": [
        r"\bheart\b",
        r"\bcardiac\b",
        r"\bmyocard",

        # generic ventricle is dangerous:
        # brain also has ventricles
        r"\bleft ventricle\b",
        r"\bright ventricle\b",
        r"\bcardiac ventricle\b",
        r"\bventricular myocard",
        r"\bsinoatrial\b",
        r"\bright atrium\b",
        r"\bleft atrium\b",
        r"\bmitral valve\b",
    ],

    "Spleen": [
        r"\bspleen\b",
        r"\bsplenic\b",
    ],

    "Muscle": [
        r"\bmuscle\b",
        r"\bskeletal muscle\b",
        r"\bstriated muscle\b",
        r"\blongissimus\b",
        r"\bgastrocnemius\b",
        r"\btibialis\b",
        r"\bsoleus\b",
        r"\bgluteus\b",
        r"\bpectoralis\b",
        r"\bbiceps femoris\b",
        r"\bsemitendinosus\b",

        # added from unclassified values
        r"\bdiaphragm\b",
        r"\blonguissimus_dorsi\b",
        #add manually
        r"\bplantaris\b",
        r"\bchicken_skeletal_muscle\b",
        ],

    "Skin": [
        r"\bskin\b",
        r"\bdermis\b",
        r"\bepiderm",
        r"\bfeather\b",
        r"\bhair follicle\b",
        r"\bwing membrane\b",
        r"\bpatagium\b",
        r"^hair$",
    ],

    "Ovary/Oocyte": [
        r"\bovary\b",
        r"\bovaries\b",
        r"\bovarian\b",
        r"\boocyte",
        # generic follicle removed
        # r"\bfollicle\b",
        r"\bovarian follicle\b",
        r"\bovary follicle\b",
        r"\bpreovulatory follicle\b",
    ],

    "Testis": [
        r"\btestis\b",
        r"\btestes\b",
        r"\btesticular\b",
        r"\bepididym",
        r"\bsperm\b",
    ],
    "Prostate": [
        r"\bprostate\b",
        r"\bprostatic\b",
    ],

    # Samples simply described as "gonad"
    "Gonad": [
        r"\bgonad\b",
        r"\bgonadal\b",
    ],

    "Intestine": [
        r"\bintestin",
        r"\bileum\b",
        r"\bileal\b",
        r"\bjejunum\b",
        r"\bduodenum\b",
        r"\bcolon\b",
        r"\bcecum\b",
        r"\bcaecum\b",
        r"\brectum\b",
        r"\bgut\b",
        r"\bpeyer'?s? patches?\b",
        r"\bpeyer'?s? patch\b",
    ],

    "Adipose": [
        r"\badipose\b",
        r"\babdominal fat\b",
        r"\bsubcutaneous fat\b",
        r"\bvisceral fat\b",
        r"\bmesenteric fat\b",
        r"\bbackfat\b",
        r"\bfat tissue\b",
        r"\bfat pad\b",
        r"\bomentum\b",
        # exact "fat" is dealt with below as well
        r"^fat$",
# add manually
        r"\btail fat\b",
        r"\bperirenal fat\b",
        r"\bpectoral fatty tissue\b",
    ],

    "Mammary": [
        r"\bmammary\b",
        r"\bbreast tissue\b",
        r"^breast$",
        r"\bmilk\b",
        r"\bcolostrm\b",
        r"^Colostrum$",
    ],

    "Uterus/Reproductive tract": [
        r"\buterus\b",
        r"\buterine\b",
        r"\bendomet",
        r"\bcervix\b",
        r"\boviduct\b",
        r"\bcaruncular\b",
        r"\bintercaruncular\b",
        r"\bmyometrium\b",
        r"\bshell gland\b",
    ],

    "Pancreas": [
        r"\bpancrea",
    ],

    "Thymus": [
        r"\bthymus\b",
        r"\bthymic\b",
    ],

    "Bone marrow": [
        r"\bbone marrow\b",
    ],

    "Placenta": [
        r"\bplacenta",
        r"\bdecidua\b",
        r"\bchorioallanto",
    ],

    "Eye/Retina": [
        r"\bretina",
        r"\bretinal\b",
        r"\beye\b",
        r"\brpe\b",
        r"\bchoroid\b",
        r"\biris\b",
        r"\bconjunctiv",
    ],

    "Lymph node": [
        r"\blymph[ _-]?nodes?\b",
    
        # abbreviations / specific descriptions
        r"\bld[-_ ]?ln\b",
        r"\bpln\b",
        r"\bmesln\b",
        r"\bretropharyngeal[ _-]?nodes?\b",
    ],

    "Tonsil": [
        r"\btonsil\b",
        r"\btonsillar\b",
    ],

    "Bursa of Fabricius": [
        r"\bbursa of fabricius\b",
        r"\bbursal\b",
    ],

    "Pituitary": [
        r"\bpituitary\b",
        r"\bpars tuberalis\b",
    ],

    "Pineal gland": [
        r"\bpineal gland\b",
        r"^pineal$",
    ],

    "Adrenal": [
        r"\badrenal\b",
    ],

    "Thyroid": [
        r"\bthyroid\b",
    ],

    "Stomach/Foregut": [
        r"\bstomach\b",
        r"\brumen\b",
        r"\babomas",
        r"\bomasum\b",
        r"\bgizzard\b",
        r"\bproventriculus\b",
    ],

    "Trachea": [
        r"\btrachea\b",
        r"\btracheal\b",
    ],

    "Nerve/PNS": [
        r"\bsciatic nerve\b",
        r"\bperipheral nerve\b",
        r"\bdorsal root gangli",
        r"\bdrg\b",
    ],

    "Vascular": [
        r"\baorta\b",
        r"\baortic\b",
        r"\bartery\b",
        r"\barterial\b",
        r"\bvein\b",
        r"\bvascular tissue\b",
    ],

    "Bladder": [
        r"\bbladder\b",
    ],

    "Synovium": [
        r"\bsynovium\b",
        r"\bsynovial\b",
    ],

    "Ear/Cochlea": [
        r"\bcochlea\b",
        r"\bcochlear\b",
        r"\binner ear\b",
        r"^ear$",
        r"^Ear$",
    ],

    "Bone/Cartilage": [
        # Don't use generic bone marrow as Bone/Cartilage.
        r"\bcortical bone\b",
        r"\btrabecular bone\b",
        r"\bbone tissue\b",
        r"\bcartilage\b",
        r"\bcartilaginous\b",
        r"\bbone\b",
        r"\btibia\b",
        r"\bmeniscus\b",
        r"\btendon\b",
        r"\bantler\b",
    ],
    
    "Esophagus": [
        r"\besophagus\b",
        r"\besophageal\b",
        r"\boesophagus\b",
        r"\boesophageal\b",
    ],

    "Harderian gland": [
        r"\bharderian gland\b",
    ],
    
    "Amniotic fluid": [
        r"\bamniotic fluid\b",
    ],
    
    "Nasal/Pharynx": [
        r"\bnasal pharynx\b",
        r"\bnasopharynx\b",
        r"\bvomeronasal\b",
        r"\bolfactory epitheli(?:um|a)\b",
        r"\bolfactory mucosa\b",
    ],

    "Oral": [
        r"\blip\b",
        r"\blips\b",
        r"\bsalivary gland\b",
        r"\bsalivary\b",
        r"\bgingiva\b",
        r"\bgingival\b",
        r"\boral\b",
        r"\bbuccal\b",
        r"\btongue\b",
        r"\bpalate\b",
]
}


############################################################
# 3. Non-tissue sample context
#
# These should NOT automatically generate "Mixed".
############################################################

CONTEXT_PATTERNS = {

    "Organoid": [
        r"\borganoid\b",
        r"\borganoids\b",
    ],

    "Tumor/Neoplasm": [
        r"\btumou?r\b",
        r"\bneoplasm\b",
        r"\bcarcinoma\b",
        r"\bmelanoma\b",
        r"\bsarcoma\b",
    ],

    "Embryonic/Developmental": [
        r"\bembryos?\b",
        r"\bembryonic\b",
        r"\bmorula\b",
        r"\bblastocyst\b",
        r"\binner cell mass\b",

        # neural crestもここ
        r"\bneural crest\b",
    ],
}

############################################################
# 4. Values treated as missing and mixed patterns
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

MIXED_PATTERNS = [
    r"^mixed$",
    r"^organs? mix$",
    r"^mixture of multiple tissues$",
    r"^mixed tissues?$",
    r"^mixed[- ]tissues?$",
    r"^mixed pool$",
    r"^mixed \d+ tissues$",
]

############################################################
# 5. Compile regex once
############################################################

COMPILED_TISSUE_PATTERNS = {
    category: [
        re.compile(pattern, flags=re.IGNORECASE)
        for pattern in patterns
    ]
    for category, patterns in TISSUE_PATTERNS.items()
}


COMPILED_CONTEXT_PATTERNS = {
    category: [
        re.compile(pattern, flags=re.IGNORECASE)
        for pattern in patterns
    ]
    for category, patterns in CONTEXT_PATTERNS.items()
}

COMPILED_MIXED_PATTERNS = [
    re.compile(pattern, flags=re.IGNORECASE)
    for pattern in MIXED_PATTERNS
]

############################################################
# 6. Helper
############################################################

def find_hits(text, pattern_dict):

    hits = []

    for category, patterns in pattern_dict.items():

        if any(pattern.search(text) for pattern in patterns):
            hits.append(category)

    return list(dict.fromkeys(hits))


############################################################
# 7. Tissue classification function
############################################################

def normalize_tissue(value):

    ########################################################
    # Missing
    ########################################################

    if pd.isna(value):
        return (
            "Unknown",
            "",
            "",
            "Missing"
        )

    text = str(value).strip()

    if text.lower() in MISSING_VALUES:
        return (
            "Unknown",
            "",
            "",
            "Missing"
        )


    ########################################################
    # Explicit mixed-tissue samples
    ########################################################

    if any(
        pattern.search(text)
        for pattern in COMPILED_MIXED_PATTERNS
    ):
        return (
            "Mixed",
            "Explicit mixed tissue",
            "",
            "Multiple_tissues"
        )

    ########################################################
    # Detect tissue and sample context separately
    ########################################################

    tissue_hits = find_hits(
        text,
        COMPILED_TISSUE_PATTERNS
    )

    context_hits = find_hits(
        text,
        COMPILED_CONTEXT_PATTERNS
    )


    ########################################################
    # Resolve known overlaps
    ########################################################

    # Bone marrow should not become Bone/Cartilage
    if (
        "Bone marrow" in tissue_hits
        and "Bone/Cartilage" in tissue_hits
    ):
        tissue_hits.remove(
            "Bone/Cartilage"
        )


    # More specific reproductive organs beat generic Gonad
    if "Gonad" in tissue_hits:

        if (
            "Ovary/Oocyte" in tissue_hits
            or "Testis" in tissue_hits
        ):
            tissue_hits.remove("Gonad")


    # Cardiac / heart muscle is Heart, not Mixed Heart + Muscle
    if (
        "Heart" in tissue_hits
        and "Muscle" in tissue_hits
    ):
    
        if re.search(
            r"\bcardiac muscle\b"
            r"|\bheart muscle\b"
            r"|\bmyocard",
            text,
            flags=re.IGNORECASE
        ):
            tissue_hits.remove("Muscle")

    if re.search(
        r"\b(?:hepatic|ileum|ileal|jejunum|jejunal|intestinal|abomasal)"
        r"[ _-]+lymph[ _-]?nodes?\b",
        text,
        flags=re.IGNORECASE
    ):
        tissue_hits = ["Lymph node"]

    # Epididymal adipose tissue is adipose, not Testis
    if (
        "Testis" in tissue_hits
        and "Adipose" in tissue_hits
        and re.search(
            r"\bepididymal adipose\b",
            text,
            flags=re.IGNORECASE
        )
    ):
        tissue_hits.remove("Testis")

    # "tumor-free" means non-tumor tissue
    if (
        "Tumor/Neoplasm" in context_hits
        and re.search(
            r"\btumou?r[- ]free\b",
            text,
            flags=re.IGNORECASE
        )
    ):
        context_hits.remove("Tumor/Neoplasm")


    ########################################################
    # No anatomical tissue detected
    ########################################################

    if len(tissue_hits) == 0:

        # e.g. "tumor" alone or "organoid" alone
        if context_hits:
            return (
                "Other",
                "",
                ";".join(context_hits),
                "Context_only"
            )

        return (
            "Other",
            "",
            "",
            "No_match"
        )


    ########################################################
    # Exactly one anatomical tissue detected
    ########################################################

    if len(tissue_hits) == 1:

        status = "Single_match"

        if context_hits:
            status = "Single_match_with_context"

        return (
            tissue_hits[0],
            tissue_hits[0],
            ";".join(context_hits),
            status
        )


    ########################################################
    # Multiple anatomical tissues detected
    #
    # Only anatomical tissue hits create Mixed.
    ########################################################

    return (
        "Mixed",
        ";".join(tissue_hits),
        ";".join(context_hits),
        "Multiple_tissues"
    )


############################################################
# 8. Read metadata
############################################################

print("Reading:", INPUT_FILE)

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
# 9. Check required column
############################################################

if "Tissue" not in df.columns:

    raise ValueError(
        "Input file does not contain a 'Tissue' column."
    )


############################################################
# 10. Normalize Tissue
#
# Priority:
# 1. Explicit Tissue field
# 2. Source_name rescue
# 3. Isolation_source rescue
############################################################

# check whether raw metadata is missing 
def is_missing_raw(value):

    text = str(value).strip().lower()

    return text in MISSING_VALUES

def merge_contexts(*contexts):

    hits = []

    for context in contexts:

        if not context:
            continue

        for x in context.split(";"):

            x = x.strip()

            if (
                x
                and x not in hits
            ):
                hits.append(x)

    return ";".join(hits)

def normalize_tissue_with_rescue(row):

    ########################################################
    # Raw values
    ########################################################

    tissue_raw = str(
        row.get("Tissue", "")
    ).strip()

    cell_type_raw = str(
        row.get("Cell_type", "")
    ).strip()


    ########################################################
    # 1. Normalize explicit Tissue field
    ########################################################

    primary = normalize_tissue(
        tissue_raw
    )

    primary_tissue = primary[0]
    primary_match = primary[1]
    primary_context = primary[2]
    primary_status = primary[3]


    ########################################################
    # Explicit Tissue produced a usable anatomical result
    #
    # Keep it.
    ########################################################

    if primary_status in {
        "Single_match",
        "Single_match_with_context",
        "Multiple_tissues"
    }:

        return (
            primary_tissue,
            primary_match,
            primary_context,
            primary_status,
            "Tissue"
        )


    ########################################################
    # Auxiliary Tissue rescue is allowed only when
    # Cell_type itself is missing.
    ########################################################

    cell_type_missing = is_missing_raw(
        cell_type_raw
    )


    if not cell_type_missing:

        return (
            primary_tissue,
            primary_match,
            primary_context,
            primary_status,
            "Tissue"
        )


    ########################################################
    # 2. Try Source_name
    ########################################################

    source_name = str(
        row.get(
            "Source_name",
            ""
        )
    ).strip()


    source_result = normalize_tissue(
        source_name
    )

    source_tissue = source_result[0]
    source_match = source_result[1]
    source_context = source_result[2]
    source_status = source_result[3]


    ########################################################
    # Source_name produced anatomical Tissue
    ########################################################

    if source_status in {
        "Single_match",
        "Single_match_with_context",
        "Multiple_tissues"
    }:

        combined_context = merge_contexts(
            primary_context,
            source_context
        )

        final_status = source_status

        if (
            source_status == "Single_match"
            and combined_context
        ):
            final_status = (
                "Single_match_with_context"
            )

        return (
            source_tissue,
            source_match,
            combined_context,
            final_status,
            "Source_name"
        )


    ########################################################
    # 3. Try Isolation_source
    ########################################################

    isolation_source = str(
        row.get(
            "Isolation_source",
            ""
        )
    ).strip()


    isolation_result = normalize_tissue(
        isolation_source
    )

    isolation_tissue = isolation_result[0]
    isolation_match = isolation_result[1]
    isolation_context = isolation_result[2]
    isolation_status = isolation_result[3]


    ########################################################
    # Isolation_source produced anatomical Tissue
    ########################################################

    if isolation_status in {
        "Single_match",
        "Single_match_with_context",
        "Multiple_tissues"
    }:

        combined_context = merge_contexts(
            primary_context,
            source_context,
            isolation_context
        )

        final_status = isolation_status

        if (
            isolation_status == "Single_match"
            and combined_context
        ):
            final_status = (
                "Single_match_with_context"
            )

        return (
            isolation_tissue,
            isolation_match,
            combined_context,
            final_status,
            "Isolation_source"
        )


    ########################################################
    # 4. No anatomical Tissue was rescued
    #
    # Context information may still be useful.
    ########################################################

    combined_context = merge_contexts(
        primary_context,
        source_context,
        isolation_context
    )


    if combined_context:

        return (
            "Other",
            "",
            combined_context,
            "Context_only",
            (
                "Tissue"
                if primary_context
                else
                "Source_name"
                if source_context
                else
                "Isolation_source"
            )
        )



    ########################################################
    # Nothing useful
    ########################################################

    return (
        primary_tissue,
        primary_match,
        primary_context,
        primary_status,
        "Tissue"
    )

############################################################
# 10. Normalize Tissue
############################################################

print("Normalizing Tissue...")

results = df.apply(
    normalize_tissue_with_rescue,
    axis=1
)

df["Tissue_normalized"] = [
    x[0] for x in results
]

df["Tissue_match"] = [
    x[1] for x in results
]

df["Sample_context"] = [
    x[2] for x in results
]

df["Tissue_status"] = [
    x[3] for x in results
]

df["Tissue_normalization_source"] = [
    x[4] for x in results
]

############################################################
# 11. Save complete table
############################################################

df.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)


############################################################
# 12. Summary
############################################################

print()
print("=== Tissue normalized ===")

print(
    df["Tissue_normalized"]
    .value_counts(dropna=False)
    .to_string()
)


print()
print("=== Tissue status ===")

print(
    df["Tissue_status"]
    .value_counts(dropna=False)
    .to_string()
)


print()
print("=== Sample context ===")

print(
    df["Sample_context"]
    .value_counts(dropna=False)
    .head(20)
    .to_string()
)


############################################################
# 13. Save unclassified values
############################################################

unclassified = (
    df.loc[
        df["Tissue_status"].isin(
            [
                "No_match",
                "Context_only"
            ]
        ),
        [
            "Tissue",
            "Sample_context"
        ]
    ]
    .value_counts()
    .reset_index(name="Count")
)

unclassified.to_csv(
    UNCLASSIFIED_FILE,
    sep="\t",
    index=False
)


############################################################
# 14. Save true Mixed values
############################################################

mixed = (
    df.loc[
        df["Tissue_status"]
        == "Multiple_tissues",
        [
            "Tissue",
            "Tissue_match",
            "Sample_context"
        ]
    ]
    .value_counts()
    .reset_index(name="Count")
)

mixed.to_csv(
    MIXED_FILE,
    sep="\t",
    index=False
)


############################################################
# 15. Finish
############################################################

print()
print("Saved:")
print(OUTPUT_FILE)
print(UNCLASSIFIED_FILE)
print(MIXED_FILE)

