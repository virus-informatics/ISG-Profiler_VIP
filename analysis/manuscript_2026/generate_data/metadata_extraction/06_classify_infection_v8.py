#!/usr/bin/env python3

import re
import pandas as pd


############################################################
# 1. Input / output
############################################################

INPUT_FILE = "biosample_metadata_vivo_vitro_v3.tsv"

OUTPUT_FILE = "biosample_metadata_infection_v8.tsv"

AGENT_COUNTS_FILE = "infection_agent_raw_counts.tsv"
UNKNOWN_FILE = "infection_classification_unknown.tsv"
REVIEW_FILE = "infection_classification_review.tsv"


############################################################
# 2. Missing values
############################################################

MISSING_VALUES = {
    "",
    "unknown",
    "missing",
    "mssing",
    "not applicable",
    "not available",
    "not collected",
    "not provided",
    "not specified",
    "na",
    "n/a",
    "none",
    "null",
}


############################################################
# 3. Attribute keys indicating an experimental infection
#
# These are relatively strong indicators of an infection
# experiment.
#
# Important:
# The presence of one of these keys means that the sample
# belongs to an experimental infection/challenge context.
#
# It does NOT necessarily mean that the particular sample
# is infected; it may be an experimental control.
############################################################

STRONG_EXPERIMENTAL_KEYS = {

    # challenge
    "challenge",
    "week_post_challenge",
    "weeks_post_challenge",
    "day_post_challenge",
    "days_post_challenge",
    "time_post_challenge",

    # infection timing
    "days_post_infection",
    "day_post_infection",
    "time_post_infection",
    "weeks_post_infection",
    "week_post_infection",

    # pathogen-specific timing seen in the dataset
    "weeks_post_siv_infection",
    "days_post_siv_infection",

    # inoculation details
    "infectious_dose",
    "infection_route",
    "virus_dose",
    "multiplicity_of_infection",
    "moi",
}


############################################################
# 4. Explicit experimental-infection phrases
#
# These phrases can also provide strong evidence even when
# stored under general keys such as "treatment".
############################################################

EXPERIMENTAL_PATTERNS = [

    r"\bexperimentally[ _-]+infected\b",
    r"\bexperimental[ _-]+infection\b",

    r"\bchallenged[ _-]+with\b",
    r"\bchallenge[ _-]+infection\b",
    r"\bviral[ _-]+challenge\b",

    r"\binoculated[ _-]+with\b",
    r"\bvirus[ _-]+inoculated\b",
    # Experimental infection
    r"\b\d+\s+(?:hours?|days?|weeks?)\s+post[ _-]+inoculation\b",
    r"\bMOI\s*=\s*[0-9.]+",
    r"\bPFU\b",
    r"\bplaque[- ]forming units?\b",
    r"\bchallenged\b.*\b(?:dpi|hpi)\b",
]


############################################################
# 5. Explicit non-experimental infection phrases
#
# Only explicit evidence is used here.
#
# Absence of experimental evidence does NOT automatically
# mean Non-experimental.
############################################################

NON_EXPERIMENTAL_PATTERNS = [
    r"\bnaturally[ _-]+infected\b",
    r"\bnatural[ _-]+infection\b",
    r"\bfield[ _-]+infected\b",
    r"\bfield[ _-]+infection\b",
    r"\bclinical[ _-]+infection\b",
]

NON_EXPERIMENTAL_KEYS = {
    "infection",
    "infection_status",
    "infected",
    "disease_status",
    "health_state",
    "sample_type",
}

EXPERIMENTAL_FREETEXT_PATTERNS = {

    "experimentally_infected":
        r"\bexperimentally[ _-]+infect(?:ed|ion)\b",

    "infected_with":
        r"\binfected[ _-]+with\b",

    "challenged_with":
        r"\bchalleng(?:ed|e)[ _-]+with\b",

    "inoculated_with":
        r"\binoculat(?:ed|ion)[ _-]+with\b",

    "post_infection":
        r"\bpost[ _-]+infection\b",

    "days_post_infection":
        r"\bdays?[ _-]+post[ _-]+infection\b",

    "dpi":
        r"\b\d+[ _-]*dpi\b",

    "post_challenge":
        r"\bpost[ _-]+challenge\b",

    "days_post_challenge":
        r"\bdays?[ _-]+post[ _-]+challenge\b",

    "dpc":
        r"\b\d+[ _-]*dpc\b",

    "moi":
        r"\b(?:moi|multiplicity[ _-]+of[ _-]+infection)\b",

    "infectious_dose":
        r"\binfectious[ _-]+dose\b",

    "mock_infected":
        r"\bmock[ _-]+infected\b",

    "uninfected_control":
        r"\buninfected[ _-]+control\b",
}


############################################################
# 5b. Tier-2 experimental-design rescue
#
# Applied only when:
#   Experimental_infection == "Unknown"
#   Infection_status == "Infected"
#
# "Infected" alone is NOT sufficient. We require an explicit
# time-point/design cue in All_characteristics.
############################################################

TIME_KEYS = {
    "time",
    "time_point",
    "timepoint",
    "day",
    "days",
    "sampling_time",
    "sampling_day",
    "time_post_infection",
    "time_post_challenge",
    "day_post_infection",
    "days_post_infection",
    "day_post_challenge",
    "days_post_challenge",
    "dpi",
    "dpc",
}

TIME_PATTERNS = [
    r"\bday[ _-]*\d+\b",                         # Day 3, day3
    r"\bd[ _-]?\d+\b",                           # D4, D10
    r"\b\d+\s*d\b",                              # 2d
    r"\b\d+\s*h\b",                              # 24h
    r"\b\d+\s*(?:hr|hrs|hour|hours)\b",          # 48 hours
    r"\b\d+\s*(?:day|days)\b",                   # 3 days
    r"\btime[ _-]?point\b",
    r"\bpost[ _-]?(?:infection|challenge)\b",
    r"\b\d+\s*(?:dpi|dpc)\b",
]

COMPILED_TIME_PATTERNS = [
    re.compile(x, flags=re.IGNORECASE)
    for x in TIME_PATTERNS
]

def detect_tier2_infected_timepoint(attributes, all_characteristics):

    """
    Rescue an unresolved Infected sample only when an
    experimental time-point/design cue is also present.

    Infection_status == Infected is checked by the caller.
    """

    text = str(all_characteristics or "")

    # Search the complete All_characteristics text so that
    # titles/values such as TPT_INFECTED_D4_rep3 or
    # Chix 22_Infected_Day 3 can be detected.
    time_match = None

    for pattern in COMPILED_TIME_PATTERNS:
        m = pattern.search(text)
        if m:
            time_match = m.group(0)
            break

    # Also accept explicit time-related attribute keys.
    time_key_evidence = ""

    for key_norm, value, key_original in attributes:
        if key_norm in TIME_KEYS and not is_missing(value):
            time_key_evidence = f"{key_original}={value}"
            break

    if not time_match and not time_key_evidence:
        return ""

    evidence = time_key_evidence or time_match

    return f"Infected + experimental timepoint: {evidence}"


############################################################
# Strong experimental infection metadata discovered by QC
############################################################

STRONG_EXPERIMENTAL_KEY_PATTERNS = [

    # e.g.
    # weeks_pre(-)_or_post(+)_experimental_infection
    r"experimental[\s_-]*infection",

    # post infection / post-infection time metadata
    r"post[\s_-]*infection",

    # post challenge time metadata
    r"post[\s_-]*challenge",

    # common compact time-point keys
    r"^dpi$",
    r"^dpc$",
]


def is_strong_experimental_key(key_norm):

    if not key_norm:
        return False

    for pattern in STRONG_EXPERIMENTAL_KEY_PATTERNS:

        if re.search(
            pattern,
            key_norm,
            flags=re.IGNORECASE
        ):
            return True

    return False

############################################################
# 5c. Tier-2 experimental-control rescue
#
# Applied only when:
#   Experimental_infection == "Unknown"
#   Infection_status == "Control/Uninfected"
#
# Goal:
# Rescue samples that look like controls of an infection
# experiment, while NOT converting ordinary healthy /
# uninfected samples into experimental controls.
############################################################

CONTROL_DESIGN_STRONG_PATTERNS = [
    # Strong infection-experiment control language
    r"\bmock[ _-]?infected\b",
    r"\bmock[ _-]?infection\b",
    r"\bmock[ _-]?challeng(?:e|ed)\b",
    r"\buninfected[ _-]+control\b",
    r"\bnon[ _-]?infected[ _-]+control\b",
    r"\binfection[ _-]+control\b",
    r"\bchallenge[ _-]+control\b",
]

CONTROL_MOCK_PATTERN = re.compile(
    r"\bmock\b",
    flags=re.IGNORECASE
)

CONTROL_OR_UNINFECTED_PATTERN = re.compile(
    r"\b(?:control|uninfected|non[ _-]?infected)\b",
    flags=re.IGNORECASE
)

CONTROL_REPLICATE_PATTERNS = [
    r"\breplicate\s*\d*\b",
    r"\brep[ _-]?\d+\b",
    r"(?:^|[_-])rep[ _-]?\d+\b",
]

CONTROL_INFECTION_DESIGN_PATTERNS = [
    r"\binfect(?:ed|ion)\b",
    r"\bchallenge(?:d)?\b",
    r"\binoculat(?:ed|ion)\b",
    r"\bpost[ _-]?(?:infection|challenge|inoculation)\b",
    r"\b(?:dpi|dpc)\b",
    r"\b(?:moi|pfu|tcid50)\b",
    r"\binfectious[ _-]?dose\b",
    r"\binfection[ _-]?route\b",
]

COMPILED_CONTROL_STRONG = [
    re.compile(x, flags=re.IGNORECASE)
    for x in CONTROL_DESIGN_STRONG_PATTERNS
]

COMPILED_CONTROL_REPLICATE = [
    re.compile(x, flags=re.IGNORECASE)
    for x in CONTROL_REPLICATE_PATTERNS
]

COMPILED_CONTROL_INFECTION_DESIGN = [
    re.compile(x, flags=re.IGNORECASE)
    for x in CONTROL_INFECTION_DESIGN_PATTERNS
]


def _first_regex_match(text, patterns):

    for pattern in patterns:
        m = pattern.search(text)
        if m:
            return m.group(0)

    return ""


def detect_tier2_control_design(attributes, all_characteristics):

    """
    Detect evidence that a Control/Uninfected sample is a
    control arm of an infection experiment.

    Returns:
        (status, evidence)

    Empty strings mean: do not rescue.

    Rules:
      1. explicit infection-control language -> rescue
      2. mock + timepoint -> rescue
      3. mock alone (within Control/Uninfected samples) -> rescue
      4. replicate + timepoint -> rescue
      5. replicate alone (within Control/Uninfected samples) -> rescue
      6. timepoint + infection/challenge design evidence -> rescue

    Important:
      This function is called only for samples already classified
      as Infection_status == "Control/Uninfected".

      Healthy/uninfected/control alone is still NOT sufficient.
      SPF alone is also NOT sufficient.
    """

    text = str(all_characteristics or "")

    # Strong explicit infection-control wording.
    strong = _first_regex_match(
        text,
        COMPILED_CONTROL_STRONG
    )

    if strong:
        return (
            "Tier2_control_explicit",
            f"Explicit experimental-control phrase: {strong}"
        )

    # Detect mock / replicate independently of time information.
    mock_match = CONTROL_MOCK_PATTERN.search(text)

    replicate = _first_regex_match(
        text,
        COMPILED_CONTROL_REPLICATE
    )

    control_match = (
        CONTROL_OR_UNINFECTED_PATTERN.search(text)
    )

    # Reuse the same time-point logic used for infected rescue.
    time_match = _first_regex_match(
        text,
        COMPILED_TIME_PATTERNS
    )

    time_key_evidence = ""

    for key_norm, value, key_original in attributes:

        if (
            key_norm in TIME_KEYS
            and not is_missing(value)
        ):
            time_key_evidence = (
                f"{key_original}={value}"
            )
            break

    time_evidence = (
        time_key_evidence
        or time_match
    )

    # Stronger mock + timepoint category.
    if mock_match and time_evidence:
        return (
            "Tier2_control_mock_timepoint",
            (
                "Mock control + experimental timepoint: "
                f"{time_evidence}"
            )
        )

    # New in v7:
    # Control/Uninfected + mock is sufficient even without a timepoint.
    if mock_match:
        return (
            "Tier2_control_mock",
            "Control/Uninfected + mock"
        )

    # Since the caller already requires Control/Uninfected, this is
    # mostly a consistency check against the raw metadata text.
    if not control_match:
        return ("", "")

    # Stronger replicate + timepoint category.
    if replicate and time_evidence:
        return (
            "Tier2_control_timepoint_replicate",
            (
                "Control/Uninfected + timepoint + replicate: "
                f"{time_evidence}; {replicate}"
            )
        )

    # New in v7:
    # Control/Uninfected + replicate is sufficient even without a timepoint.
    if replicate:
        return (
            "Tier2_control_replicate",
            (
                "Control/Uninfected + replicate: "
                f"{replicate}"
            )
        )

    # Without mock/replicate, keep requiring a timepoint before
    # using broader infection/challenge design evidence.
    if not time_evidence:
        return ("", "")

    infection_design = _first_regex_match(
        text,
        COMPILED_CONTROL_INFECTION_DESIGN
    )

    if infection_design:
        return (
            "Tier2_control_infection_design",
            (
                "Control/Uninfected + timepoint + "
                "infection/challenge evidence: "
                f"{time_evidence}; {infection_design}"
            )
        )

    return ("", "")


############################################################
# 6. Keys useful for Infection_status
############################################################

INFECTION_STATUS_KEYS = {

    "infection",
    "infection_status",
    "infected",
    "infect",
    "status",
    "health_state",
    "disease_state",
    "disease_status",

    "virus_infection",

    "sars_cov_2_infected",
    "sars_cov_2_infection",

    "h6n2_infection",

    "treatment/infection",
    "ex_vivo_infection_or_control",
}


############################################################
# 7. Patterns for infection status
############################################################

CONTROL_PATTERNS = [

    r"^control$",
    r"^-control-$",

    r"\buninfected\b",
    r"\bnon[- ]infected\b",

    r"\bmock\b",
    r"\bmock[- ]infected\b",

    r"\bhealthy\b",

]


INFECTED_PATTERNS = [

    r"^infected$",

    r"\binfected[ _-]+with\b",

    r"\bvirus[- ]?infected\b",
    r"\bviral[ _-]+infection\b",

    r"\binfection[ _-]+positive\b",

]


############################################################
# 8. Attribute keys that may contain an infection agent
#
# Higher-confidence keys come first.
############################################################

AGENT_KEYS_PRIMARY = [

    "infected_with",
    "infectionagent",
    "infection_agent",
    "infectious_agent",
    # "agent",
    "pathogen",

    "challenge",

    "virus",
    "viral_strain",
]


############################################################
# Lower-confidence fields that can sometimes contain
# the pathogen name.
############################################################

AGENT_KEYS_SECONDARY = [

    "infection",
    "infection_status",
    "virus_infection",
    # "disease",

]


############################################################
# 9. Values that should NOT become Infection_agent_raw
############################################################

GENERIC_AGENT_VALUES = {
    "",
    "infected",
    "infection",
    "uninfected",
    "control",
    "-control-",
    "healthy",
    "mock",
    "mock infected",
    "yes",
    "no",
    "positive",
    "negative",

    # add
    "inf",
    "uninf",
    "con",
    "kontrol",
    "treated",
    "pbs",
    "pbs buffer",
    "no infection",
    "no virus infection",
    "uninfected control",
    "none (uninfected control)",
}


############################################################
# Infection agent type patterns
############################################################

VIRUS_PATTERNS = [
    r"\bvirus\b",
    r"\bviral\b",
    r"\bSARS[- ]?CoV[- ]?2\b",
    r"\bSIV\b",
    r"\bSIVmac\d*\b",
    r"\bHBV\b",
    r"\bHCV\b",
    r"\bEBOV\b",
    r"\bEbola\b",
    r"\bPRRSV?\b",
    r"\bASFV?\b",
    r"\bAIV\b",
    r"\bH[1-9]N[1-9]\b",
    r"\bILTV\b",
    r"\bIBV\b",
    r"\bBDBV\b",
    r"\bZika\b",
    r"\bChikungunya\b",
    r"\bPowassan\b",
    r"\bKyasanur\b",
    r"\bDHAV[- ]?1\b",
    # manually added
    r"\bA/[A-Za-z0-9._-]+/[A-Za-z0-9._-]+/\d{2,4}\b",
    r"\bPRVABC59\b",
    r"\bVN1203\b",
    r"\bBC500\b",
    r"\bSendai(?:[ _-]+virus)?[ _-]+Cantell\b",
    r"\bpeste[ _-]+des[ _-]+petits[ _-]+ruminants\b",
    r"\bSFTSV\b",
    r"\bH1N1pdm09\b",
    r"\bALV[- ]?J\b",
    r"\bVSV\b",
    r"\bEMCV\b",
    r"\bReo(?:virus)?\b",
    r"\bPCV2\b",
    r"\bGPCMV\b",
    r"\bCoxsackievirus\b",
    r"\bPRV3M\b",
    r"\bRatPyV2\b",
    r"\bSV40\b",
    r"\bIBDV\b",
    r"\bInfluenza\b",
    r"\bBetacoronavirus\b",
    r"\bORFV\b",
    r"\bMARV371bat\b",
]


BACTERIUM_PATTERNS = [
    r"\bE\.?\s*coli\b",
    r"\bEscherichia coli\b",

    r"\bMycobacterium\b",
    r"\bM\.?\s*bovis\b",
    r"\bM\.?\s*tuberculosis\b",
    r"\bM\.?\s*leprae\b",

    r"\bMycobacterium avium\b",
    r"\bparatuberculosis\b",

    r"\bM\.?\s*gallisepticum\b",
    r"\bMycoplasma\b",

    r"\bSalmonella\b",
    r"\bS\.?\s*pullorum\b",

    r"\bF\.?\s*nucleatum\b",
    r"\bFusobacterium\b",

    r"\bBordetella\b",
    r"\bStaphylococcus\b",
    r"\bStreptococcus\b",
    r"\bCampylobacter\b",
    r"\bClostridium\b",
    # manually added
    r"\bJohne['’]?s[ _-]+Disease\b",
    r"\bTuberculosis,[ _-]*Bovine\b",
    r"\bBovine[ _-]+Tuberculosis\b",
    r"\bBorreliella\b",
    r"\bBorrelia\b",
    r"\bPseudomonas aeruginosa\b",
]


PROTOZOAN_PATTERNS = [
    r"\bPlasmodium\b",
    r"\bEimeria\b",
    r"\bLeishmania\b",
    r"\bNeospora\b",
    r"\bNc[- ]?Spain",
    r"\bToxoplasma\b",
    r"\bcoccidiosis\b", # protozoan disease
    r"\bE\.?\s*tenella\b",
]


HELMINTH_PATTERNS = [
    r"\bTeladorsagia\b",
    r"\bTrichuris\b",
    r"\bAnisakis\b",
    r"\bhelminth\b",
    r"\bnematode\b",
    r"\bAscaris suum\b",
    r"\bA\.?\s*suum\b",
]


FUNGUS_PATTERNS = [
    r"\bfungus\b",
    r"\bfungal\b",
    r"\bCandida\b",
    r"\bAspergillus\b",
]

############################################################
# 10. Compile regex
############################################################

COMPILED_EXPERIMENTAL_PATTERNS = [
    re.compile(
        x,
        flags=re.IGNORECASE
    )
    for x in EXPERIMENTAL_PATTERNS
]


COMPILED_NON_EXPERIMENTAL_PATTERNS = [
    re.compile(
        x,
        flags=re.IGNORECASE
    )
    for x in NON_EXPERIMENTAL_PATTERNS
]


COMPILED_CONTROL_PATTERNS = [
    re.compile(
        x,
        flags=re.IGNORECASE
    )
    for x in CONTROL_PATTERNS
]


COMPILED_INFECTED_PATTERNS = [
    re.compile(
        x,
        flags=re.IGNORECASE
    )
    for x in INFECTED_PATTERNS
]


############################################################
# 11. Helpers
############################################################

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


def is_missing(value):

    if pd.isna(value):
        return True

    return (
        str(value)
        .strip()
        .lower()
        in MISSING_VALUES
    )

def match_any(text, patterns):

    return any(
        re.search(
            pattern,
            text,
            flags=re.IGNORECASE
        )
        for pattern in patterns
    )


def classify_agent_type(agent):

    if is_missing(agent):
        return "Unknown"

    text = str(agent).strip()

    if match_any(
        text,
        VIRUS_PATTERNS
    ):
        return "Virus"

    if match_any(
        text,
        BACTERIUM_PATTERNS
    ):
        return "Bacterium"

    if match_any(
        text,
        PROTOZOAN_PATTERNS
    ):
        return "Protozoan"

    if match_any(
        text,
        HELMINTH_PATTERNS
    ):
        return "Helminth"

    if match_any(
        text,
        FUNGUS_PATTERNS
    ):
        return "Fungus"

    return "Other/Unclassified"


############################################################
# Parse:
#
# key=value | key=value | key=value
############################################################

def parse_characteristics(text):

    if is_missing(text):
        return []

    attributes = []

    for item in str(text).split(" | "):

        item = item.strip()

        if not item:
            continue

        if "=" not in item:
            continue

        key, value = item.split(
            "=",
            1
        )

        key = key.strip()
        value = value.strip()

        if not key:
            continue

        attributes.append(
            (
                normalize_key(key),
                value,
                key
            )
        )

    return attributes


############################################################
# Check whether any regex matches
############################################################

def any_pattern_match(
    text,
    compiled_patterns
):

    if not text:
        return False

    return any(
        pattern.search(text)
        for pattern in compiled_patterns
    )


############################################################
# 12. Experimental infection classification
############################################################

def classify_experimental_infection(
    attributes
):

    ########################################################
    # Priority 1:
    # Strong experimental key
    ########################################################

    for key_norm, value, key_original in attributes:

        if key_norm in STRONG_EXPERIMENTAL_KEYS:

            # Ignore completely missing values
            if is_missing(value):
                continue

            return (
                "Yes",
                f"{key_original}={value}",
                "Experimental_key"
            )


    ########################################################
    # Priority 2:
    # Explicit experimental phrase
    #
    # Search values rather than the attribute key.
    ########################################################

    for key_norm, value, key_original in attributes:

        if any_pattern_match(
            value,
            COMPILED_EXPERIMENTAL_PATTERNS
        ):

            return (
                "Yes",
                f"{key_original}={value}",
                "Experimental_phrase"
            )


    ########################################################
    # Priority 3:
    # Explicit natural / field infection
    #
    # Only use sample-specific metadata fields.
    ########################################################
    
    for key_norm, value, key_original in attributes:
    
        if key_norm not in NON_EXPERIMENTAL_KEYS:
            continue
    
        if any_pattern_match(
            value,
            COMPILED_NON_EXPERIMENTAL_PATTERNS
        ):
            return (
                "No",
                f"{key_original}={value}",
                "Non_experimental_phrase"
            )


    ########################################################
    # Insufficient evidence
    ########################################################

    return (
        "Unknown",
        "",
        "Insufficient_evidence"
    )

def detect_experimental_infection_from_freetext(row):

    text = str(
        row.get("All_characteristics", "")
    ).strip()

    if not text:
        return ""

    for label, pattern in EXPERIMENTAL_FREETEXT_PATTERNS.items():

        match = re.search(
            pattern,
            text,
            flags=re.IGNORECASE
        )

        if match:
            return (
                f"{label}: "
                f"{match.group(0)}"
            )

    return ""

############################################################
# 13. Infection status
#
# This is independent of Experimental_infection.
#
# Possible results:
#
# Infected
# Control/Uninfected
# Unknown
############################################################

def classify_infection_status(
    attributes
):

    ########################################################
    # First inspect explicit infection-status keys
    ########################################################

    status_attributes = [
        x
        for x in attributes
        if x[0] in INFECTION_STATUS_KEYS
    ]


    ########################################################
    # Control / uninfected has priority.
    #
    # This is important because samples from infection
    # experiments may be experimental controls.
    ########################################################

    for key_norm, value, key_original in status_attributes:

        if any_pattern_match(
            value,
            COMPILED_CONTROL_PATTERNS
        ):

            return (
                "Control/Uninfected",
                f"{key_original}={value}"
            )


    ########################################################
    # Explicit infected
    ########################################################

    for key_norm, value, key_original in status_attributes:

        if any_pattern_match(
            value,
            COMPILED_INFECTED_PATTERNS
        ):

            return (
                "Infected",
                f"{key_original}={value}"
            )


    ########################################################
    # Some experimental agent keys strongly indicate that
    # the sample is infected, but we intentionally do NOT
    # infer infection merely from challenge timing.
    #
    # Example:
    #
    # week_post_challenge=0
    #
    # could correspond to a baseline/control sample.
    ########################################################

    return (
        "Unknown",
        ""
    )


############################################################
# 14. Infection agent extraction
############################################################

def clean_agent_value(value):

    value = str(value).strip()

    if value.lower() in GENERIC_AGENT_VALUES:
        return ""

    if value.lower() in MISSING_VALUES:
        return ""

    return value


def extract_infection_agent(
    attributes
):

    ########################################################
    # Primary agent fields
    ########################################################

    for target_key in AGENT_KEYS_PRIMARY:

        for key_norm, value, key_original in attributes:

            if key_norm != target_key:
                continue

            cleaned = clean_agent_value(value)
            
            if cleaned and classify_agent_type(cleaned) != "Other/Unclassified":
                return (
                    cleaned,
                    key_original
                )


    ############################################################
    # Secondary fields
    ############################################################
    
    for target_key in AGENT_KEYS_SECONDARY:
    
        for key_norm, value, key_original in attributes:
    
            if key_norm != target_key:
                continue
    
            cleaned = clean_agent_value(value)
    
            if not cleaned:
                continue
    
            # Secondary fields are accepted only when
            # the value itself looks like an infectious agent.
            if classify_agent_type(cleaned) == "Other/Unclassified":
                continue
    
            return (
                cleaned,
                key_original
            )

    return (
        "",
        ""
    )


############################################################
# Strong experimental metadata rescue
############################################################

def rescue_strong_experimental_metadata(attributes):

    """
    Detect strong metadata showing that the sample belongs
    to an experimental infection/challenge design.

    Returns
    -------
    rescue : bool
    evidence : str or None
    """

    for key_norm, value, key_original in attributes:

        if not is_strong_experimental_key(key_norm):
            continue

        value_clean = str(value).strip()

        if not value_clean:
            continue

        # Ignore obvious missing values
        if is_missing(value_clean):
            continue

        evidence = (
            f"{key_original}={value_clean}"
        )

        return True, evidence

    return False, None

############################################################
# 15. Process one BioSample
############################################################

def classify_infection(row):

    all_characteristics = row.get(
        "All_characteristics",
        ""
    )

    attributes = parse_characteristics(
        all_characteristics
    )


    ########################################################
    # Experimental infection
    ########################################################

    (
        experimental,
        experimental_evidence,
        experimental_status
    ) = classify_experimental_infection(
        attributes
    )


    ########################################################
    # Infection status
    ########################################################

    (
        infection_status,
        infection_status_evidence
    ) = classify_infection_status(
        attributes
    )

    if (
        experimental == "Unknown"
        and infection_status in {
            "Infected",
            "Control/Uninfected"
        }
    ):
    
        rescue_evidence = detect_experimental_infection_from_freetext(
            row
        )
    
        if rescue_evidence:
    
            experimental = "Yes"
    
            experimental_status = (
                "Experimental_freetext_rescue"
            )
    
            experimental_evidence = rescue_evidence



    ########################################################
    # Tier-2 rescue:
    # unresolved Infected + explicit experimental time point
    #
    # We intentionally do NOT classify "infection=Infected"
    # alone as experimental infection because natural
    # infection remains possible.
    ########################################################

    if (
        experimental == "Unknown"
        and infection_status == "Infected"
    ):

        tier2_evidence = detect_tier2_infected_timepoint(
            attributes,
            all_characteristics
        )

        if tier2_evidence:

            experimental = "Yes"

            experimental_status = (
                "Tier2_infected_timepoint"
            )

            experimental_evidence = tier2_evidence




    ########################################################
    # Tier-2 experimental-control rescue
    #
    # Only unresolved samples already identified as
    # Control/Uninfected are considered here.
    ########################################################

    if (
        experimental == "Unknown"
        and infection_status == "Control/Uninfected"
    ):

        (
            control_rescue_status,
            control_rescue_evidence
        ) = detect_tier2_control_design(
            attributes,
            all_characteristics
        )

        if control_rescue_status:

            # Set to Yes temporarily. The final reconciliation
            # immediately below converts an experimentally
            # supported Control/Uninfected sample to No/Control.
            experimental = "Yes"
            experimental_status = control_rescue_status
            experimental_evidence = control_rescue_evidence


    ########################################################
    # Tier 3A
    # Explicit infection=Infected rescue
    ########################################################

    if (
        experimental == "Unknown"
        and infection_status == "Infected"
    ):

        for key_norm, value, key_original in attributes:

            if (
                key_norm == "infection"
                and str(value).strip().lower() == "infected"
            ):

                experimental = "Yes"

                experimental_status = (
                    "Tier3_explicit_infected_status"
                )

                experimental_evidence = (
                    f"{key_original}={value}"
                )

                break


    ########################################################
    # Tier 3B
    # Strong experimental metadata rescue
    ########################################################

    if experimental == "Unknown":

        (
            strong_experimental,
            strong_evidence
        ) = rescue_strong_experimental_metadata(
            attributes
        )

        if strong_experimental:

            if infection_status == "Infected":

                experimental = "Yes"

                experimental_status = (
                    "Tier3_strong_metadata_infected"
                )

                experimental_evidence = strong_evidence

            elif infection_status == "Control/Uninfected":

                experimental = "No/Control"

                experimental_status = (
                    "Tier3_strong_metadata_control"
                )

                experimental_evidence = strong_evidence
                
    ########################################################
    # Final reconciliation:
    # Experimental infection -> 3 classes
    #
    # Yes        : experimentally infected sample
    # No/Control : sample belongs to an experimental infection
    #              setting but is explicitly control/uninfected/mock
    # Unknown    : insufficient evidence
    #
    # Important:
    # Control/Uninfected samples are moved to No/Control only
    # when experimental infection evidence had already been detected.
    # A merely healthy/uninfected sample without experimental evidence
    # remains Unknown.
    ########################################################


    ########################################################
    # Rescue unresolved experimental controls
    ########################################################
    
    if (
        experimental == "Unknown"
        and infection_status == "Control/Uninfected"
    ):
    
        ####################################################
        # 1. Existing Tier2 control-design rescue
        ####################################################
    
        (
            control_rescue_status,
            control_rescue_evidence
        ) = detect_tier2_control_design(
            attributes,
            all_characteristics
        )
    
        if control_rescue_status:
    
            experimental = "No/Control"
            experimental_status = control_rescue_status
            experimental_evidence = control_rescue_evidence
    
    
        ####################################################
        # 2. Strong experimental-infection metadata
        #
        # Example:
        # weeks_pre(-)_or_post(+)_experimental_infection
        ####################################################
    
        if experimental == "Unknown":
    
            (
                strong_experimental,
                strong_evidence
            ) = rescue_strong_experimental_metadata(
                attributes
            )
    
            if strong_experimental:
    
                experimental = "No/Control"
    
                experimental_status = (
                    "Tier3_strong_metadata_control"
                )
    
                experimental_evidence = (
                    strong_evidence
                )

    ########################################################
    # Final 3-class reconciliation
    #
    # If a sample has evidence that it belongs to an
    # experimental infection design but its explicit
    # infection status is Control/Uninfected, the final
    # class must be No/Control rather than Yes.
    #
    # This is intentionally done AFTER all rescue tiers so
    # that Experimental_key / Experimental_phrase /
    # freetext / Tier2 evidence cannot leave experimental
    # controls mislabeled as Yes.
    ########################################################

    if (
        experimental == "Yes"
        and infection_status == "Control/Uninfected"
    ):
        experimental = "No/Control"


    ########################################################
    # Infection agent
    ########################################################

    (
        infection_agent,
        infection_agent_source
    ) = extract_infection_agent(
        attributes
    )


    return pd.Series({

        "Experimental_infection":
            experimental,

        "Experimental_infection_status":
            experimental_status,

        "Experimental_infection_evidence":
            experimental_evidence,

        "Infection_status":
            infection_status,

        "Infection_status_evidence":
            infection_status_evidence,

        "Infection_agent_raw":
            infection_agent,

        "Infection_agent_source":
            infection_agent_source,

    })


############################################################
# 16. Read metadata
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


required_columns = {
    "BioSample",
    "All_characteristics",
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


############################################################
# 17. Classify
############################################################

print(
    "Classifying infection metadata..."
)

results = df.apply(
    classify_infection,
    axis=1
)


############################################################
# 18. Add results
############################################################

for col in results.columns:

    df[col] = results[col]

df["Infection_agent_type"] = (
    df["Infection_agent_raw"]
    .apply(classify_agent_type)
)

############################################################
# 19. Save complete metadata
############################################################

df.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)


############################################################
# 20. Summary:
# Experimental infection
############################################################

print()

print(
    "=== Experimental infection ==="
)

print(
    df[
        "Experimental_infection"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)


print()

print(
    "=== Experimental infection status ==="
)

print(
    df[
        "Experimental_infection_status"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)


############################################################
# 21. Summary:
# Infection status
############################################################

print()

print(
    "=== Infection status ==="
)

print(
    df[
        "Infection_status"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)


############################################################
# 22. Summary:
# Infection agents
############################################################

print()

print(
    "=== Infection agent source ==="
)

print(
    df.loc[
        df["Infection_agent_raw"] != "",
        "Infection_agent_source"
    ]
    .value_counts()
    .to_string()
)


agent_counts = (
    df.loc[
        df["Infection_agent_raw"] != "",
        [
            "Infection_agent_raw",
            "Infection_agent_source"
        ]
    ]
    .value_counts()
    .reset_index(
        name="Count"
    )
)


agent_counts.to_csv(
    AGENT_COUNTS_FILE,
    sep="\t",
    index=False
)


print()

print(
    "=== Most common infection agents ==="
)

print(
    agent_counts
    .head(100)
    .to_string(
        index=False
    )
)


############################################################
# 23. Save Unknown cases for review
############################################################

unknown = df[
    df["Experimental_infection"]
    == "Unknown"
].copy()


review_columns = [
    "BioSample",
    "Tissue",
    "Cell_type",
    "Source_name",
    "Isolation_source",
    "Experimental_infection",
    "Experimental_infection_status",
    "Experimental_infection_evidence",
    "Infection_status",
    "Infection_status_evidence",
    "Infection_agent_raw",
    "Infection_agent_source",
    "All_characteristics",
]


review_columns = [
    x
    for x in review_columns
    if x in df.columns
]


unknown[
    review_columns
].to_csv(
    UNKNOWN_FILE,
    sep="\t",
    index=False
)


############################################################
# 24. Save all informative cases for manual review
############################################################

review = df[
    (
        df["Experimental_infection"]
        != "Unknown"
    )
    |
    (
        df["Infection_status"]
        != "Unknown"
    )
    |
    (
        df["Infection_agent_raw"]
        != ""
    )
].copy()


review[
    review_columns
].to_csv(
    REVIEW_FILE,
    sep="\t",
    index=False
)



for agent_type in [
    "Virus",
    "Bacterium",
    "Protozoan",
    "Helminth",
    "Fungus",
    "Other/Unclassified",
]:

    subset = df[
        df["Infection_agent_type"] == agent_type
    ]

    if len(subset) == 0:
        continue

    print()
    print(
        f"=== {agent_type}: most common agents ==="
    )
    print(
        subset["Infection_agent_raw"]
        .value_counts()
        .head(50)
        .to_string()
    )

############################################################
# 25. Finish
############################################################

print()

print("Saved:")
print(OUTPUT_FILE)
print(AGENT_COUNTS_FILE)
print(UNKNOWN_FILE)
print(REVIEW_FILE)

print("\n=== Experimental freetext rescue evidence ===")

rescue = df[
    df["Experimental_infection_status"] ==
    "Experimental_freetext_rescue"
].copy()

print("Rescued BioSamples:", len(rescue))

print(
    rescue["Experimental_infection_evidence"]
    .value_counts()
    .head(100)
    .to_string()
)

rescue.to_csv(
    "experimental_freetext_rescue_review.tsv",
    sep="\t",
    index=False
)

print("\n=== Experimental infection x Infection status ===")

print(
    pd.crosstab(
        df["Experimental_infection"],
        df["Infection_status"],
        margins=True
    )
)
remaining_status = df[
    (df["Experimental_infection"] == "Unknown") &
    (
        df["Infection_status"].isin([
            "Infected",
            "Control/Uninfected"
        ])
    )
].copy()

print("\n=== Infection status present but experimental still Unknown ===")
print(
    remaining_status["Infection_status"]
    .value_counts()
)

remaining_status.to_csv(
    "infection_status_but_experimental_unknown.tsv",
    sep="\t",
    index=False
)

print(
    df["Experimental_infection"]
    .value_counts(dropna=False)
)

print(
    df["Experimental_infection_status"]
    .value_counts(dropna=False)
)

print(
    pd.crosstab(
        df["Experimental_infection"],
        df["Infection_status"],
        margins=True
    )
)

