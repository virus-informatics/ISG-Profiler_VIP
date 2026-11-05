#!/usr/bin/env python3

import os
import sys
import csv
import time
import json
import math
import requests
import pandas as pd

from concurrent.futures import ThreadPoolExecutor, as_completed


# ============================================================
# Settings
# ============================================================

#SRA Run accession ID list without header
INPUT_FILE = "logan_168438.list"

RUN_MAPPING_FILE = "run_to_biosample.tsv"
BIOSAMPLE_CACHE_FILE = "biosample_characteristics.tsv"
FINAL_OUTPUT_FILE = "SRA_characteristics.tsv"

ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
BIOSAMPLE_URL = "https://www.ebi.ac.uk/biosamples/samples"

# ENA batch size
RUN_BATCH_SIZE = 200

# BioSamples API parallelism
N_THREADS = 10

TIMEOUT = 60
MAX_RETRY = 5

# Save BioSample progress every N completed samples
SAVE_EVERY = 500

def request_with_retry(url, params=None):

    for attempt in range(1, MAX_RETRY + 1):

        try:
            r = requests.get(
                url,
                params=params,
                timeout=TIMEOUT
            )

            # retry server/rate-limit errors
            if r.status_code in {429, 500, 502, 503, 504}:
                wait = min(60, 2 ** attempt)
                print(
                    f"HTTP {r.status_code}; retry in {wait}s",
                    file=sys.stderr
                )
                time.sleep(wait)
                continue

            r.raise_for_status()
            return r

        except Exception as e:

            if attempt == MAX_RETRY:
                raise

            wait = min(60, 2 ** attempt)

            print(
                f"Request failed ({attempt}/{MAX_RETRY}): "
                f"{e}; retry in {wait}s",
                file=sys.stderr
            )

            time.sleep(wait)

    return None


# ============================================================
# Read Run accession list
# ============================================================

def read_run_list(path):

    with open(path) as f:
        runs = [
            line.strip().split()[0]
            for line in f
            if line.strip()
        ]

    # preserve order, remove duplicates
    runs = list(dict.fromkeys(runs))

    return runs


# ============================================================
# Step 1:
# Run accession -> BioSample mapping using ENA
# ============================================================

def query_ena_batch(run_list):

    query = " OR ".join(
        f'run_accession="{run}"'
        for run in run_list
    )

    params = {
        "result": "read_run",
        "query": query,
        "fields": (
            "run_accession,"
            "sample_accession,"
            "secondary_sample_accession,"
            "scientific_name"
        ),
        "format": "json",
        "limit": 0
    }

    r = request_with_retry(
        ENA_URL,
        params=params
    )

    data = r.json()

    return data


def get_run_mapping(runs):

    # --------------------------------------------------------
    # Resume if mapping file already exists
    # --------------------------------------------------------

    if os.path.exists(RUN_MAPPING_FILE):

        old = pd.read_csv(
            RUN_MAPPING_FILE,
            sep="\t",
            dtype=str
        )

        done_runs = set(
            old["Run_accession"]
            .dropna()
            .tolist()
        )

        print(
            f"Existing Run mapping: {len(done_runs)}",
            file=sys.stderr
        )

    else:

        old = pd.DataFrame()

        done_runs = set()


    remaining = [
        x for x in runs
        if x not in done_runs
    ]

    print(
        f"Runs remaining for ENA query: {len(remaining)}",
        file=sys.stderr
    )

    rows = []

    batches = [
        remaining[i:i + RUN_BATCH_SIZE]
        for i in range(
            0,
            len(remaining),
            RUN_BATCH_SIZE
        )
    ]

    for batch_i, batch in enumerate(batches, 1):

        print(
            f"ENA batch {batch_i}/{len(batches)}",
            file=sys.stderr
        )

        try:

            data = query_ena_batch(batch)

            found = set()

            for x in data:

                run = x.get("run_accession", "")

                found.add(run)

                rows.append({
                    "Run_accession": run,
                    "Sample_accession":
                        x.get("sample_accession", ""),
                    "BioSample":
                        x.get("secondary_sample_accession", ""),
                    "Species":
                        x.get("scientific_name", ""),
                    "Run_status": "OK"
                })

            # Keep missing runs too
            missing = set(batch) - found

            for run in missing:

                rows.append({
                    "Run_accession": run,
                    "Sample_accession": "",
                    "BioSample": "",
                    "Species": "",
                    "Run_status": "Not_found"
                })

        except Exception as e:

            print(
                f"ENA batch failed: {e}",
                file=sys.stderr
            )

            for run in batch:

                rows.append({
                    "Run_accession": run,
                    "Sample_accession": "",
                    "BioSample": "",
                    "Species": "",
                    "Run_status":
                        f"ERROR: {e}"
                })


        # append periodically
        if (
            batch_i % 10 == 0
            or batch_i == len(batches)
        ):

            new_df = pd.DataFrame(rows)

            if len(old) > 0:
                combined = pd.concat(
                    [old, new_df],
                    ignore_index=True
                )
            else:
                combined = new_df

            combined = (
                combined
                .drop_duplicates(
                    subset="Run_accession",
                    keep="last"
                )
            )

            combined.to_csv(
                RUN_MAPPING_FILE,
                sep="\t",
                index=False
            )

            old = combined
            rows = []

    return pd.read_csv(
        RUN_MAPPING_FILE,
        sep="\t",
        dtype=str
    )


# ============================================================
# Step 2:
# Fetch all biosample characteristics
# ============================================================

def extract_text_values(values):

    texts = []

    if not isinstance(values, list):
        return texts

    for value in values:

        if (
            isinstance(value, dict)
            and value.get("text")
        ):
            texts.append(
                str(value["text"]).strip()
            )

    return texts

def characteristics_to_string(characteristics):

    output = []

    for key, values in characteristics.items():

        texts = extract_text_values(values)

        if texts:

            output.append(
                f"{key}={'; '.join(texts)}"
            )

    return " | ".join(output)


def fetch_one_biosample(accession):

    if not accession:

        return {
            "BioSample": accession,
            "All_characteristics": "",
            "BioSample_status":
                "Missing_accession"
        }

    url = f"{BIOSAMPLE_URL}/{accession}"

    try:

        r = request_with_retry(url)

        bs = r.json()

        characteristics = bs.get(
            "characteristics",
            {}
        )


        all_characteristics = (
            characteristics_to_string(
                characteristics
            )
        )

        return {
            "BioSample": accession,
            "All_characteristics":
                all_characteristics,
            "BioSample_status": "OK"
        }

    except Exception as e:

        return {
            "BioSample": accession,
            "All_characteristics": "",
            "BioSample_status":
                f"ERROR: {e}"
        }


# ============================================================
# BioSample cache / resume
# ============================================================

def get_biosample_metadata(
    biosamples
):

    # --------------------------------------------------------
    # Existing cache
    # --------------------------------------------------------

    if os.path.exists(
        BIOSAMPLE_CACHE_FILE
    ):

        cache = pd.read_csv(
            BIOSAMPLE_CACHE_FILE,
            sep="\t",
            dtype=str
        )

        done = set(
            cache["BioSample"]
            .dropna()
            .tolist()
        )

        print(
            f"Cached BioSamples: {len(done)}",
            file=sys.stderr
        )

    else:

        cache = pd.DataFrame()

        done = set()


    remaining = [
        x for x in biosamples
        if x not in done
    ]

    print(
        f"BioSamples remaining: {len(remaining)}",
        file=sys.stderr
    )

    pending_rows = []

    with ThreadPoolExecutor(
        max_workers=N_THREADS
    ) as executor:

        future_map = {
            executor.submit(
                fetch_one_biosample,
                accession
            ): accession
            for accession in remaining
        }

        for i, future in enumerate(
            as_completed(future_map),
            1
        ):

            accession = future_map[future]

            try:

                row = future.result()

            except Exception as e:

                row = {
                    "BioSample": accession,
                    "All_characteristics": "",
                    "BioSample_status":
                        f"ERROR: {e}"
                }

            pending_rows.append(row)

            if (
                i % 100 == 0
                or i == len(remaining)
            ):

                print(
                    f"BioSamples: "
                    f"{i}/{len(remaining)}",
                    file=sys.stderr
                )


            # ------------------------------------------------
            # Save progress
            # ------------------------------------------------

            if (
                len(pending_rows)
                >= SAVE_EVERY
                or i == len(remaining)
            ):

                tmp = pd.DataFrame(
                    pending_rows
                )

                if len(cache) > 0:

                    combined = pd.concat(
                        [cache, tmp],
                        ignore_index=True
                    )

                else:

                    combined = tmp

                combined = (
                    combined
                    .drop_duplicates(
                        subset="BioSample",
                        keep="last"
                    )
                )

                combined.to_csv(
                    BIOSAMPLE_CACHE_FILE,
                    sep="\t",
                    index=False
                )

                cache = combined

                pending_rows = []


    return pd.read_csv(
        BIOSAMPLE_CACHE_FILE,
        sep="\t",
        dtype=str
    )


# ============================================================
# Main
# ============================================================

def main():

    runs = read_run_list(
        INPUT_FILE
    )

    print(
        f"Unique Run accessions: "
        f"{len(runs)}",
        file=sys.stderr
    )


    # ========================================================
    # Run -> BioSample
    # ========================================================

    run_df = get_run_mapping(runs)


    # ========================================================
    # Unique BioSamples
    # ========================================================

    biosamples = (
        run_df["BioSample"]
        .dropna()
        .loc[
            lambda x:
            (x != "")
            & (x != "nan")
        ]
        .drop_duplicates()
        .tolist()
    )

    print(
        f"Unique BioSamples: "
        f"{len(biosamples)}",
        file=sys.stderr
    )


    # ========================================================
    # BioSample metadata
    # ========================================================

    bio_df = get_biosample_metadata(
        biosamples
    )


    # ========================================================
    # Merge back to Run level
    # ========================================================

    final = run_df.merge(
        bio_df,
        on="BioSample",
        how="left"
    )


    # ========================================================
    # Save final table
    # ========================================================

    final.to_csv(
        FINAL_OUTPUT_FILE,
        sep="\t",
        index=False
    )


    # ========================================================
    # Summary
    # ========================================================

    total = len(final)


    print(
        "\n=================================",
        file=sys.stderr
    )

    print(
        f"Total Runs: {total}",
        file=sys.stderr
    )

    print(
        f"Final output: "
        f"{FINAL_OUTPUT_FILE}",
        file=sys.stderr
    )

    print(
        "=================================",
        file=sys.stderr
    )


if __name__ == "__main__":
    main()
