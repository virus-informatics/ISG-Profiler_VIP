#!/usr/bin/env python3

import argparse
import time
from Bio import Entrez

parser = argparse.ArgumentParser()
parser.add_argument("input", help="accession list, one accession per line")
parser.add_argument("output", help="output FASTA")
parser.add_argument(
    "--email",
    required=True,
    help="email address used for NCBI Entrez"
)
parser.add_argument(
    "--batch-size",
    type=int,
    default=200
)
args = parser.parse_args()

# NCBI Entrez settings
Entrez.email = args.email
Entrez.tool = "ISG_profiler_db_builder"

############################################################
# Read accession list
############################################################

with open(args.input) as f:
    accessions = [
        x.strip()
        for x in f
        if x.strip()
    ]

print("Total accessions:", len(accessions))


############################################################
# Download in batches
############################################################

with open(args.output, "w") as out:

    for start in range(
        0,
        len(accessions),
        args.batch_size
    ):

        batch = accessions[
            start:start + args.batch_size
        ]

        end = start + len(batch)

        print(
            f"Downloading {start + 1}-{end} "
            f"/ {len(accessions)}"
        )

        for attempt in range(5):

            try:

                handle = Entrez.efetch(
                    db="nuccore",
                    id=",".join(batch),
                    rettype="fasta",
                    retmode="text"
                )

                fasta = handle.read()
                handle.close()

                out.write(fasta)
                out.flush()

                break

            except Exception as e:

                print(
                    f"Attempt {attempt + 1} failed:",
                    e
                )

                if attempt == 4:
                    raise

                time.sleep(
                    5 * (attempt + 1)
                )

        # Conservative interval between requests
        time.sleep(0.4)


print("Finished.")
