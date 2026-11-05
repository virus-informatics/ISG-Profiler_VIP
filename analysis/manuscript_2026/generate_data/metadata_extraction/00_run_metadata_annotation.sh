#!/usr/bin/env bash

set -euo pipefail

echo "Step 1: Retrieve BioSample metadata"
python3 01_get_biosample_info.py

echo "Step 2: Extract tissue and cell information"
python3 02_get_tissue_cell_info_ver2.py

echo "Step 3: Normalize tissue annotations"
python3 03_norm_tissue_ver3.py

echo "Step 4: Normalize cell-type annotations"
python3 04_norm_cell_ver2.py

echo "Step 5: Classify in vivo / in vitro status"
python3 05_classify_vivo_vitro_ver2.py

echo "Step 6: Classify infection status"
python3 06_classify_infection_v8.py

echo "Step 7: Merge annotations with SRA metadata"
python3 07_bind_sra_bio.py

echo "Done: SRA_metadata_annotated_v2.tsv"
