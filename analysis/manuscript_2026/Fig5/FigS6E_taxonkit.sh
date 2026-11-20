#!/bin/bash
# Fig. S6E: Detection of ISG non-inducing viral families by ISG-VIP
# Add taxonomy to BLASTn hits for `FigS6E.ipynb`; then resume the notebook.

basedir="/path/to/project"


path_out="$basedir/output/revise_final/FigS6E/"
db_path="$basedir/../../db/taxdmp"

cd "$path_out"

cat taxids_blastn.txt | taxonkit lineage --data-dir "$db_path" | taxonkit reformat -f "{d}" --lineage-field 2 --data-dir "$db_path" -F | taxonkit reformat -f "{f};{g};{s}" --lineage-field 2 --data-dir "$db_path" -F | cut -f 1,3,4 > ./taxids_lineage_blastn.txt
cat taxids_blastn.txt | taxonkit lineage --data-dir "$db_path" | taxonkit reformat -f "{d};{K};{p};{c};{o};{f};{g};{s}" --lineage-field 2 --data-dir "$db_path" -F | cut -f 1,3 > taxids_long_blastn.txt
