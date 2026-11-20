# Fig. 6B: ISG induction associated with Chaphamaparvovirus infection
# Add taxonomy to the two batches prepared by `Fig6B_tophit.ipynb`; then resume the notebook.

basedir="/path/to/project"


output_base="$basedir/output/revise_final/Fig6B/blastx"
db_path="$basedir/../../db/taxdmp"

for batch in 260209 260302; do
    path_out="$output_base/${batch}_blastx/"
    cd "$path_out"

    cat taxids.txt | taxonkit lineage --data-dir "$db_path" | taxonkit reformat -f "{d}" --lineage-field 2 --data-dir "$db_path" -F | taxonkit reformat -f "{f};{g};{s}" --lineage-field 2 --data-dir "$db_path" -F | cut -f 1,3,4 > ./taxids_lineage.txt
    cat taxids.txt | taxonkit lineage --data-dir "$db_path" | taxonkit reformat -f "{d};{K};{p};{c};{o};{f};{g};{s}" --lineage-field 2 --data-dir "$db_path" -F | cut -f 1,3 > taxids_long.txt
done
