# Fig. 5: Detection of previously overlooked viral infections and flagging of non-genuine detections
# Add taxonomy to BLASTx top hits from `Fig5_table.ipynb`; then resume the notebook.

basedir="/path/to/project"


path_out="$basedir/output/revise_final/Fig5/"
db_path="$basedir/../../db/taxdmp"

cd "$path_out"

cat taxids.txt | taxonkit lineage --data-dir "$db_path" | taxonkit reformat -f "{d}" --lineage-field 2 --data-dir "$db_path" -F | taxonkit reformat -f "{f};{g};{s}" --lineage-field 2 --data-dir "$db_path" -F | cut -f 1,3,4 > ./taxids_lineage.txt
cat taxids.txt | taxonkit lineage --data-dir "$db_path" | taxonkit reformat -f "{d};{K};{p};{c};{o};{f};{g};{s}" --lineage-field 2 --data-dir "$db_path" -F | cut -f 1,3 > taxids_long.txt