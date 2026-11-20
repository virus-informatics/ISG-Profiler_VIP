# Fig. S7D: Gene expression changes associated with Hepatovirus infection in chicken liver samples, assessed by GSEA
# Add taxonomy to BLASTx top hits from `FigS7D_tophit.ipynb`; then resume the notebook.

basedir="/path/to/project"


path_out="$basedir/output/revise_final/FigS7D/blastx/"
db_path="$basedir/../../db/taxdmp"

cd "$path_out"

cat taxids.txt | taxonkit lineage --data-dir "$db_path" | taxonkit reformat -f "{d}" --lineage-field 2 --data-dir "$db_path" -F | taxonkit reformat -f "{f};{g};{s}" --lineage-field 2 --data-dir "$db_path" -F | cut -f 1,3,4 > ./taxids_lineage.txt
cat taxids.txt | taxonkit lineage --data-dir "$db_path" | taxonkit reformat -f "{d};{K};{p};{c};{o};{f};{g};{s}" --lineage-field 2 --data-dir "$db_path" -F | cut -f 1,3 > taxids_long.txt
