# Fig. 5H: Enrichment of events not attributable to genuine infection in the FN fraction
# Prepare viral and host taxonomy tables for `Fig5HI.ipynb` and `FigS6F.ipynb`.

basedir="/path/to/project"

mkdir -p "$basedir/output/revise_final/Fig5HI"
cd "$basedir/output/revise_final/Fig5HI"

DB="$basedir/../../db/taxdmp"
IN="$basedir/../../db/virushostdb_260821/virushostdb.tsv"

# Virus taxonomy from column 1 of Virus–Host DB
{
    echo -e "virus_taxid\tvirus_domain\tvirus_phylum\tvirus_class\tvirus_order\tvirus_family\tvirus_genus\tvirus_species"

    tail -n +2 "$IN" \
    | cut -f1 \
    | sort -u \
    | taxonkit lineage --data-dir "$DB" \
    | taxonkit reformat \
        -f "{d}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
        --lineage-field 2 \
        --data-dir "$DB" \
        -F \
    | cut -f1,3-9
} > virus_taxonomy.tsv

# Host taxonomy from column 8 of Virus–Host DB
{
    printf 'host_taxid\thost_domain\thost_phylum\thost_class\thost_order\thost_family\thost_genus\thost_species\n'

    tail -n +2 "$IN" \
    | cut -f8 \
    | LC_ALL=C sort -u \
    | taxonkit lineage --data-dir "$DB" \
    | taxonkit reformat \
        -f '{d}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' \
        --lineage-field 2 \
        --data-dir "$DB" \
        -F \
    | cut -f1,3-9
} > host_taxonomy.tsv
