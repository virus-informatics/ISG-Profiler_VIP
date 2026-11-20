#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 8
#$ -l lmem
#$ -l s_vmem=60G,mem_req=60G
#$ -o .
#$ -e .
# Fig. 5: Detection of previously overlooked viral infections and flagging of non-genuine detections
# Run BLASTn on the HPC for contigs exported by `Fig5_table.ipynb`.
# Return the results before resuming the notebook.

base_path="/path/to/project"
cd ${base_path}

# Run BLASTn
module load /usr/local/package/modulefiles/blast+/2.15.0

blastn -query ${base_path}/data/250909_virus_hit/virus_hit_all_contigs.fasta \
	-db /usr/local/db/blast/ncbi/core_nt \
	-out ${base_path}/output/250909_blastn/virus_hit_all_contigs_megablast.txt \
	-evalue 1e-4 \
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe staxids stitle qlen slen qcovs qcovhsp" \
	-num_threads 8
