#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 8
#$ -l lmem
#$ -l s_vmem=60G,mem_req=60G
#$ -o .
#$ -e .
# Fig. S6E: Detection of ISG non-inducing viral families by ISG-VIP
# Run BLASTn on the HPC for additional contigs exported by `FigS6E.ipynb`.
# Return the results before resuming the notebook.

base_path="/path/to/project"
cd ${base_path}

# Run BLASTn
module load /usr/local/package/modulefiles/blast+/2.15.0

blastn -query ${base_path}/output/260818_revise/4474_samples_rev_v2/figS5C/virus_hit_all_contigs_w_bl.fasta \
	-db /usr/local/db/blast/ncbi/core_nt \
	-out ${base_path}/output/260818_revise/4474_samples_rev_v2/figS5C/virus_hit_all_contigs_w_bl.txt \
	-evalue 1e-4 \
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe staxids stitle qlen slen qcovs qcovhsp" \
	-num_threads 8

