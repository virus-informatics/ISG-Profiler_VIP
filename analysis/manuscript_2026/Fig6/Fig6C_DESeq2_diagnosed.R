# Fig. 6C: Gene expression changes associated with Chaphamaparvovirus infection in chicken liver samples
# Run DESeq2 for the diagnosed comparison using tables from `Fig6C_concat_featurecounts.ipynb`.
# Analyze the results with `Fig6C_GSEA.ipynb`.

rm( list=ls(all=TRUE) )
library("tidyverse")
library("DESeq2")
packageVersion("DESeq2")

# Input and output paths
basedir <- "/path/to/project"

output_dir <- file.path(basedir, "output/revise_final/Fig6C/deseq2")
# Input files
count_dir <- file.path(basedir, "output/revise_final/Fig6C/counts")
count_file <- file.path(count_dir, "merged_four_project_all_counts_ver2.txt")
sample_info <- file.path(count_dir, "sample_metadata.txt")
go_file <- file.path(
  basedir,
  "data/251215_rna_seq/To_kurihara_downloaded260209_chaphama/go/gene_orthologs_hs_gallus_table.txt"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

positive_samples <- c("Chap.miya", "ChickenLiver1to3")
negative_projects <- c("PRJNA577590", "PRJNA612882", "PRJNA622813")
comparison_name <- "Fujino"

count <- read.table(count_file, sep = "\t", header = TRUE,
                    row.names = 1, check.names = FALSE)
sample <- read.table(sample_info, sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE)

# Compare the two specified positive samples with negative samples from three projects
sample <- sample %>%
  filter(Run %in% positive_samples |
           (BioProject %in% negative_projects & Virus_status == "Negative")) %>%
  mutate(con = ifelse(Run %in% positive_samples, "Positive", "Negative")) %>%
  arrange(con)

stopifnot(
  all(positive_samples %in% sample$Run),
  !anyDuplicated(sample$Run)
)

out <- sample %>%
  group_by(con) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(ratio = 100 * n / sum(n))
print(comparison_name)
print(out)

# Compare positive and negative samples

count <- count %>% select(all_of(sample$Run))
count <- as.matrix(count)

group <- data.frame(con = factor(sample$con, levels = c("Negative", "Positive")))
rownames(group) <- sample$Run
stopifnot(identical(colnames(count), rownames(group)))

dds <- DESeqDataSetFromMatrix(countData = count, colData = group, design = ~ con)
dds <- DESeq(dds)

res <- results(dds, contrast = c("con", "Positive", "Negative"))
head(res)

out.pdfname <- file.path(output_dir, paste0("plotMA_", comparison_name, "_Positive_vs_Negative.pdf"))
pdf(out.pdfname, width = 5, height = 5)
plotMA(res, alpha = 0.01)
dev.off()

f.name <- file.path(output_dir, paste0("result_", comparison_name, "_Positive_vs_Negative.txt"))
write.table(res, file = f.name, row.names = T, col.names = T, sep = "\t", quote=F)

ortholog <- read.table(go_file, sep = "\t", header = T)

res <- res %>% as.data.frame() %>% mutate(Symbol.2 = rownames(res))
res.hs.Id <- res %>% left_join(ortholog, by = "Symbol.2")
res.hs.Id <- res.hs.Id %>% select(Symbol.1, Symbol.2, 1:6)
res.hs.Id <- res.hs.Id %>% dplyr::rename(gene_name_hs = Symbol.1, gene_name_gallus = Symbol.2)

f.name <- file.path(output_dir, paste0("result_", comparison_name, "_Positive_vs_Negative_hs_gallus_and_Gene.txt"))
write.table(res.hs.Id, file = f.name, row.names = F, col.names = T, sep = "\t", quote=F)
