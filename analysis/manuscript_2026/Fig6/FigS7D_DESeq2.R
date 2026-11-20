# Fig. S7D: Gene expression changes associated with Hepatovirus infection in chicken liver samples, assessed by GSEA
# Run DESeq2 using sample labels from `FigS7D_tophit.ipynb`.
# Analyze the results with `FigS7D_GSEA.ipynb`.

rm( list=ls(all=TRUE) )
library("tidyverse")
library("DESeq2")
packageVersion("DESeq2")

# Define paths
basedir <- "/path/to/project"

output_dir <- file.path(basedir, "output/revise_final/FigS7D/deseq2")
# Input files
count_dir <- file.path(basedir, "data/251215_rna_seq/To_kurihara_hepatovirus/featureCount")
input_dir <- file.path(basedir, "output/revise_final/FigS7D/metadata")
sample_metadata <- file.path(input_dir, "hepatovirus_metadata_genus.tsv")
go_file <- file.path(basedir, "data/251215_rna_seq/To_kurihara_hepatovirus/go/gene_orthologs_hs_gallus_table.txt")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (i in c("PRJEB63475", "PRJNA774885")) {
  count <- read.table(file.path(count_dir, paste("all_counts_",i,"_ver2.txt", sep = "")),
                      sep = "\t",  header = T, row.names = 1)
  
  sample <- read.table(sample_metadata, sep = "\t",
                       header = TRUE, stringsAsFactors = FALSE)
  
  sample <- sample %>%
    filter(BioProject == i) %>%
    dplyr::rename(Run = Sample) %>%
    filter((!is.na(Positive) & Positive == 1) | (!is.na(Negative) & Negative == 1))
  
  sample <- sample %>%
    mutate(con = ifelse(!is.na(Positive) & Positive == 1, "Positive",
                        ifelse(!is.na(Negative) & Negative == 1, "Negative", NA))) %>%
    filter(!is.na(con))
  
  sample <- sample %>% arrange(con)
  
  out <- sample %>%
    group_by(con) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(ratio = 100 * n / sum(n))
  print(i)
  print(out)
  
  # Compare positive and negative samples
  
  count <- count %>% select(one_of(sample$Run))
  count <- as.matrix(count)
  
  group <- data.frame(con = factor(sample$con, levels = c("Negative", "Positive")))
  rownames(group) <- sample$Run
  
  dds <- DESeqDataSetFromMatrix(countData = count, colData = group, design = ~ con)
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("con", "Positive", "Negative"))
  head(res)
  
  out.pdfname <- file.path(output_dir, paste("plotMA_", i, "_Positive_vs_Negative.pdf", sep = ""))
  pdf(out.pdfname, width = 5, height = 5)
  plotMA(res, alpha = 0.01)
  dev.off()
  
  f.name <- file.path(output_dir, paste("result_",i, "_Positive_vs_Negative.txt", sep = ""))
  write.table(res, file = f.name, row.names = T, col.names = T, sep = "\t", quote=F)
  
  ortholog <- read.table(go_file, sep = "\t",  header = T)
  
  res <- res %>% as.data.frame() %>% mutate(Symbol.2 = rownames(res))
  res.hs.Id <- res %>% left_join(ortholog, by = "Symbol.2")
  res.hs.Id <- res.hs.Id %>% select(Symbol.1, Symbol.2, 1:6)
  res.hs.Id <- res.hs.Id %>% dplyr::rename(gene_name_hs = Symbol.1, gene_name_gallus = Symbol.2)
  
  f.name <- file.path(output_dir, paste("result_",i, "_Positive_vs_Negative_hs_gallus_and_Gene.txt", sep = ""))
  write.table(res.hs.Id, file = f.name, row.names = F, col.names = T, sep = "\t", quote=F)
}

