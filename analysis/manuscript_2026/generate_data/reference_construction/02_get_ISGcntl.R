#!/usr/bin/env Rscript
rm( list=ls(all=TRUE) ) # clean up R workspace
library("tidyverse")
library("dplyr")
library(data.table)
setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "\nUsage:\n",
      "Rscript get_ISGcntl.R ",
      "<ISGcntl_geneID> ",
      "<output_gene2refseq> ",
      "<output_isoform_list>\n"
    )
  )
}


isg_file <- args[1]
outfile_table    <- args[2]
outfile_isoform  <- args[3]

# Fixed reference
gene2refseq_file <- "gene2refseq_Amniota_all_human_orthologs.list"

cat("start\n")
# import gene2refseq&ortholog list
ID_iso_info_Amni_ortho <- fread(gene2refseq_file, sep = "\t", header = TRUE, nThread = 1)
cat("fread OK\n")

# Extract isoform information for ISGs and control genes
ISGcntl <- read.table(isg_file, header = T, sep = "\t")
ID_iso_info_ISGcntl <- ID_iso_info_Amni_ortho %>%
  inner_join(
    ISGcntl,
    by = "hum_geneID"
  )%>%
  select(
    tax_id,
    Other_GeneID,
    Isoform,
    Other_symbol,
    hum_geneID,
    hum_symbol,
    type
)

write.table(
  ID_iso_info_ISGcntl,
  outfile_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ID_iso_info_ISGcntl %>%
  select(Isoform) %>%
  write.table(
    file = outfile_isoform,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

cat("finished\n")
