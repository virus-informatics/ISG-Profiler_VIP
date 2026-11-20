#!/usr/bin/env R
# Fig. S7B: Gene expression changes associated with Chaphamaparvovirus infection in chicken liver samples
# Plot representative-gene expression and membership tables prepared by `FigS7B_table.ipynb`.

library(tidyverse)
library(circlize)
library(patchwork)

gs.interest.v <- c("WIELAND_UP_BY_HBV_INFECTION",
                   "Interferon Alpha Response",
                   "IL_6/JAK/STAT3 Signaling",
                   "Inflammatory Response",
                   "Allograft Rejection",
                   "G2_M Checkpoint",
                   "Epithelial Mesenchymal Transition",
                   "Xenobiotic Metabolism",
                   "Adipogenesis",
                   "Bile Acid Metabolism")

project_Id.v <- c("Fujino", "PRJNA622813","PRJNA612882","PRJNA577590")


# Define paths
basedir <- "/path/to/project"

output_dir <- file.path(basedir, "output/revise_final/FigS7B")
input_dir <- output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

# Read the input tables
mat_df  <- read.csv(file.path(input_dir, "gene_project_log2FC_matrix.csv"), row.names = 1, check.names = FALSE)
term_bin_df <- read.csv(file.path(input_dir, "term_anno.csv"), row.names = 1, check.names = FALSE)

mat_df.long <- mat_df %>% mutate(symbol = rownames(mat_df)) %>% gather(key = project_Id, value = log2FC, -symbol)
term_bin_df.long <- term_bin_df %>% mutate(symbol = rownames(term_bin_df)) %>% gather(key = gs, value = membership, -symbol) %>% filter(membership == 1)

term_bin_df.long <- term_bin_df.long %>% filter(gs %in% gs.interest.v) %>% mutate(gs = factor(gs, levels = gs.interest.v)) %>% arrange(gs, symbol)

symbol_order.v <- term_bin_df.long %>% pull(symbol) %>% unique()

term_bin_df.long <- term_bin_df.long %>% mutate(symbol = factor(symbol, levels = rev(symbol_order.v)))



mat_df.long <- mat_df.long %>% filter(symbol %in% symbol_order.v) %>%
  mutate(symbol = factor(symbol, levels = rev(symbol_order.v)),
         project_Id = factor(project_Id, levels = project_Id.v))


col_fun <- colorRamp2(c(-2, 0, 2), c("#1B9E9E", "white", "#C51B7D"))

g1 <- ggplot(mat_df.long, aes(x = project_Id, y = symbol, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = col_fun(seq(-2, 2, length.out = 100)),
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

g2 <- ggplot(term_bin_df.long, aes(x = gs, y = symbol, fill = membership)) +
  geom_tile() +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

pdf("stat_log2FC_heatmap_four_projects.pdf", width = 10, height = 10)
print(g1 + g2)
dev.off()

