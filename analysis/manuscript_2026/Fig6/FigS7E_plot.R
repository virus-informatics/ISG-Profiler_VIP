# Fig. S7E: Expression changes of representative genes for each GSEA term shown in Fig. S7D
# Plot representative-gene expression and membership tables prepared by `FigS7E_table.ipynb`.

library(circlize)
library(tidyverse)
library(patchwork)

gs.interest.v <- c("WIELAND_UP_BY_HBV_INFECTION",
                   "Interferon Alpha Response",
                   "IL_6/JAK/STAT3 Signaling",
                   "Xenobiotic Metabolism",
                   "Adipogenesis")

project_Id.v <- c("PRJNA774885","PRJEB63475")

# Define paths
basedir <- "/path/to/project"

output_dir <- file.path(basedir, "output/revise_final/FigS7E")
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



mat_df.long <- mat_df.long %>% filter(symbol %in% symbol_order.v, project_Id %in% project_Id.v) %>%
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

# Allow space for gene rows and the vertical term labels (sizes in inches).
plot_width <- 5.0 + 0.25 * n_distinct(mat_df.long$project_Id)
plot_height <- max(5.0, 2.8 + 0.18 * length(symbol_order.v))

pdf("stat_log2FC_heatmap.pdf", width = plot_width, height = plot_height)
print(g1 + g2)
dev.off()
