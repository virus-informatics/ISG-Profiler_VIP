rm(list = ls(all = TRUE))

# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(scales)
set.seed(1234)

# Working directory
setwd(
  "SOMEWHERE"
)


# Read ISG data
ISG_meta <- fread(
  "Aves_Mam_mbio_250613.ISGcntl_norm_logan.txt",
  sep = "\t",
  select = c(
    "ID",
    "hum_symbol",
    "type",
    "ISG_score",
    "species",
    "Infection"
  )
)


# Keep gene-level ISG data
ISG_gene_only <- ISG_meta %>%
  filter(
    hum_symbol != "ISG score",
    type == "ISG"
  )


# Read Core ISG list
core_isgs <- read.table(
  "isg_59.list",
  header = FALSE,
  stringsAsFactors = FALSE
)$V1


# Select top 10 host species by sample number
top10_species <- ISG_gene_only %>%
  select(
    ID,
    species
  ) %>%
  distinct() %>%
  count(
    species,
    name = "n_sample"
  ) %>%
  arrange(
    desc(n_sample)
  ) %>%
  slice_head(
    n = 10
  )

print(top10_species)


# Keep top 10 host species
ISG_top10 <- ISG_gene_only %>%
  filter(
    species %in% top10_species$species
  )


# Calculate sample-level ISG score

sample_ISG_summary <- ISG_top10 %>%
  group_by(ID) %>%
  summarise(
    ISG_score_mean = mean(
      ISG_score,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

ISG_top10_with_score <- ISG_top10 %>%
  left_join(
    sample_ISG_summary,
    by = "ID"
  )


# Calculate correlation between each Core ISG and ISG score

core_ISG_cor <- ISG_top10_with_score %>%
  filter(
    hum_symbol %in% core_isgs
  ) %>%
  group_by(
    species,
    hum_symbol
  ) %>%
  summarise(
    n = sum(
      complete.cases(
        ISG_score,
        ISG_score_mean
      )
    ),
    rho = cor(
      ISG_score,
      ISG_score_mean,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )


# Order species by sample size

species_order <- top10_species %>%
  arrange(
    desc(n_sample)
  ) %>%
  pull(species)

core_ISG_cor <- core_ISG_cor %>%
  mutate(
    species = factor(
      species,
      levels = species_order
    )
  )

# Core ISG vs ISG score boxplot
p_core_box <- ggplot(
  core_ISG_cor,
  aes(
    x = species,
    y = rho
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    linewidth = 0.4
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.3,
    size = 0.7
  ) +
  scale_y_continuous(
    limits = c(-0.5, 1),
    breaks = seq(
      -0.5,
      1,
      0.5
    )
  ) +
  labs(
    x = "Host species",
    y = "Spearman's correlation coefficient"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    panel.grid = element_blank()
  )

print(p_core_box)

ggsave(
  "ISGscore_boxplot_260930.pdf",
  p_core_box,
  width = 10,
  height = 8
)


# Calculate pairwise Core ISG correlations within each species

cor_list <- list()

for (sp in species_order) {
  
  message(
    "Calculating correlation: ",
    sp
  )
  
  dat_sp <- ISG_top10 %>%
    filter(
      species == sp,
      hum_symbol %in% core_isgs
    ) %>%
    select(
      ID,
      hum_symbol,
      ISG_score
    ) %>%
    pivot_wider(
      names_from = hum_symbol,
      values_from = ISG_score
    )
  
  mat_sp <- dat_sp %>%
    select(-ID) %>%
    as.matrix()
  
  cor_sp <- cor(
    mat_sp,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  
  diag(cor_sp) <- NA
  
  cor_list[[sp]] <- cor_sp
}


# Create common Core ISG correlation array

gene_names <- sort(core_isgs)

cor_array <- array(
  NA_real_,
  dim = c(
    length(gene_names),
    length(gene_names),
    length(cor_list)
  ),
  dimnames = list(
    gene_names,
    gene_names,
    names(cor_list)
  )
)

for (i in seq_along(cor_list)) {
  
  tmp <- cor_list[[i]]
  
  genes_present <- intersect(
    gene_names,
    rownames(tmp)
  )
  
  cor_array[
    genes_present,
    genes_present,
    i
  ] <- tmp[
    genes_present,
    genes_present
  ]
}


# Calculate mean correlation across species

mean_cor <- apply(
  cor_array,
  c(1, 2),
  mean,
  na.rm = TRUE
)

diag(mean_cor) <- NA


# Cluster Core ISGs by mean correlation

mean_cor_cluster <- mean_cor

mean_cor_cluster[
  is.na(mean_cor_cluster)
] <- 0

diag(mean_cor_cluster) <- 1

dist_gene <- as.dist(
  1 - mean_cor_cluster
)

hc_gene <- hclust(
  dist_gene,
  method = "average"
)

gene_order_cluster <- hc_gene$labels[
  hc_gene$order
]


# Prepare correlation matrix for plotting

mean_cor_df <- as.data.frame(
  as.table(mean_cor)
)

colnames(mean_cor_df) <- c(
  "Gene1",
  "Gene2",
  "Correlation"
)

plot_cor_cluster <- mean_cor_df %>%
  mutate(
    Gene1 = factor(
      Gene1,
      levels = gene_order_cluster
    ),
    Gene2 = factor(
      Gene2,
      levels = rev(gene_order_cluster)
    )
  )


# Plot clustered Core ISG correlation matrix
p_gene_cluster <- ggplot(
  plot_cor_cluster,
  aes(
    x = Gene1,
    y = Gene2,
    fill = Correlation
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  coord_equal() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-0.8, 0.8),
    oob = scales::squish,
    name = "Mean\nSpearman's rho"
  ) +
  labs(
    x = "Core ISG",
    y = "Core ISG"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7
    ),
    axis.text.y = element_text(size = 7),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

print(p_gene_cluster)

ggsave(
  "gene_cor_10sp_clust_260930.pdf",
  p_gene_cluster,
  width = 6,
  height = 6
)
