rm(list = ls(all = TRUE))

# Packages
library(dplyr)
library(tidyr)
library(stringr)
library(broom)
library(ggplot2)
library(patchwork)
library(data.table)


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

ISG_meta <- ISG_meta %>%
  filter(type == "ISG")


# Calculate sample-level ISG score
ISG_mean_df <- ISG_meta %>%
  group_by(ID) %>%
  summarise(
    ISG_mean = mean(ISG_score, na.rm = TRUE),
    species = first(species),
    Infection = first(Infection),
    .groups = "drop"
  ) %>%
  transmute(
    ID,
    hum_symbol = "ISG score",
    type = "ISG score",
    ISG_score = ISG_mean,
    species,
    Infection
  )

ISG_meta <- bind_rows(
  ISG_meta,
  ISG_mean_df
)


# Read geNomad data
Aves_Mam_mbio_genomad_filt <- read.table(
  "Aves_Mam_mbio_genomad_filt_ML_250630.txt",
  sep = "\t",
  header = TRUE
)


# Prepare virus metadata
virus_meta <- Aves_Mam_mbio_genomad_filt %>%
  select(
    ID,
    taxonomy,
    nucleotide
  )

taxonomy_map <- virus_meta %>%
  select(
    taxonomy,
    nucleotide
  ) %>%
  distinct() %>%
  separate(
    taxonomy,
    into = c(
      "Virus",
      "Realm",
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family"
    ),
    sep = ";",
    fill = "right",
    remove = FALSE
  ) %>%
  mutate(
    Family = na_if(Family, ""),
    Family = case_when(
      is.na(Family) & nucleotide == "DNA" ~ "DNA Unclassified",
      is.na(Family) & nucleotide == "RNA" ~ "RNA Unclassified",
      TRUE ~ Family
    )
  ) %>%
  select(
    taxonomy,
    Family,
    nucleotide
  ) %>%
  distinct()


# Join virus annotation
ISG_pos <- ISG_meta %>%
  filter(Infection == "Positive") %>%
  inner_join(
    virus_meta,
    by = "ID",
    relationship = "many-to-many"
  )

ISG_neg <- ISG_meta %>%
  filter(Infection == "Negative") %>%
  mutate(
    taxonomy = "Negative",
    nucleotide = "Negative"
  )

ISG_all <- bind_rows(
  ISG_pos,
  ISG_neg
) %>%
  filter(
    !is.na(ISG_score),
    !is.na(species),
    !is.na(taxonomy)
  )


# Factorize

ISG_all$taxonomy <- relevel(
  factor(ISG_all$taxonomy),
  ref = "Negative"
)

ISG_all$species <- factor(
  ISG_all$species
)

# Split data by ISG

gene_list <- split(
  ISG_all,
  ISG_all$hum_symbol,
  drop = TRUE
)

rm(
  ISG_pos,
  ISG_neg,
  ISG_all
)

gc()


# Estimate viral taxonomy effects
result_file <- "virus_effect_taxonomy_260820_2.tsv"

if (!file.exists(result_file)) {

  results <- vector(
    "list",
    length(gene_list)
  )

  names(results) <- names(gene_list)

  for (i in seq_along(gene_list)) {

    gene <- names(gene_list)[i]

    message(
      "[",
      i,
      "/",
      length(gene_list),
      "] ",
      gene
    )

    dat <- gene_list[[i]]

    if (
      nlevels(droplevels(dat$taxonomy)) < 2 ||
      nlevels(droplevels(dat$species)) < 2
    ) {
      next
    }

    fit <- tryCatch(
      lm(
        ISG_score ~ taxonomy + species,
        data = dat
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) next

    results[[i]] <- tidy(fit) %>%
      filter(
        str_starts(
          term,
          "taxonomy"
        )
      ) %>%
      mutate(
        ISG = gene,
        taxonomy = str_remove(
          term,
          "^taxonomy"
        )
      )
  }

  virus_effect_df <- bind_rows(results)

  fwrite(
    virus_effect_df,
    result_file,
    sep = "\t"
  )

} else {

  virus_effect_df <- fread(
    result_file
  )
}


# Add viral family information

virus_effect_df <- virus_effect_df %>%
  left_join(
    taxonomy_map,
    by = "taxonomy"
  )


# FDR correction

virus_effect_df <- virus_effect_df %>%
  group_by(ISG) %>%
  mutate(
    p_adj = p.adjust(
      p.value,
      method = "BH"
    )
  ) %>%
  ungroup()


# Define viral family order

family_order <- virus_effect_df %>%
  filter(
    ISG == "ISG score",
    !Family %in% c(
      "Negative",
      "DNA Unclassified",
      "RNA Unclassified"
    )
  ) %>%
  arrange(estimate) %>%
  pull(Family)


# Define ISG order

gene_order <- virus_effect_df %>%
  filter(
    ISG != "ISG score"
  ) %>%
  group_by(ISG) %>%
  summarise(
    mean_est = mean(
      estimate,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(mean_est)
  ) %>%
  pull(ISG)

gene_order <- c(
  gene_order,
  "ISG score"
)


# Prepare heatmap data

plot_df <- virus_effect_df %>%
  filter(
    Family %in% family_order
  ) %>%
  mutate(
    Family = factor(
      Family,
      levels = family_order
    ),
    ISG = factor(
      ISG,
      levels = rev(gene_order)
    )
  )


# Plot viral family × ISG effect heatmap

p_heatmap <- ggplot(
  plot_df,
  aes(
    x = Family,
    y = ISG,
    fill = estimate
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Estimate"
  ) +
  labs(
    x = "",
    y = "ISG"
  ) +
  theme_minimal(
    base_size = 10
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(
      size = 7
    ),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )


# Read viral family classification

virus_level <- read.table(
  "low_ISG_virus_260109.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

family_annotation <- data.frame(
  Family = family_order
) %>%
  left_join(
    virus_level %>%
      select(
        Family,
        ISG_level
      ),
    by = "Family"
  ) %>%
  mutate(
    ISG_level = case_when(
      ISG_level == "low" ~ "No ISG-inducing",
      is.na(ISG_level) ~ "ISG-inducing"
    )
  )


# Plot viral family annotation

p_annotation <- ggplot(
  family_annotation,
  aes(
    x = Family,
    y = "Viral family\ncategory",
    fill = ISG_level
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "ISG-inducing" = "#C51B2A",
      "No ISG-inducing" = "grey75"
    ),
    name = "Viral family\ncategory"
  ) +
  scale_x_discrete(
    limits = family_order,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(
    base_size = 10
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8
    ),
    axis.text.y = element_text(
      size = 10
    ),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )


# Combine and save viral family × ISG heatmap

p_final <- p_heatmap / p_annotation +
  plot_layout(
    heights = c(8, 1),
    guides = "collect"
  )

print(p_final)

ggsave(
  "VirusFamily_ISG_heatmap_260930.pdf",
  p_final,
  width = 10,
  height = 8
)


# Prepare viral family effect matrix

effect_df <- virus_effect_df %>%
  filter(
    Family %in% family_order
  ) %>%
  mutate(
    Family = factor(
      Family,
      levels = family_order
    )
  ) %>%
  arrange(Family) %>%
  select(
    Family,
    ISG,
    estimate
  ) %>%
  pivot_wider(
    names_from = ISG,
    values_from = estimate
  )

family_names <- effect_df$Family

effect_mat <- effect_df %>%
  select(-Family) %>%
  as.matrix()

rownames(effect_mat) <- family_names

effect_mat[
  is.na(effect_mat)
] <- 0


# Calculate viral family correlations

cor_mat <- cor(
  t(effect_mat),
  method = "spearman"
)

diag(cor_mat) <- NA

cor_df <- as.data.frame(
  as.table(cor_mat)
)

colnames(cor_df) <- c(
  "Family1",
  "Family2",
  "Correlation"
)

cor_df <- cor_df %>%
  mutate(
    Family1 = factor(
      Family1,
      levels = family_order
    ),
    Family2 = factor(
      Family2,
      levels = family_order
    )
  )


# Plot viral family correlation matrix
p_cor <- ggplot(
  cor_df,
  aes(
    x = Family1,
    y = Family2,
    fill = Correlation
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.55
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-0.8, 0.8),
    oob = scales::squish
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    x = "",
    y = "Family"
  )


# Combine and save viral family correlation matrix
p_final_cor <- p_cor / p_annotation +
  plot_layout(
    heights = c(8, 1),
    guides = "collect"
  )

print(p_final_cor)

ggsave(
  "VirusFamily_cor_matrix_260930.pdf",
  p_final_cor,
  width = 10,
  height = 8
)
