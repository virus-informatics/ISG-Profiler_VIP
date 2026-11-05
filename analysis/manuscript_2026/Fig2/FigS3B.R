rm(list = ls(all = TRUE))
# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(data.table)
library(scales)

# Working directory / input
setwd(
  "SOMEWHERE"
)

metadata_file <- "SRA_metadata_annotated_v2.tsv"

isg_file <-
  "Aves_Mam_mbio_250613_ISGscore_logan.txt"


# Parameters
MIN_N <- 5
TOP_N <- 10

order_levels <- c(
  "Anseriformes",
  "Galliformes",
  "Apodiformes",
  "Columbiformes",
  "Passeriformes",
  "Dasyuromorphia",
  "Didelphimorphia",
  "Primates",
  "Lagomorpha",
  "Rodentia",
  "Eulipotyphla",
  "Chiroptera",
  "Artiodactyla",
  "Carnivora",
  "Perissodactyla"
)


# Read metadata
dat2 <- fread(
  metadata_file,
  sep = "\t",
  header = TRUE,
  fill = TRUE
)

dat2_sub <- dat2[, .(
  Run_accession,
  Tissue_normalized,
  Tissue_status
)]

dat2_sub[dat2_sub == ""] <- NA


# Read ISG data and keep infection-negative samples
ISG_meta <- read.table(
  isg_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

neg_data <- ISG_meta %>%
  filter(Infection == "Negative")

# Merge and restrict to 15 manuscript orders
meta <- merge(
  neg_data,
  dat2_sub,
  by.x = "ID",
  by.y = "Run_accession",
  all.x = TRUE
) %>%
  as_tibble() %>%
  filter(
    order %in% order_levels,
    !is.na(ISG_mean)
  ) %>%
  mutate(
    order = factor(
      order,
      levels = order_levels
    )
  )

# Valid tissue annotations
tissue_dat <- meta %>%
  filter(
    Tissue_status %in% c(
      "Single_match",
      "Single_match_with_context"
    ),
    !is.na(Tissue_normalized),
    !Tissue_normalized %in% c(
      "Unknown",
      "Other",
      "Mixed"
    )
  )

# Define top 10 tissues
major_tissues <- tissue_dat %>%
  count(
    Tissue_normalized,
    sort = TRUE
  ) %>%
  slice_head(n = TOP_N) %>%
  pull(Tissue_normalized)

# Overall basal ISG by order
order_isg <- meta %>%
  group_by(order) %>%
  summarise(
    N = n(),
    Median_ISG = median(
      ISG_mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Overall basal ISG by tissue
tissue_overall_isg <- tissue_dat %>%
  filter(
    Tissue_normalized %in% major_tissues
  ) %>%
  group_by(Tissue_normalized) %>%
  summarise(
    N = n(),
    Median_ISG = median(
      ISG_mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Basal ISG by order × tissue
tissue_isg <- tissue_dat %>%
  filter(
    Tissue_normalized %in% major_tissues
  ) %>%
  group_by(
    order,
    Tissue_normalized
  ) %>%
  summarise(
    N = n(),
    Median_ISG = median(
      ISG_mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  complete(
    order = factor(
      order_levels,
      levels = order_levels
    ),
    Tissue_normalized = major_tissues
  )

# Build combined table
# Main matrix = order × tissue
# Top row    = overall by order
# Left col   = overall by tissue
heat_main <- tissue_isg %>%
  transmute(
    x = as.character(order),
    y = as.character(Tissue_normalized),
    Median_ISG = Median_ISG,
    N = N
  )

heat_order <- order_isg %>%
  transmute(
    x = as.character(order),
    y = "Overall",
    Median_ISG = Median_ISG,
    N = N
  )

heat_tissue <- tissue_overall_isg %>%
  transmute(
    x = "Overall",
    y = as.character(Tissue_normalized),
    Median_ISG = Median_ISG,
    N = N
  )

heat_corner <- tibble(
  x = "Overall",
  y = "Overall",
  Median_ISG = NA_real_,
  N = NA_real_
)

# Combine and set axis order
heat_all <- bind_rows(
  heat_main,
  heat_order,
  heat_tissue,
  heat_corner
) %>%
  mutate(
    x = factor(
      x,
      levels = c(
        "Overall",
        order_levels
      )
    ),
    y = factor(
      y,
      levels = c(
        rev(major_tissues),
        "Overall"
      )
    )
  )


# N >= 5 only
# Bubble size = log10(N)
heat_all_plot <- heat_all %>%
  mutate(
    N_plot = if_else(
      !is.na(Median_ISG) & N >= MIN_N,
      log10(N),
      NA_real_
    ),
    Median_ISG_plot = if_else(
      !is.na(N) & N >= MIN_N,
      Median_ISG,
      NA_real_
    )
  )

p_tissue_combined_more5 <- ggplot(
  heat_all_plot,
  aes(
    x = x,
    y = y
  )
) +

  geom_point(
    aes(
      fill = Median_ISG_plot,
      size = N_plot
    ),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +

  scale_fill_gradientn(
    colours = c(
      "#313695",
      "#FFFFBF",
      "#F46D43",
      "#A50026"
    ),
    values = scales::rescale(
      c(
        -0.5,
        0,
        0.25,
        0.5
      )
    ),
    limits = c(
      -0.5,
      0.5
    ),
    oob = scales::squish,
    na.value = "grey90"
  ) +

  scale_size_continuous(
    range = c(2, 12),
    name = "log10(N)"
  ) +

  scale_x_discrete(
    drop = FALSE
  ) +

  scale_y_discrete(
    drop = FALSE
  ) +

  labs(
    x = "Host order",
    y = "Tissue",
    fill = "Median basal\nISG score"
  ) +

  theme_minimal() +

  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(
      size = 10
    ),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p_tissue_combined_more5)
ggsave(
  "basal_ISG_order_tissue_combined_bubble_260930.pdf",
  p_tissue_combined_more5,
  width = 11,
  height = 7
)
