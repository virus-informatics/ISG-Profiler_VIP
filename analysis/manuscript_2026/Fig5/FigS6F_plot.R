# Fig. S6F: Viral genera enriched in the FN fraction
# Plot the genus-level summary prepared by `FigS6F.ipynb`.

library(tidyverse)
library(patchwork)
library(RColorBrewer)

# Define paths
basedir <- "/path/to/project"
output_dir <- file.path(basedir, "output/revise_final/FigS6F")
input_dir <- output_dir
input_file <- file.path(input_dir, "infection_event_prop_ictv_virus_host_db.tsv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

top_n <- 50

proportion_cols <- c("TP", "FN", "TN", "FP")

# Select genera by FN enrichment
top_genera <- read_tsv(input_file, show_col_types = FALSE) %>%
  slice_max(FN_minus_TP, n = top_n, with_ties = FALSE) %>%
  arrange(desc(FN_minus_TP)) %>%
  mutate(
    genus = factor(genus, levels = rev(genus)),
    unclassified_genus = if_else(
      str_starts(as.character(genus), "unclassified"),
      "Yes",
      "No"
    ),
    vertebrate = factor(
      case_when(
        tolower(as.character(is_vertebrate)) == "true" ~ "Yes",
        tolower(as.character(is_vertebrate)) == "false" ~ "No",
        TRUE ~ "Unknown"
      ),
      levels = c("Yes", "No", "Unknown")
    ),
    Aves_Mammalia = factor(
      replace_na(Aves_Mammalia, "Unknown"),
      levels = c("Yes", "No", "Unknown")
    )
  )

taxon_labels <- setNames(
  paste(top_genera$family, top_genera$genus, sep = " | "),
  top_genera$genus
)

# Prepare host annotations
host_annotation <- top_genera %>%
  select(
    genus,
    Vertebrate = vertebrate,
    Aves_Mammalia
  ) %>%
  pivot_longer(
    -genus,
    names_to = "annotation",
    values_to = "status"
  ) %>%
  mutate(annotation = factor(
    annotation,
    levels = c("Vertebrate", "Aves_Mammalia")
  ))

# Prepare heatmap data
heatmap_data <- top_genera %>%
  pivot_longer(
    all_of(proportion_cols),
    names_to = "category",
    values_to = "proportion"
  ) %>%
  mutate(category = factor(category, levels = proportion_cols))

# Set colors and plot theme
rmax <- max(heatmap_data$proportion, na.rm = TRUE)
cmax <- max(top_genera$total_count, na.rm = TRUE)
pal_proportion <- c("white", brewer.pal(9, "BuGn")[-1])
pal_count <- c("white", brewer.pal(9, "YlGnBu")[-1])

base_theme <- theme_classic(base_size = 11) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.4),
    axis.line.y.left = element_line(color = "black", linewidth = 0.4),
    plot.margin = margin(0, 0, 0, 0)
  )

# Plot unclassified-genus annotation
p_unclassified <- ggplot(
  top_genera,
  aes("Unclassified\ngenus", genus, fill = unclassified_genus)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_y_discrete(labels = taxon_labels) +
  scale_fill_manual(
    values = c(Yes = "grey60", No = "white"),
    breaks = "Yes",
    name = "Unclassified genus"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks = element_blank()
  )

# Plot host annotations
p_host <- ggplot(host_annotation, aes(annotation, genus, fill = status)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_x_discrete(labels = c(
    Vertebrate = "Vertebrate\n(ICTV)",
    Aves_Mammalia = "Aves/Mammalia\n(virushostdb)"
  )) +
  scale_fill_manual(
    values = c(Yes = "#87cefa", No = "#ffb6c1", Unknown = "grey80"),
    name = "Host status"
  ) +
  base_theme +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks = element_blank()
  )

# Plot infection-event counts
p_total <- ggplot(top_genera, aes("Total n", genus, fill = total_count)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colours = pal_count, limits = c(0, cmax), na.value = "white",
    name = "Total n"
  ) +
  base_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# Plot infection-event proportions
p_heatmap <- ggplot(
  heatmap_data,
  aes(category, genus, fill = proportion)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colours = pal_proportion, limits = c(0, rmax), na.value = "white",
    name = "Proportion (%)"
  ) +
  base_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# Combine panels
plot <- p_unclassified + p_host + p_total + p_heatmap +
  plot_layout(
    widths = c(1, 2, 1, length(proportion_cols)),
    guides = "collect"
  ) &
  theme(legend.position = "right")

# Save figures
ggsave(
  file.path(output_dir, "FN_top50_genus_heatmap.pdf"),
  plot, width = 8, height = 6.5
)
ggsave(
  file.path(output_dir, "FN_top50_genus_heatmap.png"),
  plot, width = 8, height = 6.5, dpi = 300
)

