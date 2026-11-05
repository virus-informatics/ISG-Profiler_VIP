rm(list = ls(all = TRUE))
# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(car)
library(effectsize)

# Working directory
setwd(
  "SOMEWHERE"
)
# Read metadata
dat2 <- fread(
  "SRA_metadata_annotated_v2.tsv",
  sep = "\t",
  header = TRUE,
  fill = TRUE
)

# Only columns required for this analysis
dat2_sub <- dat2[, .(
  Run_accession,
  Tissue_normalized,
  Tissue_status
)]

dat2_sub[dat2_sub == ""] <- NA

# Read ISG data
ISG_meta <- read.table(
  "Aves_Mam_mbio_250613_ISGscore_logan.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Keep infection-negative samples
neg_data <- ISG_meta %>%
  filter(
    Infection == "Negative"
  )

# Merge metadata
meta <- merge(
  neg_data,
  dat2_sub,
  by.x = "ID",
  by.y = "Run_accession",
  all.x = TRUE
) %>%
  as_tibble()


# Define QC-approved tissue annotation
meta <- meta %>%
  mutate(
    Tissue_QC = case_when(

      Tissue_status %in% c(
        "Single_match",
        "Single_match_with_context"
      ) &
        !is.na(Tissue_normalized) &
        !Tissue_normalized %in% c(
          "Unknown",
          "Other",
          "Mixed"
        ) ~ as.character(Tissue_normalized),

      TRUE ~ NA_character_
    )
  )


# Dataset for main model
dat_main <- meta %>%
  filter(
    !is.na(ISG_mean),
    !is.na(order),
    !is.na(Tissue_QC)
  ) %>%
  droplevels()

cat(
  "Main model N =",
  nrow(dat_main),
  "\n"
)

cat(
  "Host orders =",
  n_distinct(dat_main$order),
  "\n"
)

cat(
  "Tissues =",
  n_distinct(dat_main$Tissue_QC),
  "\n"
)


# Main linear model
fit_main <- lm(
  ISG_mean ~
    order +
    Tissue_QC,
  data = dat_main
)


# Type II ANOVA
anova_main <- car::Anova(
  fit_main,
  type = 2
)

print(anova_main)


# Partial eta squared
eta_main <- effectsize::eta_squared(
  anova_main,
  partial = TRUE
)

print(eta_main)


# Delta R2
fit_main_no_order <- update(
  fit_main,
  . ~ . - order
)

fit_main_no_tissue <- update(
  fit_main,
  . ~ . - Tissue_QC
)

r2_main <- summary(fit_main)$r.squared

delta_r2_main <- tibble(
  Variable = c(
    "Host order",
    "Tissue"
  ),
  Delta_R2 = c(
    r2_main -
      summary(fit_main_no_order)$r.squared,

    r2_main -
      summary(fit_main_no_tissue)$r.squared
  )
)

print(delta_r2_main)


# Prepare partial eta2 for plotting
plot_eta_main <- as.data.frame(
  eta_main
) %>%
  as_tibble() %>%
  transmute(
    Variable = case_when(
      Parameter == "order" ~ "Host order",
      Parameter == "Tissue_QC" ~ "Tissue",
      TRUE ~ as.character(Parameter)
    ),
    Metric = "Partial \u03b7\u00b2",
    Effect_size = Eta2_partial
  )


# Prepare Delta R2 for plotting

plot_delta_main <- delta_r2_main %>%
  transmute(
    Variable = Variable,
    Metric = "\u0394R\u00b2",
    Effect_size = Delta_R2
  )


# Combine effect sizes
effect_plot_main <- bind_rows(
  plot_eta_main,
  plot_delta_main
) %>%
  mutate(
    Variable = factor(
      Variable,
      levels = c(
        "Host order",
        "Tissue"
      )
    ),
    Metric = factor(
      Metric,
      levels = c(
        "Partial \u03b7\u00b2",
        "\u0394R\u00b2"
      )
    )
  )


# Plot
p_effect_main <- ggplot(
  effect_plot_main,
  aes(
    x = Variable,
    y = Effect_size
  )
) +

  geom_col(
    width = 0.65
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        Effect_size
      )
    ),
    vjust = -0.5,
    size = 4
  ) +

  facet_wrap(
    ~ Metric,
    nrow = 1
  ) +

  scale_y_continuous(
    limits = c(0, 0.22),
    expand = expansion(
      mult = c(0, 0.08)
    )
  ) +

  labs(
    x = NULL,
    y = "Effect size"
  ) +

  theme_classic() +

  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 11
    ),
    axis.text.x = element_text(
      size = 10
    )
  )

print(p_effect_main)
ggsave(
  "effect_size_order_tissue_main_model_260930.pdf",
  p_effect_main,
  width = 6,
  height = 4
)
