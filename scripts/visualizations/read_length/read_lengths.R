
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(broom) # For regression stats

# Set directories
read_length_dir <- "/Users/matt/Documents/GitHub/mhc/scripts/visualizations/read_length/"
dotplot_dir <- "/Users/matt/Documents/GitHub/mhc/scripts/visualizations/dotplot/"

# Define sample groups
hprc_samples <- c("HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258",
                  "HG01928", "HG02055", "HG02630", "HG03492", "HG03579",
                  "NA19240", "NA20129", "NA21309", "NA24694", "NA24695")
ihw_samples <- c("IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118",
                 "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200",
                 "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364",
                 "IHW09409")

# Process box plot data
setwd(read_length_dir)
results_data <- read.csv("results2.csv")

process_boxplot_data <- function(data, variable, rename_to) {
  data %>%
    rename(!!rename_to := !!sym(variable)) %>%
    mutate(sample_group = case_when(
      Sample %in% hprc_samples & Platform == "Revio" ~ "PacBio\nHPRC",
      Sample %in% ihw_samples & Platform == "Revio" ~ "PacBio\nIHW",
      Sample %in% hprc_samples & Platform == "promethion" ~ "ONT\nHPRC",
      Sample %in% ihw_samples & Platform == "promethion" ~ "ONT\nIHW",
      TRUE ~ "Other"
    )) %>%
    filter(sample_group != "Other") %>%
    drop_na(!!sym(rename_to)) %>%
    mutate(sample_group = factor(sample_group, levels = c("PacBio\nHPRC", "PacBio\nIHW", "ONT\nHPRC", "ONT\nIHW")))
}

read_length_data <- process_boxplot_data(results_data, "Median.Read.Length..BED.Overlapping.Reads.", "median_read_length")
depth_data <- process_boxplot_data(results_data, "Average.On.Target.Depth", "avg_on_target_depth")
proportion_data <- process_boxplot_data(results_data, "Proportion.On.Target", "proportion_on_target")

# Create box plots
create_boxplot <- function(data, y_var, y_label, y_scale = scales::comma) {
  ggplot(data, aes(x = sample_group, y = !!sym(y_var), fill = sample_group)) +
    geom_boxplot(alpha = 0.3, width = 0.5, outlier.shape = NA, color = "black") +
    geom_jitter(aes(color = sample_group), width = 0.15, size = 0.75, alpha = 0.9) +
    scale_fill_viridis_d(option = "D", begin = 0.1, end = 0.9) +
    scale_color_viridis_d(option = "D", begin = 0.1, end = 0.9) +
    scale_y_continuous(labels = y_scale) +
    labs(x = "", y = y_label) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8, face = "bold", color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title.y = element_text(size = 8, face = "bold", lineheight = 0.8),
      legend.position = "none"
    )
}

median_read_length_plot <- create_boxplot(read_length_data, "median_read_length", "Median Read Length")
proportion_on_target_plot <- create_boxplot(
  proportion_data, 
  "Enrichment", 
  "Fold Enrichment"
)
avg_on_target_depth_plot <- create_boxplot(depth_data, "avg_on_target_depth", "Average On-Target Depth")

# Process scatter plot data
setwd(dotplot_dir)
pacbio_data <- read.csv("pacbio_genes.csv")
promethion_data <- read.csv("promethion_genes.csv")

pacbio_long <- pacbio_data %>%
  pivot_longer(
    cols = starts_with("HG") | starts_with("IHW"), 
    names_to = "sample", 
    values_to = "pacbio_coverage"
  )

promethion_long <- promethion_data %>%
  pivot_longer(
    cols = starts_with("HG") | starts_with("IHW"), 
    names_to = "sample", 
    values_to = "promethion_coverage"
  )

merged_data <- pacbio_long %>%
  inner_join(promethion_long, by = c("gene_name", "sample"))

# Split merged data into HPRC and IHW
hprc_data <- merged_data %>% filter(sample %in% hprc_samples)
ihw_data <- merged_data %>% filter(sample %in% ihw_samples)

# Add residuals and compute regression stats
compute_stats <- function(data) {
  model <- lm(promethion_coverage ~ pacbio_coverage, data = data)
  stats <- summary(model)
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  r_squared <- stats$r.squared
  p_value <- coef(summary(model))["pacbio_coverage", "Pr(>|t|)"]
  equation <- paste0("y = ", round(slope, 3), "x + ", round(intercept, 3))
  list(model = model, r_squared = r_squared, p_value = p_value, equation = equation)
}

add_residuals <- function(data, stats) {
  data %>%
    mutate(
      predicted = predict(stats$model, newdata = data),
      residual = promethion_coverage - predicted
    )
}

hprc_stats <- compute_stats(hprc_data)
ihw_stats <- compute_stats(ihw_data)

hprc_data <- add_residuals(hprc_data, hprc_stats)
ihw_data <- add_residuals(ihw_data, ihw_stats)

# Format p-values
format_pval <- function(p_value) {
  if (p_value < .Machine$double.eps) {
    "< 2.2e-16"
  } else {
    format.pval(p_value, scientific = TRUE, digits = 3)
  }
}

# Create scatter plots
scatter_plot_hprc <- ggplot(hprc_data, aes(x = pacbio_coverage, y = promethion_coverage)) +
  geom_point(aes(color = abs(residual)), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  annotate(
    "text",
    x = 0.1 * max(hprc_data$pacbio_coverage, na.rm = TRUE),
    y = 0.85 * max(hprc_data$promethion_coverage, na.rm = TRUE),
    label = paste0("R² = ", round(hprc_stats$r_squared, 2), "\n",
                   "p = ", format_pval(hprc_stats$p_value), "\n",
                   hprc_stats$equation),
    hjust = 0, color = "black", size = 3
  ) +
  scale_color_viridis_c(name = "Deviation") +
  labs(title = "HPRC", x = "PacBio Coverage", y = "Promethion Coverage") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

scatter_plot_ihw <- ggplot(ihw_data, aes(x = pacbio_coverage, y = promethion_coverage)) +
  geom_point(aes(color = abs(residual)), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  annotate(
    "text",
    x = 0.1 * max(ihw_data$pacbio_coverage, na.rm = TRUE),
    y = 0.85 * max(ihw_data$promethion_coverage, na.rm = TRUE),
    label = paste0("R² = ", round(ihw_stats$r_squared, 2), "\n",
                   "p = ", format_pval(ihw_stats$p_value), "\n",
                   ihw_stats$equation),
    hjust = 0, color = "black", size = 3
  ) +
  scale_color_viridis_c(name = "Deviation") +
  labs(title = "IHW", x = "PacBio Coverage", y = "Promethion Coverage") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

# Combine all plots
combined_plot <- (
  (median_read_length_plot / proportion_on_target_plot / avg_on_target_depth_plot) | 
    (scatter_plot_hprc / scatter_plot_ihw)
)

# Display the combined plot
print(combined_plot)

# Save the figure
ggsave(filename = "figure_1.pdf", plot = combined_plot)
ggsave(filename = "figure_1.png", plot = combined_plot)

