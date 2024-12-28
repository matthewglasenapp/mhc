# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(broom) # For extracting p-values
library(viridis) # For color scale

# Set working directory
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/dotplot/")

# Read the CSV files
pacbio_data <- read.csv("pacbio_genes.csv")
promethion_data <- read.csv("promethion_genes.csv")

# Reshape the data to long format
pacbio_long <- pacbio_data %>%
  pivot_longer(
    cols = where(is.numeric), # Pivot only numeric columns (coverage values)
    names_to = "sample",
    values_to = "pacbio_coverage"
  )

promethion_long <- promethion_data %>%
  pivot_longer(
    cols = where(is.numeric), # Pivot only numeric columns (coverage values)
    names_to = "sample",
    values_to = "promethion_coverage"
  )

# Merge the datasets
merged_data <- pacbio_long %>%
  inner_join(promethion_long, by = c("gene_name", "sample"))

# Define HPRC and IHW samples
hprc_samples <- c("HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", 
                  "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", 
                  "NA19240", "NA20129", "NA21309", "NA24694", "NA24695")
ihw_samples <- c("IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", 
                 "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", 
                 "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", 
                 "IHW09409")

# Filter the data for HPRC and IHW samples
hprc_data <- merged_data %>%
  filter(sample %in% hprc_samples)

ihw_data <- merged_data %>%
  filter(sample %in% ihw_samples)

# Calculate mean coverage for each gene across all samples
calculate_means <- function(data) {
  data %>%
    group_by(gene_name) %>%
    summarize(
      mean_pacbio_coverage = mean(pacbio_coverage, na.rm = TRUE),
      mean_promethion_coverage = mean(promethion_coverage, na.rm = TRUE)
    )
}

# Compute mean coverage for HPRC and IHW datasets
hprc_means <- calculate_means(hprc_data)
ihw_means <- calculate_means(ihw_data)

# Function to compute R-squared and p-value for means
compute_stats_means <- function(data) {
  model <- lm(mean_promethion_coverage ~ mean_pacbio_coverage, data = data)
  stats <- summary(model)
  r_squared <- stats$r.squared
  p_value <- coef(summary(model))["mean_pacbio_coverage", "Pr(>|t|)"] # Extract p-value for slope
  list(r_squared = r_squared, p_value = p_value)
}

# Compute R-squared and p-value for the means
hprc_stats <- compute_stats_means(hprc_means)
ihw_stats <- compute_stats_means(ihw_means)

# Create scatter plots for means
plot_hprc <- ggplot(hprc_means, aes(x = mean_pacbio_coverage, y = mean_promethion_coverage)) +
  geom_point(alpha = 0.7, color = "blue") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  annotate(
    "text",
    x = 0.1 * max(hprc_means$mean_pacbio_coverage, na.rm = TRUE),
    y = 0.9 * max(hprc_means$mean_promethion_coverage, na.rm = TRUE),
    label = paste0("R² = ", round(hprc_stats$r_squared, 2), "\n",
                   "p = ", format.pval(hprc_stats$p_value, scientific = TRUE, digits = 3)),
    hjust = 0, color = "black", size = 4
  ) +
  labs(
    title = "Mean Coverage Comparison: HPRC Samples",
    x = "Mean PacBio Coverage",
    y = "Mean PromethION Coverage"
  ) +
  theme_minimal()

plot_ihw <- ggplot(ihw_means, aes(x = mean_pacbio_coverage, y = mean_promethion_coverage)) +
  geom_point(alpha = 0.7, color = "red") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  annotate(
    "text",
    x = 0.1 * max(ihw_means$mean_pacbio_coverage, na.rm = TRUE),
    y = 0.9 * max(ihw_means$mean_promethion_coverage, na.rm = TRUE),
    label = paste0("R² = ", round(ihw_stats$r_squared, 2), "\n",
                   "p = ", format.pval(ihw_stats$p_value, scientific = TRUE, digits = 3)),
    hjust = 0, color = "black", size = 4
  ) +
  labs(
    title = "Mean Coverage Comparison: IHW Samples",
    x = "Mean PacBio Coverage",
    y = "Mean PromethION Coverage"
  ) +
  theme_minimal()

# Combine the two plots using patchwork
combined_plot <- plot_hprc + plot_ihw + plot_layout(ncol = 1)

# Print the combined plot
print(combined_plot)
