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

# Function to compute regression stats and equation
compute_stats <- function(data) {
  model <- lm(promethion_coverage ~ pacbio_coverage, data = data)
  stats <- summary(model)
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  r_squared <- stats$r.squared
  p_value <- coef(summary(model))["pacbio_coverage", "Pr(>|t|)"]
  equation <- paste0("y = ", round(slope, 3), "x + ", round(intercept, 3))
  list(r_squared = r_squared, p_value = p_value, equation = equation)
}

# Add residuals to the dataset
add_residuals <- function(data) {
  model <- lm(promethion_coverage ~ pacbio_coverage, data = data)
  data <- data %>%
    mutate(
      predicted = predict(model, newdata = data),
      residual = promethion_coverage - predicted
    )
  return(data)
}

# Compute stats and add residuals for HPRC and IHW datasets
hprc_stats <- compute_stats(hprc_data)
ihw_stats <- compute_stats(ihw_data)
hprc_data <- add_residuals(hprc_data)
ihw_data <- add_residuals(ihw_data)

# Format p-values properly
format_pval <- function(p_value) {
  if (p_value < .Machine$double.eps) {
    return("< 2.2e-16") # Smallest value R can reliably compute
  } else {
    return(format.pval(p_value, scientific = TRUE, digits = 3))
  }
}

# Create scatter plots with regression equation and residual-based color
plot_hprc <- ggplot(hprc_data, aes(x = pacbio_coverage, y = promethion_coverage)) +
  geom_point(aes(color = abs(residual)), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  annotate(
    "text",
    x = 0.1 * max(hprc_data$pacbio_coverage, na.rm = TRUE),
    y = 0.85 * max(hprc_data$promethion_coverage, na.rm = TRUE),
    label = paste0("R² = ", round(hprc_stats$r_squared, 2), "\n",
                   "p = ", format_pval(hprc_stats$p_value), "\n",
                   hprc_stats$equation),
    hjust = 0, color = "black", size = 4
  ) +
  scale_color_viridis_c(name = "Deviation") +
  labs(
    title = "A. HPRC",
    x = "PacBio Coverage",
    y = "Promethion\nCoverage"
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

plot_ihw <- ggplot(ihw_data, aes(x = pacbio_coverage, y = promethion_coverage)) +
  geom_point(aes(color = abs(residual)), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  annotate(
    "text",
    x = 0.1 * max(ihw_data$pacbio_coverage, na.rm = TRUE),
    y = 0.85 * max(ihw_data$promethion_coverage, na.rm = TRUE),
    label = paste0("R² = ", round(ihw_stats$r_squared, 2), "\n",
                   "p = ", format_pval(ihw_stats$p_value), "\n",
                   ihw_stats$equation),
    hjust = 0, color = "black", size = 4
  ) +
  scale_color_viridis_c(name = "Deviation") +
  labs(
    title = "B. IHW",
    x = "PacBio Coverage",
    y = "Promethion\nCoverage"
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

# Combine the two plots using patchwork
combined_plot <- plot_hprc + plot_ihw + plot_layout(ncol = 1)

# Print the combined plot
print(combined_plot)



ggsave(filename = "dotplot.pdf", plot = combined_plot, width=169, units = "mm")
ggsave(filename = "dotplot.png", plot = combined_plot, width=169, units = "mm")




