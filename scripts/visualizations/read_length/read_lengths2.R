# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(readr)
library(broom)  # For regression stats

# Set directories
read_length_dir <- "/Users/matt/Documents/GitHub/mhc/scripts/visualizations/read_length/"
dotplot_dir <- "/Users/matt/Documents/GitHub/mhc/scripts/visualizations/dotplot/"
mapping_dir <- "/Users/matt/Documents/GitHub/mhc/scripts/visualizations/mapping/"

# Define sample groups
hprc_samples <- c("HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258",
                  "HG01928", "HG02055", "HG02630", "HG03492", "HG03579",
                  "NA19240", "NA20129", "NA21309", "NA24694", "NA24695")
ihw_samples <- c("IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118",
                 "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200",
                 "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364",
                 "IHW09409")

# Function to process data for boxplots
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

# Process data for existing box plots
setwd(read_length_dir)
results_data <- read.csv("results2.csv")
read_length_data <- process_boxplot_data(results_data, "Median.Read.Length..BED.Overlapping.Reads.", "median_read_length")
depth_data <- process_boxplot_data(results_data, "Average.On.Target.Depth", "avg_on_target_depth")
enrichment_data <- process_boxplot_data(results_data, "Enrichment", "fold_enrichment")

# Process data for the new "Unique Reads" box plot
setwd(mapping_dir)
revio_file <- "revio_flagstat_results.csv"
promethion_file <- "promethion_flagstat_results.csv"

revio_data <- read_csv(revio_file) %>%
  mutate(
    unique_reads = as.numeric(gsub(",", "", total)) - as.numeric(gsub(",", "", duplicates)),
    Sample = sample,  # Rename to match expected column name
    Platform = "Revio"
  )

promethion_data <- read_csv(promethion_file) %>%
  mutate(
    unique_reads = as.numeric(gsub(",", "", total)) - as.numeric(gsub(",", "", duplicates)),
    Sample = sample,  # Rename to match expected column name
    Platform = "promethion"
  )

unique_reads_data <- bind_rows(
  process_boxplot_data(revio_data, "unique_reads", "unique_reads"),
  process_boxplot_data(promethion_data, "unique_reads", "unique_reads")
)

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
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 8, face = "bold", color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title.y = element_text(size = 8, face = "bold", lineheight = 0.8),
      legend.position = "none"
    )
}

median_read_length_plot <- create_boxplot(read_length_data, "median_read_length", "Median Read Length")
avg_on_target_depth_plot <- create_boxplot(depth_data, "avg_on_target_depth", "Average On-Target Depth")
fold_enrichment_plot <- create_boxplot(enrichment_data, "fold_enrichment", "Fold Enrichment")
unique_reads_plot <- create_boxplot(unique_reads_data, "unique_reads", "Unique Reads")

# Combine all plots
compound_plot <- (
  (unique_reads_plot / avg_on_target_depth_plot / scatter_plot_hprc) | ( median_read_length_plot/ fold_enrichment_plot / scatter_plot_ihw)
)

# Display the combined plot
print(compound_plot)

# Save the figure
ggsave(filename = "/Users/matt/Documents/GitHub/mhc/scripts/visualizations/read_length/figure_1.pdf", plot = compound_plot)
ggsave(filename = "/Users/matt/Documents/GitHub/mhc/scripts/visualizations/read_length/figure_1.png", plot = compound_plot)


