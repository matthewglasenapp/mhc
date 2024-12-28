# Load required libraries
library(ggplot2)
library(dplyr)
library(patchwork)

# Set working directory
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/histogram/")

# Define maximum value for capping
max_value <- 400

# Read and Prepare Data for PacBio
pacbio_data <- read.table("pacbio_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pacbio_sample_data <- pacbio_data[, -1]
pacbio_sample_data[pacbio_sample_data > max_value] <- max_value

# Read and Prepare Data for PromethION
promethion_data <- read.csv("promethion_genes.csv", header = TRUE, stringsAsFactors = FALSE)
promethion_sample_data <- promethion_data[, -1]
promethion_sample_data[promethion_sample_data > max_value] <- max_value

# Define HPRC and IHW samples
hprc_samples <- c("HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258",
                  "HG01928", "HG02055", "HG02630", "HG03492", "HG03579",
                  "NA19240", "NA20129", "NA21309", "NA24694", "NA24695")
ihw_samples <- c("IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118",
                 "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200",
                 "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364",
                 "IHW09409")

# Subset data for HPRC and IHW samples
pacbio_hprc_data <- pacbio_sample_data[, colnames(pacbio_sample_data) %in% hprc_samples]
pacbio_ihw_data <- pacbio_sample_data[, colnames(pacbio_sample_data) %in% ihw_samples]
promethion_hprc_data <- promethion_sample_data[, colnames(promethion_sample_data) %in% hprc_samples]
promethion_ihw_data <- promethion_sample_data[, colnames(promethion_sample_data) %in% ihw_samples]

# Calculate Number of Bins
num_bins <- ceiling(log2(length(pacbio_hprc_data[[1]])) + 1)
breaks <- seq(0, max_value, length.out = num_bins + 1)

# Function to Calculate Histogram Counts
calculate_histogram <- function(column, breaks) {
  hist(column, breaks = breaks, plot = FALSE)$counts
}

# Function to Create Plot
create_plot <- function(sample_group_data, group_name) {
  hist_counts <- apply(sample_group_data, 2, calculate_histogram, breaks = breaks)
  mean_counts <- rowMeans(hist_counts, na.rm = TRUE)
  sd_counts <- apply(hist_counts, 1, sd, na.rm = TRUE)
  n_samples_per_bin <- rowSums(!is.na(hist_counts))
  sem_counts <- sd_counts / sqrt(n_samples_per_bin)
  
  ymin <- pmax(mean_counts - sem_counts, 0)
  ymax <- mean_counts + sem_counts
  
  bin_labels <- paste0(
    as.integer(head(breaks, -1)), "-", as.integer(tail(breaks, -1) - 1)
  )
  bin_labels[length(bin_labels)] <- paste0(">", as.integer(breaks[length(breaks) - 1]))
  
  plot_data <- data.frame(
    BinLabel = factor(bin_labels, levels = bin_labels),
    MeanCount = mean_counts,
    SEM = sem_counts
  )
  
  ggplot(plot_data, aes(x = BinLabel, y = MeanCount)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "black") +
    labs(
      title = group_name,
      x = "Histogram Bin",
      y = "Mean Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, hjust = 0),      # Reduced plot title size
      axis.title.x = element_text(size = 9),               # Reduced x-axis title size
      axis.title.y = element_text(size = 9),               # Reduced y-axis title size
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # Reduced x-axis tick label size
      axis.text.y = element_text(size = 8),                # Reduced y-axis tick label size
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

# Generate Plots
pacbio_hprc_plot <- create_plot(pacbio_hprc_data, "A. PacBio HPRC")
pacbio_ihw_plot <- create_plot(pacbio_ihw_data, "B. PacBio IHW")
promethion_hprc_plot <- create_plot(promethion_hprc_data, "C. PromethION HPRC")
promethion_ihw_plot <- create_plot(promethion_ihw_data, "D. PromethION IHW")

# Combine Plots using Patchwork
combined_plot <- (pacbio_hprc_plot | pacbio_ihw_plot) / (promethion_hprc_plot | promethion_ihw_plot)

# Display the Combined Plot
print(combined_plot)

# Save the Combined Plot
ggsave("histogram4.png", combined_plot, width = 169, units = "mm")
ggsave("histogram4.pdf", combined_plot, width = 169, units = "mm")
