# Load required libraries
library(ggplot2)
library(dplyr)

# Define maximum value for capping
max_value <- 400

# Step 1: Read and Prepare Data
data <- read.table("pacbio_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract only the sample columns (excluding gene_name)
sample_data <- data[, -1]

# Cap all values larger than max_value
sample_data[sample_data > max_value] <- max_value

# Step 2: Calculate Number of Bins
# Use Sturges' rule based on the number of observations in a single sample
num_bins <- ceiling(log2(length(sample_data[[1]])) + 1)

# Step 3: Define Breaks
# Create breaks for the histogram using the capped range [0, max_value]
breaks <- seq(0, max_value, length.out = num_bins + 1)

# Step 4: Calculate Histogram Counts for Each Sample
# Function to compute histogram counts
calculate_histogram <- function(column, breaks) {
  hist(column, breaks = breaks, plot = FALSE)$counts
}

# Apply the histogram function to each column of the sample data
hist_counts <- apply(sample_data, 2, calculate_histogram, breaks = breaks)

# Step 5: Calculate Mean Counts and Error
mean_counts <- rowMeans(hist_counts, na.rm = TRUE)
sd_counts <- apply(hist_counts, 1, sd, na.rm = TRUE)
n_samples_per_bin <- rowSums(!is.na(hist_counts))
sem_counts <- sd_counts / sqrt(n_samples_per_bin)

# Ensure no negative lower bounds for error bars
ymin <- pmax(mean_counts - sem_counts, 0)
ymax <- mean_counts + sem_counts

# Step 6: Create Bin Labels
bin_labels <- paste0(
  as.integer(head(breaks, -1)), "-", as.integer(tail(breaks, -1) - 1)
)
# Correct the last bin label to ">last_bin_start"
bin_labels[length(bin_labels)] <- paste0(">", as.integer(breaks[length(breaks) - 1]))

# Create a data frame for plotting
plot_data <- data.frame(
  BinLabel = factor(bin_labels, levels = bin_labels), # Ensure levels are ordered
  MeanCount = mean_counts,
  SEM = sem_counts
)

# Step 7: Generate the Plot
figure <- ggplot(plot_data, aes(x = BinLabel, y = MeanCount)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "black") +
  labs(
    title = "Mean Coverage Depth Distribution for Targeted Genes",
    x = "Histogram Bin",
    y = "Mean Count for Bin"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# Display the plot
figure

# Save the figure
ggsave(filename = "histogram_bar_plot.png", plot = figure, width = 169, units = "mm")
ggsave(filename = "histogram_bar_plot.pdf", plot = figure, width = 169, units = "mm")
