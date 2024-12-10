# Load required library
library(ggplot2)

# Set the working directory
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/histogram/")

# Read the data without setting row.names
data <- read.table("pacbio_genes2.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Convert the data to a long format
long_data <- as.numeric(unlist(data[, -1])) # Skip the first column if it contains non-numerical data

# Remove non-finite values
long_data <- long_data[is.finite(long_data)]

# Create a data frame for ggplot2
plot_data <- data.frame(values = long_data)

# Count how many values fall in the first bin (0 <= x < 30)
first_bin_count <- sum(plot_data$values >= 0 & plot_data$values < 30)

# Create a histogram with enhancements
figure <- ggplot(plot_data, aes(x = values)) +
  geom_histogram(binwidth = 30, fill = "skyblue", color = "black", boundary = 0) + 
  labs(x = "Mean Gene-Wide Coverage Depth\nTargeted Class III Genes, All Samples", y = "Frequency") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) + 
  # Add a vertical red dotted line at x = 30
  geom_vline(xintercept = 30, linetype = "dashed", color = "red", size = 1) +
  # Add a label for the first bin
  annotate(
    "text",
    x = 15,                     # Center of the first bin (0 to 30)
    y = first_bin_count + 10,   # Slightly above the bin height
    label = paste(first_bin_count),
    color = "black",
    size = 5
  ) +
  # Add a label for the vertical line
  annotate(
    "text",
    x = 32,                      # Slightly to the right of the line for visibility
    y = 245,
    label = "X = 30",
    hjust = 0,                   # Left-align the text
    color = "red",
    size = 4
  )

figure

ggsave(filename = "histogram.png", plot = figure, width=169, units = "mm")

