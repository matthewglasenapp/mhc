# Load necessary libraries
library(tidyverse)
library(ggplot2)

# Define sample groups
hprc_samples <- c("HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258",
                  "HG01928", "HG02055", "HG02630", "HG03492", "HG03579",
                  "NA19240", "NA20129", "NA21309", "NA24694", "NA24695")

ihw_samples <- c("IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118",
                 "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200",
                 "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364",
                 "IHW09409")

# File paths
revio_file <- "revio_flagstat_results.csv"
promethion_file <- "promethion_flagstat_results.csv"

# Function to process and clean data (no duplicate adjustment)
process_data <- function(file_path, platform_name) {
  cat("\nDEBUG: Reading data for platform:", platform_name, "\n")
  data <- read_csv(file_path, col_names = TRUE) %>%
    mutate(
      total = as.numeric(gsub(",", "", total)),  # Remove commas and convert to numeric
      unique_reads = total,  # No adjustment for duplicates
      platform = platform_name  # Add platform information
    ) %>%
    select(sample, unique_reads, platform) %>%
    drop_na()  # Remove any rows with NA in sample
  
  cat("\nDEBUG: Processed data (first 5 rows):\n")
  print(head(data, 5))
  return(data)
}

# Process Revio and PromethION data
revio_data <- process_data(revio_file, "Revio")
promethion_data <- process_data(promethion_file, "PromethION")

# Combine datasets
combined_data <- bind_rows(revio_data, promethion_data)

# Debug: Print combined data
cat("\nDEBUG: Combined data (first 10 rows):\n")
print(head(combined_data, 10))

# Calculate total reads per platform
platform_totals <- combined_data %>%
  group_by(platform) %>%
  summarise(total_reads = sum(unique_reads))

# Debug: Print total reads per platform
cat("\nDEBUG: Total reads per platform:\n")
print(platform_totals)

# Add percentage of total reads for each sample
combined_data <- combined_data %>%
  left_join(platform_totals, by = "platform") %>%
  mutate(percent_of_total = (unique_reads / total_reads) * 100)

# Debug: Check percentages
cat("\nDEBUG: Percentages (first 10 rows):\n")
print(head(combined_data, 10))

# Use the `sample_order` from the plot that accounts for duplicates
# Ensure this `sample_order` matches the duplicate-accounted plot
sample_order <- c(
  "HG002", "NA20129", "NA24694", "HG01258", "HG01106", "HG003", "NA21309", 
  "HG01928", "HG02055", "HG004", "HG02630", "HG03492", "HG005", "NA24695", 
  "NA19240", "HG03579", "IHW09122", "IHW09359", "IHW09117", "IHW09049", 
  "IHW09125", "IHW09224", "IHW09200", "IHW09198", "IHW09409", "IHW09251", 
  "IHW09021", "IHW09175", "IHW09071", "IHW09118", "IHW09364", "IHW09245"
)

# Update sample factor levels to enforce custom order
combined_data <- combined_data %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    platform = factor(platform, levels = c("Revio", "PromethION"))  # Ensure Revio is first
  ) %>%
  filter(!is.na(sample))  # Remove any rows with NA in sample

# Debug: Check levels of sample factor
cat("\nDEBUG: Levels of sample factor:\n")
print(levels(combined_data$sample))

# Define colorblind-friendly palette
custom_colors <- c("Revio" = "#7B3294", "PromethION" = "#FFC20F")  # Purple and Yellow

# Update the plot
combined_plot <- ggplot(combined_data, aes(x = sample, y = percent_of_total, fill = platform)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 1, width = 0.8) +  # Adjust dodge and bar width
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "Sample",
    y = "Percentage of Total Reads",  # Wrap Y-axis title
    fill = "Platform",
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.ticks.x = element_line(linewidth = 0.3, color = "grey"),  # Add ticks for x-axis
    axis.ticks.length.x = unit(0.2, "cm"),  # Make x-axis ticks noticeable
    panel.grid.major.y = element_line(linewidth = 0.3, color = "lightgrey"),  # Subtle gridlines
    panel.grid.major.x = element_blank(),  # No vertical gridlines
    panel.background = element_rect(fill = "white", color = NA),  # Ensure clean background
    legend.position = "top",
    plot.title = element_text(size = 14, face = "bold")
  )

# Display the updated plot
print(combined_plot)

# Save the updated plot
ggsave(filename = "mapping4.pdf", plot = combined_plot, width = 169, units = "mm")
ggsave(filename = "mapping4.png", plot = combined_plot, width = 169, units = "mm")

