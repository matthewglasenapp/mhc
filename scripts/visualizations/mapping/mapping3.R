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

# Function to process and clean data
process_data <- function(file_path, platform_name) {
  cat("\nDEBUG: Reading data for platform:", platform_name, "\n")
  data <- read_csv(file_path, col_names = TRUE) %>%
    mutate(
      total = as.numeric(gsub(",", "", total)),  # Remove commas and convert to numeric
      duplicates = as.numeric(gsub(",", "", duplicates)),  # Remove commas and convert to numeric
      unique_reads = total - duplicates,  # Calculate unique reads
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

# Calculate total unique reads per platform
platform_totals <- combined_data %>%
  group_by(platform) %>%
  summarise(total_unique_reads = sum(unique_reads))

# Debug: Print total reads per platform
cat("\nDEBUG: Total reads per platform:\n")
print(platform_totals)

# Add percentage of total unique reads for each sample
combined_data <- combined_data %>%
  left_join(platform_totals, by = "platform") %>%
  mutate(percent_of_total = (unique_reads / total_unique_reads) * 100)

# Debug: Check percentages
cat("\nDEBUG: Percentages (first 10 rows):\n")
print(head(combined_data, 10))

# Separate data into HPRC and IHW groups
combined_data <- combined_data %>%
  mutate(group = case_when(
    sample %in% hprc_samples ~ "HPRC",
    sample %in% ihw_samples ~ "IHW",
    TRUE ~ "Unknown"
  ))

# Debug: Check group assignment
cat("\nDEBUG: Group assignment counts:\n")
print(table(combined_data$group))

# Sort samples within each group (HPRC/IHW) by Revio values in descending order
hprc_order <- combined_data %>%
  filter(group == "HPRC" & platform == "Revio") %>%
  arrange(desc(unique_reads)) %>%
  pull(sample)

ihw_order <- combined_data %>%
  filter(group == "IHW" & platform == "Revio") %>%
  arrange(desc(unique_reads)) %>%
  pull(sample)

# Debug: Print ordered samples
cat("\nDEBUG: Ordered samples for HPRC:\n")
print(hprc_order)
cat("\nDEBUG: Ordered samples for IHW:\n")
print(ihw_order)

# Combine the custom sample order
sample_order <- c(hprc_order, ihw_order)

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
custom_colors <- c("Revio" = "#7B3294", "PromethION" = "#ffc20F")  # Muted Purple and Bright Yellow/Gold


# Update the plot with axis ticks for better readability
combined_plot <- ggplot(combined_data, aes(x = sample, y = percent_of_total, fill = platform)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 1, width = 0.8) +  # Adjust dodge and bar width
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "Sample",
    y = "Percentage of Total Unique Reads",  # Wrap Y-axis title
    fill = "Platform",
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),  # Smaller font for Y-axis title
    axis.ticks.x = element_line(linewidth = 0.3, color = "grey"),  # Add ticks for x-axis
    axis.ticks.length.x = unit(0.2, "cm"),  # Make x-axis ticks noticeable
    panel.grid.major.y = element_line(linewidth = 0.3, color = "lightgrey"),  # Subtle gridlines
    panel.grid.major.x = element_blank(),  # No vertical gridlines
    panel.background = element_rect(fill = "white", color = NA),  # Ensure clean background
    legend.position = "top",
  )

# Display the updated plot
print(combined_plot)

# Save the updated plot
ggsave(filename = "mapping3.pdf", plot = combined_plot, width = 169, units = "mm")
ggsave(filename = "mapping3.png", plot = combined_plot, width = 169, units = "mm")




