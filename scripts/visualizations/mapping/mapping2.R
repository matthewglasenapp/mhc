# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(viridis)  # For color palette

# Set working directory
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/mapping/")

# File paths
revio_file <- "revio_flagstat_results.csv"
promethion_file <- "promethion_flagstat_results.csv"

# Function to process and clean data
process_data <- function(file_path, platform) {
  # Read data
  data <- read_csv(file_path, col_names = TRUE) %>%
    mutate(
      total = as.numeric(gsub(",", "", total)),  # Remove commas and convert to numeric
      duplicates = as.numeric(gsub(",", "", duplicates)),  # Remove commas and convert to numeric
      unique_reads = total - duplicates,  # Calculate unique reads
      platform = platform  # Add platform column
    ) %>%
    select(sample, unique_reads, platform) %>%  # Keep necessary columns
    mutate(
      sample_group = case_when(
        sample %in% hprc_samples ~ "HPRC",
        sample %in% ihw_samples ~ "IHW",
        TRUE ~ "Other"
      )
    ) %>%
    filter(sample_group != "Other") %>%  # Keep only HPRC and IHW
    drop_na(unique_reads)  # Remove rows with NA in unique_reads
  
  return(data)
}

# Define sample groups
hprc_samples <- c("HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", 
                  "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", 
                  "NA19240", "NA20129", "NA21309", "NA24694", "NA24695")
ihw_samples <- c("IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", 
                 "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", 
                 "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", 
                 "IHW09409")

# Process data
revio_data <- process_data(revio_file, "PacBio")
promethion_data <- process_data(promethion_file, "ONT")

# Combine datasets
combined_data <- bind_rows(revio_data, promethion_data) %>%
  mutate(
    x_group = paste(platform, sample_group),  # Create combined group for x-axis
    x_group = factor(
      x_group,
      levels = c("PacBio HPRC", "PacBio IHW", "ONT HPRC", "ONT IHW")  # Ensure correct x-axis order
    )
  )

# Create the box plot using viridis colors
combined_plot <- ggplot(combined_data, aes(x = x_group, y = unique_reads, fill = x_group)) +
  geom_boxplot(
    alpha = 0.3, width = 0.3, outlier.shape = NA, color = "black"
  ) +
  geom_jitter(
    aes(color = x_group), 
    width = 0.15, height = 0, size = 2, alpha = 0.8  # Adjust jitter width to match box width
  ) +
  scale_y_continuous(labels = scales::scientific) +
  scale_fill_viridis_d(option = "D", begin = 0.1, end = 0.9) +  # Viridis color for fill
  scale_color_viridis_d(option = "D", begin = 0.1, end = 0.9) +  # Viridis color for jitter
  labs(x = "", y = "Unique\nReads") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", lineheight = 0.8),
    axis.line.y = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    panel.grid.major.y = element_line(color = "gray80", size = 0.5),  # Feint grid lines
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Display the plot
print(combined_plot)


ggsave(filename = "mapping_box.pdf", plot = combined_plot, width=169, units = "mm")
ggsave(filename = "mapping_box.png", plot = combined_plot, width=169, units = "mm")
