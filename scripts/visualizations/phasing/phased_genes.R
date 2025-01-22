# Load necessary libraries
library(tidyverse)

# Set working directory
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/phasing/")

# File path for the new data
phased_genes_file <- "phased_genes.csv"

# Function to process and clean the new data
process_phased_genes <- function(file_path) {
  # Read the data with correct delimiter
  data <- read_tsv(file_path, col_names = TRUE) %>%  # Use read_tsv for tab-delimited files
    mutate(
      platform = case_when(
        platform == "Revio" ~ "PacBio",
        platform == "PromethION" ~ "ONT",
        TRUE ~ platform
      ),
      sample_group = case_when(
        sample %in% hprc_samples ~ "HPRC",
        sample %in% ihw_samples ~ "IHW",
        TRUE ~ "Other"
      )
    ) %>%
    filter(sample_group != "Other") %>%  # Keep only HPRC and IHW
    drop_na(num_genes)  # Remove rows with NA in num_genes
  
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

# Process the new data
phased_genes_data <- process_phased_genes(phased_genes_file)

# Manually set the x-axis order and labels
phased_genes_data <- phased_genes_data %>%
  mutate(
    x_group = factor(
      paste(platform, sample_group),
      levels = c("PacBio HPRC", "PacBio IHW", "ONT HPRC", "ONT IHW")  # Correct order
    )
  )

# Create the box plot using the processed data
phased_genes_plot <- ggplot(phased_genes_data, aes(x = x_group, y = num_genes, fill = x_group)) +
  geom_boxplot(
    alpha = 0.3, width = 0.3, outlier.shape = NA, color = "black"
  ) +
  geom_jitter(
    aes(color = x_group),
    width = 0.15, height = 0, size = 2, alpha = 0.8  # Adjust jitter width to match box width
  ) +
  scale_y_continuous(labels = scales::label_number()) +  # Adjust y-axis labels to standard numbers
  scale_fill_viridis_d(option = "D", begin = 0.1, end = 0.9) +  # Viridis color for fill
  scale_color_viridis_d(option = "D", begin = 0.1, end = 0.9) +  # Viridis color for jitter
  labs(x = "", y = "Number of Phased\nGenes") +
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
print(phased_genes_plot)

ggsave(filename = "phased_genes.pdf", plot = phased_genes_plot)
ggsave(filename = "phased_genes.png", plot = phased_genes_plot)
