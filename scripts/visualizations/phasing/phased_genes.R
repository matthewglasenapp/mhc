# Load necessary libraries
library(tidyverse)

# Set working directory
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/phasing/")

# File path for the new data
phased_genes_file <- "phased_genes2.csv"

# Function to process and clean the new data
process_phased_genes <- function(file_path) {
  # Read the data with correct delimiter
  data <- read_csv(file_path, col_names = TRUE) %>%  # Use read_tsv for tab-delimited files
    mutate(
      platform = case_when(
        platform == "pacbio" ~ "PacBio",
        platform == "ont" ~ "ONT",
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

#===============================================================================
library(tidyr)
library(ggplot2)
library(dplyr)

# Load main heatmap data
data <- read.csv("phase_map.hiphase.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(data) <- gsub("\\.", "-", colnames(data))  # Replace dots with dashes

# Define HPRC samples and Class I genes
hprc_samples <- c("HG002", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "NA19240", "NA20129", "NA21309")
class_I <- c("HLA-A", "HLA-B", "HLA-C")

# Convert samples and genes into factors to preserve order
sample_order <- rev(data$sample)
gene_order <- colnames(data)[-1]

# Reshape data
df_long <- pivot_longer(data, cols = -sample, names_to = "Gene", values_to = "Value")
df_long$sample <- factor(df_long$sample, levels = sample_order)
df_long$Gene <- factor(df_long$Gene, levels = gene_order)

# ✅ Load incomplete.csv FIRST
incomplete <- read.csv("incomplete.hiphase.csv", header = TRUE, stringsAsFactors = FALSE)

# Merge and calculate percent covered
df_long <- df_long %>%
  left_join(incomplete, by = c("sample", "Gene" = "gene")) %>%
  mutate(
    largest_haploblock = as.numeric(gsub("%", "", largest_haploblock)),
    percent_covered = ifelse(Value == 1, 100, largest_haploblock)
  )

# Ensure order is still preserved after join
df_long$sample <- factor(df_long$sample, levels = sample_order)
df_long$Gene <- factor(df_long$Gene, levels = gene_order)

# Final plot — no annotations, just fill
test4 <- ggplot(df_long, aes(x = Gene, y = sample, fill = percent_covered)) +
  geom_tile(color = "white") +
  #scale_fill_gradient(low = "white", high = "red", na.value = "white",
  #name = "Percent Gene\nCovered by\nLargest\nOverlapping/nHaploblock") 
  scale_fill_viridis_c(option = "D", direction = -1, na.value = "white", 
                       name = "Percent Gene\nCovered by\nLargest\nOverlapping\nHaploblock")
theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Gene",
    y = "Sample"
  )

print(test4)

#===============================================================================

# Define HPRC samples and class I genes
hprc_samples <- c("HG002", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "NA19240", "NA20129", "NA21309")
class_I <- c("HLA-A", "HLA-B", "HLA-C")

# Make long format from df_long with edit distances and alignment lengths
df_classI_long <- df_long %>%
  filter(sample %in% hprc_samples, Gene %in% class_I) %>%
  select(sample, Gene, hap1, hap2, alignment_length_1, alignment_length_2) %>%
  mutate(
    hap1 = ifelse(hap1 == "NA", NA, hap1),
    hap2 = ifelse(hap2 == "NA", NA, hap2),
    alignment_length_1 = as.numeric(alignment_length_1),
    alignment_length_2 = as.numeric(alignment_length_2)
  ) %>%
  pivot_longer(
    cols = c(hap1, hap2),
    names_to = "hap",
    values_to = "edit_distance"
  ) %>%
  mutate(
    alignment_length = ifelse(hap == "hap1", alignment_length_1, alignment_length_2),
    edit_distance = as.numeric(edit_distance),
    similarity_score = ifelse(!is.na(edit_distance) & !is.na(alignment_length),
                              1 - (edit_distance / alignment_length), NA),
    Gene_hap = paste0(Gene, "_", ifelse(hap == "hap1", "1", "2")),
    label = case_when(
      is.na(similarity_score) ~ NA_character_,
      similarity_score == 1 ~ "1",
      TRUE ~ sprintf("%.3f", similarity_score)
    ),
    text_color = ifelse(is.na(similarity_score) | similarity_score < 0.9, "black", "white")
  )

# Plot
similarity_plot <- ggplot(df_classI_long, aes(x = Gene_hap, y = sample, fill = similarity_score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label, color = text_color), size = 4, fontface = "bold", na.rm = TRUE) +
  #scale_fill_gradient(high = "red", low = "white", na.value = "grey", name = #"Percent\nSequence\nIdentity") +
  scale_fill_viridis_c(option = "D", direction = -1, na.value = "grey", name = "Percent\nSequence\nIdentity") + 
  scale_color_identity() +
  labs(
    title = "Normalized Edit Distance to HPRC HLA Haplotypes",
    x = "Gene Haplotype",
    y = "Sample"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(similarity_plot)
