library(tidyr)
library(ggplot2)
library(dplyr)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/phase_heat_map/")

# Load heatmap data
#data <- read.csv("phase_map.hiphase.csv", header = TRUE, stringsAsFactors = FALSE)
#data <- read.csv("phase_map.whatshap.csv", header = TRUE, stringsAsFactors = FALSE)
data <- read.csv("phase_map.longphase.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(data) <- gsub("\\.", "-", colnames(data))  # Replace dots with dashes

# Convert samples and genes into factors to preserve order
sample_order <- rev(data$sample)  # Reverse sample order
gene_order <- colnames(data)[-1]  # Keep gene order

df_long <- pivot_longer(
  data,
  cols = -sample,       
  names_to = "Gene",    
  values_to = "Value"   
)

df_long$sample <- factor(df_long$sample, levels = sample_order)
df_long$Gene <- factor(df_long$Gene, levels = gene_order)

# Load incomplete.csv and format labels
#incomplete <- read.csv("incomplete.csv", header = TRUE, stringsAsFactors = FALSE)
#incomplete$label <- paste0(incomplete$num_haploblocks, " (", incomplete$prop_phased, ")")

# Merge `incomplete.csv` with `df_long` while preserving order
#df_long <- df_long %>%
  #left_join(incomplete, by = c("sample", "Gene" = "gene")) %>%
  #mutate(label = ifelse(Value == 0, label, NA))  # Only label missing/0 tiles

# Ensure factor levels remain correct after merge
df_long$sample <- factor(df_long$sample, levels = sample_order)
df_long$Gene <- factor(df_long$Gene, levels = gene_order)

# Plot heatmap
figure <- ggplot(df_long, aes(x = Gene, y = sample, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  #geom_text(aes(label = label), size = 3, na.rm = TRUE) +  # Add text only to missing/0 tiles
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = "Gene",
    y = "Sample"
  )

print(figure)

# Save the figure
#ggsave(filename = "hiphase_heat_map.pdf", plot = figure)
#ggsave(filename = "hiphase_heat_map.png", plot = figure)
#ggsave(filename = "whatshap_heat_map.pdf", plot = figure)
#ggsave(filename = "whatshap_heat_map.png", plot = figure)
#ggsave(filename = "longphase_whatshap_heat_map.pdf", plot = figure)
ggsave(filename = "longphase_heat_map.png", plot = figure)

# Load necessary library
library(dplyr)

# Example data (assuming 'data' is already loaded)
# Convert all non-sample columns to numeric if needed
data_numeric <- data %>% select(-sample)

# Count zeros per row
zero_counts <- rowSums(data_numeric == 0)

# Count occurrences
all_ones <- sum(zero_counts == 0)  # Samples with all 1s
one_zero <- sum(zero_counts == 1)  # Samples with exactly one zero
two_zeros <- sum(zero_counts == 2)  # Samples with exactly two zeros

# Print results
cat("Samples with all 1s:", all_ones, "\n")
cat("Samples with exactly one 0:", one_zero, "\n")
cat("Samples with exactly two 0s:", two_zeros, "\n")



