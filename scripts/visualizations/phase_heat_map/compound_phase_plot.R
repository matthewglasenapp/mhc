library(tidyverse)
library(viridis)
library(patchwork)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/phase_heat_map/")

exclude_drb5 <- TRUE

# =====================
# 1. BOX PLOT: NUM GENES
# =====================
hprc_samples <- c("HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258",
                  "HG01928", "HG02055", "HG02630", "HG03492", "HG03579",
                  "NA19240", "NA20129", "NA21309", "NA24694", "NA24695")
ihw_samples <- c("IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118",
                 "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200",
                 "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364",
                 "IHW09409")

phased_data <- read_csv("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/phasing/phased_genes2.csv") %>%
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
  filter(sample_group != "Other") %>%
  mutate(x_group = factor(paste(platform, sample_group),
                          levels = c("PacBio HPRC", "PacBio IHW", "ONT HPRC", "ONT IHW")))

p1 <- ggplot(phased_data, aes(x = x_group, y = num_genes, fill = x_group)) +
  geom_boxplot(alpha = 0.3, width = 0.3, outlier.shape = NA, color = "black") +
  geom_jitter(aes(color = x_group), width = 0.15, height = 0, size = 1.5, alpha = 0.8) +
  scale_y_continuous(labels = scales::label_number()) +
  scale_fill_viridis_d(option = "D", begin = 0.1, end = 0.9) +
  scale_color_viridis_d(option = "D", begin = 0.1, end = 0.9) +
  labs(x = "", y = "Number of Phased\nGenes") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8, face = "bold"),
        legend.position = "none")

# ==============================
# 2. HEATMAP: PERCENT COVERED (All Samples, Ordered)
# ==============================

# Load coverage matrix
coverage_data <- read_csv("phase_map.hiphase_revio.csv")
colnames(coverage_data) <- gsub("\\.", "-", colnames(coverage_data))

# Preserve ordering
sample_order <- rev(coverage_data$sample)         # top → bottom
gene_order <- colnames(coverage_data)[-1]         # left → right

# Reshape
df_coverage <- coverage_data %>%
  pivot_longer(cols = -sample, names_to = "Gene", values_to = "Value")

# Load % coverage info
incomplete <- read_csv("incomplete.hiphase_revio.csv") %>%
  mutate(largest_haploblock = as.numeric(gsub("%", "", largest_haploblock)))

# Join and compute percent_covered
df_coverage <- df_coverage %>%
  left_join(incomplete, by = c("sample", "Gene" = "gene")) %>%
  mutate(percent_covered = ifelse(Value == 1, 100, largest_haploblock)) %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    Gene = factor(Gene, levels = gene_order)
  )

df_coverage_filtered <- if (exclude_drb5) {
  df_coverage %>% filter(Gene != "HLA-DRB5")
} else {
  df_coverage
}

# Plot
p2 <- ggplot(df_coverage_filtered, aes(x = Gene, y = sample, fill = percent_covered)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "D", direction = -1, na.value = "white",
                       name = "Percent Gene\nCovered by\nLargest\nOverlapping\nHaploblock") +
  labs(x = "Gene", y = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6, angle = 30, hjust = 1, vjust = 1.2),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        panel.grid = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        axis.ticks = element_blank())

# ===========================
# 3. HEATMAP: NORMALIZED EDIT DISTANCE (Per Haplotype, Class I only)
# ===========================

# Define HPRC sample subset and Class I genes
hprc_samples <- c("HG002", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630",
                  "HG03492", "HG03579", "NA19240", "NA20129", "NA21309")
class_I_genes <- c("HLA-A", "HLA-B", "HLA-C")

# Load and prepare edit distance data
edit_data <- read_csv("edit_distance.csv") %>%
  mutate(
    hap1 = ifelse(hap1 == "NA", NA, hap1),
    hap2 = ifelse(hap2 == "NA", NA, hap2),
    alignment_length_1 = as.numeric(alignment_length_1),
    alignment_length_2 = as.numeric(alignment_length_2)
  )

# Join with coverage data to align gene/sample rows
df_edit <- df_coverage %>%
  left_join(edit_data, by = c("sample", "Gene" = "gene"))

# Filter to subset + class I and reshape
df_haplo <- df_edit %>%
  filter(sample %in% hprc_samples, Gene %in% class_I_genes) %>%
  select(sample, Gene, hap1, hap2, alignment_length_1, alignment_length_2) %>%
  pivot_longer(cols = c(hap1, hap2), names_to = "hap", values_to = "edit_distance") %>%
  mutate(
    alignment_length = ifelse(hap == "hap1", alignment_length_1, alignment_length_2),
    edit_distance = as.numeric(edit_distance),
    similarity = ifelse(!is.na(edit_distance) & !is.na(alignment_length),
                        1 - (edit_distance / alignment_length), NA),
    Gene_hap = paste0(Gene, "_", ifelse(hap == "hap1", "1", "2")),
    label = case_when(
      is.na(similarity) ~ NA_character_,
      similarity == 1 ~ "1",
      TRUE ~ sprintf("%.3f", similarity)
    ),
    text_color = ifelse(is.na(similarity) | similarity < 0.9, "black", "white")
  ) %>%
  mutate(
    sample = factor(sample, levels = sample_order),        # preserve vertical order
    Gene_hap = factor(Gene_hap, levels = unique(Gene_hap)) # preserve horizontal order
  )

# Plot
p3 <- ggplot(df_haplo, aes(x = Gene_hap, y = sample, fill = similarity)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label, color = text_color), size = 2.5, fontface = "bold", na.rm = TRUE) +
  scale_fill_viridis_c(option = "D", direction = -1, na.value = "grey", name = "Percent\nSequence\nIdentity") +
  scale_color_identity() +
  labs(x = "Gene Haplotype", y = "Sample") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )
# =====================
# COMBINE AND SAVE
# =====================
compound <- p1 / p2 / p3 + plot_layout(ncol = 1, heights = c(1, 2, 1))
ggsave("compound_plot_fixed.pdf", compound, width = 6.5, height = 9)
ggsave("compound_plot_fixed.png", compound, width = 6.5, height = 9)

print(compound)
