library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/downsample/")
csv_path <- "all_concordance_results.csv"

# 🔧 SELECT MHC CLASS TO PLOT
#selected_class <- "MHC_Class_II"
selected_class <- "MHC_Class_II"
#selected_class <- "MHC_Class_III"

# 🔧 TOGGLE to include 0.5% (0.005) data or not
include_low_proportions <- FALSE

# Read and filter data
df <- read_csv(csv_path) %>%
  filter(MHC_Class == selected_class)

if (!include_low_proportions) {
  df <- df %>% filter(Proportion >= 0.01)
}

# Get rounded mean depth per Proportion
depth_labels <- df %>%
  group_by(Proportion) %>%
  summarise(Depth = round(mean(Depth)), .groups = "drop")

# Summarize metrics
summary_df <- df %>%
  group_by(Proportion, Variant, Metric) %>%
  summarise(
    mean_value = mean(Value),
    sd_value = sd(Value),
    sem = sd(Value) / sqrt(n()),
    .groups = "drop"
  )

# Define color mapping and plot order
metric_colors <- c(
  "Precision" = "grey50",
  "Recall" = "grey70",
  "F1" = "#4169E1"  # Royal blue
)

summary_df$Metric <- factor(summary_df$Metric, levels = c("Precision", "Recall", "F1"))

plot_metric <- function(variant_type, show_legend = TRUE) {
  filtered <- summary_df %>% 
    filter(Variant == variant_type) %>%
    arrange(Metric)  # Ensures F1 is on top
  
  min_y <- min(filtered$mean_value - filtered$sd_value, na.rm = TRUE)
  
  ggplot(filtered, aes(x = Proportion * 100, y = mean_value, color = Metric, shape = Metric)) +
    geom_errorbar(
      aes(ymin = mean_value - sem, ymax = mean_value + sem),
      width = 0.15,
      linewidth = 0.6
    ) +
    geom_point(size = 2) +
    geom_text(
      data = depth_labels,
      aes(x = Proportion * 100, y = min_y - 0.05, label = paste0(Depth, "X")),
      inherit.aes = FALSE,
      color = "black", size = 3, vjust = 1
    ) +
    scale_color_manual(values = metric_colors) +
    scale_shape_manual(values = c(16, 17, 15)) +
    scale_x_continuous(
      breaks = seq(1, 10, 1),
      labels = paste0(seq(1, 10, 1), "%"),
      expand = expansion(mult = c(0.04, 0.02))  # ⬅️ Increased left padding
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.05))) +
    labs(
      x = "Percent Reads Retained",
      y = "Metric Value",
      title = variant_type,
      color = "Metric"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = ifelse(show_legend, "top", "none"),
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.justification = "center",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.background = element_rect(fill = "white", color = NA),  # explicitly set background
      axis.line = element_line(color = "black")
    )
}

# Create plots
# Show SNP and INDEL plots separately, each with its own legend
plot_snp <- plot_metric("SNP", show_legend = TRUE)
plot_indel <- plot_metric("INDEL", show_legend = TRUE)

print(plot_snp)
print(plot_indel)


# Combine vertically
final_plot <- (plot_snp / plot_indel) +
  plot_annotation(title = paste("Genotype Concordance for", gsub("_", " ", selected_class)))

print(final_plot)
print(plot_snp)

ggsave("mhc_III_downsample.pdf", plot = final_plot)
ggsave("mhc_III_snp.pdf", plot = plot_snp)
ggsave("mhc_III_indel.pdf", plot = plot_indel)

# =============================================================================
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/downsample/")
csv_path <- "all_concordance_results.csv"

# TOGGLE
include_low_proportions <- FALSE

# Load and filter
df <- read_csv(csv_path)

if (!include_low_proportions) {
  df <- df %>% filter(Proportion >= 0.01)
}

# Calculate depth labels
depth_labels <- df %>%
  group_by(MHC_Class, Proportion) %>%
  summarise(Depth = round(mean(Depth)), .groups = "drop")

# Summarize metrics
summary_df <- df %>%
  group_by(MHC_Class, Proportion, Variant, Metric) %>%
  summarise(
    mean_value = mean(Value),
    sd_value = sd(Value),
    sem = sd(Value) / sqrt(n()),
    .groups = "drop"
  )

# Color palette
metric_colors <- c(
  "Precision" = "grey50",
  "Recall" = "grey70",
  "F1" = "#4169E1"
)

summary_df$Metric <- factor(summary_df$Metric, levels = c("Precision", "Recall", "F1"))

# Plotting function
plot_metric <- function(class_name, variant_type, show_legend = TRUE) {
  filtered <- summary_df %>%
    filter(MHC_Class == class_name, Variant == variant_type) %>%
    arrange(Metric)
  
  if (nrow(filtered) == 0) {
    warning("No data for: ", class_name, " - ", variant_type)
    return(ggplot() + ggtitle(paste("Missing:", class_name, variant_type)))
  }
  
  min_y <- min(filtered$mean_value - filtered$sd_value, na.rm = TRUE)
  
  # Join depth labels safely
  depth_joined <- filtered %>%
    distinct(MHC_Class, Proportion) %>%
    left_join(depth_labels, by = c("MHC_Class", "Proportion"))
  
  ggplot(filtered, aes(x = Proportion * 100, y = mean_value, color = Metric, shape = Metric)) +
    geom_errorbar(
      aes(ymin = mean_value - sem, ymax = mean_value + sem),
      width = 0.15,
      linewidth = 0.3,
      show.legend = FALSE
    ) + 
    geom_point(size = 0.75) +
    geom_text(
      data = depth_joined,
      aes(x = Proportion * 100, y = min_y - 0.05, label = paste0(Depth, "X")),
      inherit.aes = FALSE,
      color = "black", size = 1.5, vjust = 1
    ) +
    scale_color_manual(values = metric_colors) +
    scale_shape_manual(values = c(16, 17, 15)) +
    scale_x_continuous(
      breaks = seq(1, 10, 1),
      labels = paste0(seq(1, 10, 1), "%"),
      expand = expansion(mult = c(0.04, 0.02))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.05))) +
    labs(
      x = "Percent Reads Retained",
      y = "Metric Value",
      title = paste(variant_type, "-", gsub("_", " ", class_name)),  # ASCII dash only
      color = "Metric"
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      legend.position = ifelse(show_legend, "bottom", "none"),
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 6),       # Smaller text
      legend.justification = "center",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      plot.title = element_text(size = 8),
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6)
    )
}

# Plot generation
classes <- c("MHC_Class_I", "MHC_Class_II", "MHC_Class_III")
plots <- list()

for (cls in classes) {
  message("Generating plots for: ", cls)
  plots[[paste0(cls, "_SNP")]]   <- plot_metric(cls, "SNP", show_legend = (cls == "MHC_Class_I"))
  plots[[paste0(cls, "_INDEL")]] <- plot_metric(cls, "INDEL", show_legend = FALSE)
}

# Combine in 3 rows × 2 columns
final_plot <- (
  wrap_plots(plots, ncol = 2, byrow = TRUE, guides = "collect")
) &
  theme(legend.position = "bottom") &
  guides(
    color = guide_legend(override.aes = list(size = 2)),
    shape = guide_legend(override.aes = list(size = 2))
  )


# Show in RStudio
print(final_plot)

# Export cleanly
ggsave("compound_downsample.pdf", plot = final_plot, width = 6.5, height = 9, units = "in")
ggsave("compound_downsample.png", plot = final_plot, width = 6.5, height = 9, units = "in", dpi = 300)
