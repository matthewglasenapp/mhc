library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)  # Load patchwork for side-by-side plots

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/downsample/")

# Load the CSV file
#data <- read.csv("MHC_Class_I_downsample_concordance.csv")
#data <- read.csv("MHC_Class_II_downsample_concordance.csv")
data <- read.csv("MHC_Class_III_downsample_concordance.csv")

# Function to create the plot for a given variant type
plot_variant <- function(variant_to_plot, show_legend = TRUE) {
  # Filter data for the selected Variant
  filtered_data <- subset(data, Variant == variant_to_plot)
  
  # Get unique Proportion-Depth pairs for annotation
  depth_labels <- filtered_data %>%
    group_by(Proportion) %>%
    summarise(Depth = round(unique(Depth)))  # ✅ Round Depth to nearest integer
  
  # Generate the plot
  p <- ggplot(filtered_data, aes(x = Proportion, y = Value, color = Metric, shape = Metric, group = Metric)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_color_viridis_d(option = "viridis") +  
    scale_x_continuous(
      breaks = seq(0.01, 0.1, 0.01),
      labels = function(x) paste0(x * 100, "%")
    ) + 
    theme_minimal() +
    labs(y = "Metric Value", x = "Percent Reads Retained", title = variant_to_plot) +
    
    # ✅ Add Depth annotations below x-axis (rounded)
    geom_text(data = depth_labels, aes(x = Proportion, y = min(filtered_data$Value) - 0.05, 
                                       label = paste0(Depth, "X")), 
              color = "black", size = 3, vjust = 0.9, inherit.aes = FALSE) +
    
    theme(
      legend.position = ifelse(show_legend, "right", "none"),  # ✅ Hide legend for SNP panel
      legend.title = element_blank(),
      text = element_text(size = 14),
      
      # Add axis lines
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.line = element_line(color = "black", size = 0.8)
    )
  
  return(p)
}

# Create plots for SNP and INDEL (Hide legend for SNP to avoid redundancy)
plot_snp <- plot_variant("SNP", show_legend = TRUE)
plot_indel <- plot_variant("INDEL", show_legend = TRUE)

# Combine them side by side with patchwork & ADD OVERALL TITLE
final_plot <- (plot_snp / plot_indel) + 
  plot_annotation(title = "HLA Class III Genotype Concordance")

# Display the final combined plot
print(final_plot)


#ggsave("mhc_I_downsample.pdf", plot = final_plot, width = 8, height = 10)
#ggsave("mhc_I_downsample.png", plot = final_plot, width = 8, height = 10, dpi = 300)

#ggsave("mhc_II_downsample.pdf", plot = final_plot, width = 8, height = 10)
#ggsave("mhc_II_downsample.png", plot = final_plot, width = 8, height = 10, dpi = 300)

ggsave("mhc_III_downsample.pdf", plot = final_plot, width = 8, height = 10)
ggsave("mhc_III_downsample.png", plot = final_plot, width = 8, height = 10, dpi = 300)



