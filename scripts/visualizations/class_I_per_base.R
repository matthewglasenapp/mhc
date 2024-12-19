library(tidyverse)
library(ggplot2)
library(dplyr)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/")

# per-base coverage file
hla_per_base <- "/Users/matt/Documents/GitHub/mhc/clean_data/revio_mean_std_depth.rds"

data <- readRDS(hla_per_base)

# Collapse data into 100-base-pair windows
collapsed_data <- data %>%
  mutate(window = floor(base / 100) * 100 + 50) %>% # Create windows centered at midpoints
  group_by(window) %>%
  summarize(
    mean_depth = mean(mean_depth),
    std_depth = mean(std_depth)
  )

# Rename 'window' to 'base' for consistency
collapsed_data <- collapsed_data %>%
  rename(base = window)

# MHC Class I Part 1
data_filtered <- collapsed_data %>% filter(base >= 29715000 & base <= 30019999)

# MHC Class II
#data_filtered <- data %>% filter(base >= 32438877 & base <= 33144325)

# MHC Class III
#data_filtered <- data %>% filter(base >= 31383288 & base <= 32438877)

figure <- ggplot(data_filtered, aes(x = base, y = mean_depth)) +
  geom_line(size = 0.75) +
  xlab("Position on Chromosome 6 (Mb)") +
  ylab("Mean Coverage Depth") +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_x_continuous(
    labels = function(x) sprintf("%.2f", x / 1e6),
    breaks = pretty(data_filtered$base, n = 10)
  ) +
  scale_y_continuous(
    limits = c(-220, max(data_filtered$mean_depth)),  # Adjusted y-axis limits for 3x lower positioning
    expand = c(0, 0)
  ) +
  geom_rect(
    aes(xmin = 29792233, ymin = -150, xmax = 29793136, ymax = -125),  # Shifted ymin and ymax downward
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29792233, 29793136)), y = -175,  # Text centered below the rectangle
    label = "HLA-V", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29941259, ymin = -150, xmax = 29949572, ymax = -125),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29941259, 29949572)), y = -175,
    label = "HLA-A", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29722774, ymin = -150, xmax = 29738528, ymax = -125),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29722774, 29738528)), y = -175,
    label = "HLA-F", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29800414, ymin = -150, xmax = 29802425, ymax = -125),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29800414, 29802425)), y = -175,
    label = "HLA-P", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29826966, ymin = -150, xmax = 29831125, ymax = -125),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29826966, 29831125)), y = -175,
    label = "HLA-G", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29887751, ymin = -150, xmax = 29890482, ymax = -125),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29887751, 29890482)), y = -175,
    label = "HLA-H", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29926458, ymin = -150, xmax = 29929232, ymax = -125),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29926458, 29929232)), y = -175,
    label = "HLA-K", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29956595, ymin = -150, xmax = 29958570, ymax = -125),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29956595, 29958570)), y = -175,
    label = "HLA-W", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 30006605, ymin = -150, xmax = 30009539, ymax = -125),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(30006605, 30009539)), y = -175,
    label = "HLA-J", color = "black", vjust = 0, size = 3
  )

# Display the figure
print(figure)


ggsave(filename = "HLA_Class_I.pdf", plot = figure, width=169, units = "mm")
#ggsave(filename = "HLA_Class_I.png", plot = figure, width=169, units = "mm")