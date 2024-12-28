# Load ggplot2 library
library(ggplot2)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/histogram/")


# Define the values
values <- c(125, 108, 104, 100, 73, 113, 136, 86, 107, 91, 98, 102, 74, 85, 100, 84, 94, 104, 122, 111, 98, 88, 109, 92, 97, 96, 126, 112, 134, 108, 81, 86, 137, 90, 88, 54, 111, 132, 134, 88, 101, 106, 74, 56, 92, 125, 104, 133, 125, 83, 47, 46, 163, 119, 87, 64, 61, 60, 75, 95, 122, 122, 127, 127, 132, 94, 107, 91, 113)


# Create a data frame for ggplot2
data <- data.frame(values = values)

# Plot the histogram
figure <- ggplot(data, aes(x = values)) +
  geom_histogram(breaks = seq(0, max(values) + 20, by = 20), # Define exact breaks
                 fill = "skyblue", color = "black") + 
  labs(x = "Mean Gene-Wide Coverage Depth\nTargeted Class III Genes\nAveraged Across Samples", y = "Frequency") +
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) + 
  geom_vline(xintercept = 40, linetype = "dashed", color = "red", size = 1) +
  annotate(
    "text",
    x = 42,                      # Slightly to the right of the line for visibility
    y = 15,
    label = "X = 40",
    hjust = 0,                   # Left-align the text
    color = "red",
    size = 4
  )

figure

ggsave(filename = "histogram2.pdf", plot = figure, width=169, units = "mm")
ggsave(filename = "histogram2.png", plot = figure, width=169, units = "mm")


