# Load ggplot2 library
library(ggplot2)

# Define the values
values <- c(248.125, 282.875, 257.21875, 204.21875, 143.65625, 230.5625, 288.40625, 
            166.65625, 195.90625, 156.71875, 169.09375, 206.875, 149.625, 171.65625, 
            193.4375, 170.21875, 196.90625, 219.90625, 262.3125, 225.625, 181.46875, 
            185.375, 239.34375, 183.78125, 198.1875, 208.34375, 268.40625, 236.5625, 
            281.40625, 215.28125, 144.09375, 156.03125, 283, 167.5, 221.625, 80.28125, 
            193.78125, 258.4375, 294.75, 170.15625, 191.53125, 207.03125, 133.84375, 
            104.375, 183.71875, 260.625, 218.28125, 269.28125, 272.1875, 178.4375, 
            122.6875, 123.8125, 368.15625, 241.125, 146.78125, 100.15625, 103.40625, 
            92.03125, 128.21875, 179.625, 245.09375, 262.71875, 279.625, 264.25, 
            276.71875, 192.84375, 240.5625, 174.21875, 232.96875)

# Create a data frame for ggplot2
data <- data.frame(values = values)

# Plot the histogram
figure <- ggplot(data, aes(x = values)) +
  geom_histogram(breaks = seq(0, max(values) + 30, by = 30), # Define exact breaks
                 fill = "skyblue", color = "black") + 
  labs(x = "Mean Gene-Wide Coverage Depth\nTargeted Class III Genes\nAveraged Across Samples", y = "Frequency") +
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) + 
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1) +
  annotate(
    "text",
    x = 62,                      # Slightly to the right of the line for visibility
    y = 15,
    label = "X = 60",
    hjust = 0,                   # Left-align the text
    color = "red",
    size = 4
  )

figure

ggsave(filename = "histogram2.pdf", plot = figure, width=169, units = "mm")


