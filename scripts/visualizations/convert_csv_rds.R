# Read the CSV file
#data <- read.csv("/Users/matt/Documents/GitHub/mhc/clean_data/revio_mean_std_depth.tsv", sep = "\t", stringsAsFactors = FALSE)

# Read the CSV file
data <- read.csv("/Users/matt/Documents/GitHub/mhc/clean_data/promethion_mean_std_depth.tsv", sep = "\t", stringsAsFactors = FALSE)

# Save the data to an RDS file
#saveRDS(data, "/Users/matt/Documents/GitHub/mhc/clean_data/revio_mean_std_depth.rds")

# Save the data to an RDS file
saveRDS(data, "/Users/matt/Documents/GitHub/mhc/clean_data/promethion_mean_std_depth.rds")
