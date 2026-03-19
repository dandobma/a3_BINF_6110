####Setup####
# Load libraries
library(dplyr)
library(readr)
library(tidyr)

# Set data directory to bracken output location
data_dir <- "data/processed/bracken"

####Read bracken files####
files <- list.files(data_dir, pattern = "_bracken.txt", full.names = TRUE)

# Read and store in a list
bracken_list <- lapply(files, function(file) {
  df <- read_tsv(file)
  
  # Extract sample name from filename
  sample_name <- gsub("_bracken.txt", "", basename(file))
  
  df$sample <- sample_name
  return(df)
})

# Combine all into one dataframe
bracken_df <- bind_rows(bracken_list)

head(bracken_df)

####Create abundance matrix####
# rows - species, columns = samples, values = counts
abundance_table <- bracken_df %>%
  select(name, sample, new_est_reads) %>%
  pivot_wider(names_from = sample, values_from = new_est_reads, values_fill = 0)

# Check structure
dim(abundance_table)
head(abundance_table)

# Set rownames
abundance_mat <- as.data.frame(abundance_table)
rownames(abundance_mat) <- abundance_mat$name
abundance_mat$name <- NULL

write.csv(abundance_mat, "data/processed/species_abundance_matrix.csv")
