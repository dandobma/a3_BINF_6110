####Setup####
#install.packages("BiocManager")
#install.packages(c("remotes", "CVXR"))
#BiocManager::install("ANCOMBC", ask = FALSE, update = TRUE)
#BiocManager::install("microbiome")
# Load libraries
library(dplyr)
library(readr)
library(tidyr)
library(phyloseq)
library(ggplot2)
library(vegan)
library(ANCOMBC)
library(microbiome)


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

####Create metadata####
metadata <- data.frame(
  sample_id = c(
    "SRR8146935", "SRR8146936", "SRR8146938",
    "SRR8146951", "SRR8146952", "SRR8146954"
  ),
  diet = c(
    "omnivore", "omnivore", "omnivore",
    "vegan", "vegan", "vegan"
  )
)

rownames(metadata) <- metadata$sample_id
metadata$sample_id <- NULL

write.csv(metadata, "data/processed/sample_metadata.csv")

# Check metadata to abundance matrix
all(colnames(abundance_mat) == rownames(metadata))

####Build phyloseq object####
otu <- otu_table(as.matrix(abundance_mat), taxa_are_rows = TRUE)
samples <- sample_data(metadata)

ps <- phyloseq(otu, samples)

# Verify
ps

####Transform to relative abundance####
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

####Melt phyloseq object to long format####
plot_df <- psmelt(ps_rel)

# Check columns
colnames(plot_df)
head(plot_df)

####Plot stacked bar chart####
# Keep only top 20 most abundant species
top_taxa <- plot_df %>%
  group_by(OTU) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:20) %>%
  pull(OTU)

plot_top <- plot_df %>%
  mutate(Taxon = ifelse(OTU %in% top_taxa, OTU, "Other"))

# Create bar chart
ggplot(plot_top, aes(x = Sample, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ diet, scales = "free_x", space = "free_x") +
  labs(
    title = "Top 20 species by relative abundance",
    x = "Sample",
    y = "Relative abundance"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save figure
ggsave(
  filename = "results/species_relative_abundance_barplot.png",
  width = 12,
  height = 6,
  dpi = 300
)

####Alpha diversity####
alpha_df <- estimate_richness(ps, measures = c("Shannon", "Simpson"))

# Add metadata
alpha_df$diet <- metadata[rownames(alpha_df), "diet"]

alpha_df

# Plot Shannon diversity
ggplot(alpha_df, aes(x = diet, y = Shannon, fill = diet)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_bw() +
  labs(
    title = "Shannon Diversity by Diet",
    x = "Diet",
    y = "Shannon Diversity"
  )

# Save figure
ggsave(
  filename = "results/Shannon_diversity_by_diet.png",
  width = 12,
  height = 6,
  dpi = 300
)

# Test output
wilcox.test(Shannon ~ diet, data = alpha_df)

# Plot Simpson diversity
ggplot(alpha_df, aes(x = diet, y = Simpson, fill = diet)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_bw() +
  labs(
    title = "Simpson diversity by diet",
    x = "Diet",
    y = "Simpson diversity"
  )

# Save figure
ggsave(
  filename = "results/Simpson_diversity_by_diet.png",
  width = 12,
  height = 6,
  dpi = 300
)

# Test output
wilcox.test(Simpson ~ diet, data = alpha_df)

####Beta diversity####
# Calculate distance
dist_mat <- distance

# Ordination (PCoA)
ord_df <- as.data.frame(ord$vectors)
ord_df$Sample <- rownames(ord_df)
ord_df$diet <- metadata[ord_df$Sample, "diet"]
ord_df

# Plot
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = diet)) +
  geom_point(size = 4) +
  theme_bw() +
  labs(
    title = "PCoA (Bray-Curtis) by diet",
    x = "PCoA1",
    y = "PCoA2"
  )

# Save
ggsave(
  filename = "results/PCoA_by_diet.png",
  width = 12,
  height = 6,
  dpi = 300
)

# PERMANOVA
sample_df <- data.frame(sample_data(ps))

adonis2(dist_mat ~ diet, data = sample_df)

####ANCOMBC####
set.seed(123)

ancom_out <- ancombc2(
  data = ps,
  fix_formula = "diet",
  rand_formula = NULL,
  p_adj_method = "BH",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "diet",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 1,
  verbose = FALSE,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(tol = 1e-2, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100)
)

# Extract results
res_df <- ancom_out$res

head(res_df)
colnames(res_df)

# Get significant taxa
sig <- res_df[res_df$diff_dietvegan == TRUE, ]
sig