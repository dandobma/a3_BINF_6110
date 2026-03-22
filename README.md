# Assignment 3

## Matei Dan-Dobre

### 1407506

# Introduction

The human gut microbiome is a highly diverse and dynamic microbial community that plays a critical role in host physiology, including metabolism, immune regulation, and disease susceptibility (Lynch & Pedersen, 2016). Growing evidence has demonstrated that the composition and function of the gut microbiome are strongly influenced by environmental factors, among which diet is one of the most significant modulators (Sonnenburg & Bäckhed, 2016).

Both short-term dietary interventions and long-term dietary patterns have been shown to alter gut microbial communities. For example, rapid shifts in diet can lead to measurable changes in microbiome composition within days, highlighting the responsiveness of the microbiome to dietary inputs (David et al., 2014). Over longer timescales, habitual dietary patterns such as vegan and omnivorous diets have been associated with distinct microbial signatures. Diets rich in plant-based foods tend to promote the growth of fiber-degrading and short-chain fatty acid–producing taxa, whereas diets higher in animal-derived products are associated with different microbial profiles and metabolic outputs (De Filippis et al., 2019).

Advances in sequencing technologies, particularly shotgun metagenomics, have enabled detailed characterization of microbial communities at high taxonomic resolution. Unlike amplicon-based approaches, shotgun metagenomics allows for species-level identification and functional inference. Tools such as Kraken2 and Bracken have facilitated this process by providing rapid and accurate taxonomic classification and abundance estimation from sequencing reads (Wood et al., 2019; Lu et al., 2017).

Despite these advances, identifying consistent diet-associated differences in the gut microbiome remains challenging due to substantial inter-individual variability and limitations in sample size across studies. While some studies report clear distinctions between dietary groups, others find that individual variation can obscure broader patterns, particularly in smaller cohorts (Lynch & Pedersen, 2016; Sonnenburg & Bäckhed, 2016). This highlights the importance of applying robust analytical pipelines and cautious interpretation when examining microbiome data.

In this context, the present analysis examines publicly available shotgun metagenomic data from the SRA project SRP126540 to explore differences in gut microbiome composition between vegan and omnivorous individuals. By integrating taxonomic classification, diversity analyses, and differential abundance testing, this analysis aims to evaluate whether dietary patterns are associated with detectable differences in microbial communities within a small sample set.

# Methods

Shotgun metagenomic sequencing data were obtained from the NCBI Sequence Read Archive (SRA) under project accession SRP126540, associated with the study by De Filippis et al. (2019). Six samples were selected for analysis, consisting of three vegan and three omnivore gut microbiome samples. Raw sequencing reads were downloaded using the SRA Toolkit (fasterq-dump) and processed as paired-end FASTQ files. Data retrieval was automated using a Bash script (scripts/download_data.sh) to ensure reproducibility.

Taxonomic classification of sequencing reads was performed using Kraken2 (Wood et al., 2019) with a pre-built standard reference database (~8 GB). Kraken2 assigns taxonomic labels to sequencing reads by matching k-mers to a reference database. For each sample, paired-end reads were classified using multithreading, and classification reports were generated for downstream analysis. This step was executed using a custom script (scripts/run_kraken2.sh).

Species-level abundance estimation was performed using Bracken (Lu et al., 2017), which refines Kraken2 classifications by re-estimating taxon abundances based on k-mer distributions. Bracken was applied to Kraken2 output reports using a read length of 150 bp and species-level classification. This step produced corrected abundance estimates for each taxon in each sample and was automated using the script (scripts/run_bracken.sh).

All downstream analyses were conducted in R (version 4.5.2). Bracken output files were imported and combined into a single abundance matrix, where rows corresponded to taxa and columns corresponded to samples. A sample metadata table was constructed to indicate dietary group (vegan or omnivore) for each sample. A phyloseq object was then created to facilitate analysis and visualization (McMurdie & Holmes, 2013). These steps were implemented in an R script (scripts/analysis.R).

Relative abundance data were calculated by normalizing counts within each sample. Taxonomic composition was visualized using stacked bar plots of the most abundant taxa. Alpha diversity was assessed using Shannon and Simpson diversity indices, and differences between dietary groups were evaluated using the Wilcoxon rank-sum test.

Beta diversity was assessed using Bray–Curtis dissimilarity, followed by principal coordinates analysis (PCoA) to visualize differences in community composition between samples. Statistical significance of group differences was evaluated using permutational multivariate analysis of variance (PERMANOVA) implemented via the adonis2 function in the vegan package (Oksanen et al., 2020).

Differential abundance analysis was performed using ANCOMBC2 (Lin & Peddada, 2020), which accounts for compositional bias in microbiome data. The model tested for differences in taxon abundance between vegan and omnivore groups, with p-values adjusted using the Benjamini–Hochberg method. Sensitivity analysis was enabled to assess the robustness of detected differences.


