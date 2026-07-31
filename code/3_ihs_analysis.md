# iHS analysis

This page details the generation of integrated haplotype scores and subsequent enrichment analyses, which are presented in Figures 3 and 4; Supplementary Figures 6 and 7; and Supplementary Tables 3, 4, and 5.

## iHS per-population

**Step 1) Starting with the Level 4 VCF, for each population generate the IMPUTE format `.hap` file and `.map` files and calculate standardized iHS values using hapbin. Reference 1kG genetic maps were concatenated across all chromosomes.**

`hapbin_ihs.sbatch ${pop}` → `${pop}.ihs`

```bash
#!/bin/bash

source activate vcftools

# Subset Level 4 VCF to each population
bcftools view -S ${pop}.txt -Oz -o ${pop}.vcf.gz
bcftools index ${pop}.vcf.gz

# Create IMPUTE format files
vcftools --gzvcf ${pop}.vcf.gz --IMPUTE --out ${pop}
plink --vcf ${pop}.vcf.gz --cm-map project/lbarreiro/USERS/bridget/references/maps/impute2-1000GP_Phase3/genetic_map_chr@_combined_b37.txt --recode --out ${pop}

# iHS calculated for all variants with a minor allele frequency greater than 5% and up to the point at which EHH drops below 0.1
ihsbin --hap ${pop}.impute.hap --map ${pop}.map --out ${pop}_iHS --minmaf 0.05 --cutoff 0.1
```

**Step 2) Generate iHS stats per-gene by averaging absolute iHS scores within 100kb windows centered on each gene.**

`new_iHS_per_gene_all_pops_BCpipeline.R` → `new_iHS_per_gene_all_pops_lowconf_anc_annot_BCpipeline_gff3_df.RData`

```r
#########################################################################################
#####                                   READ ME                                     #####
#########################################################################################

# computation of standardised iHS values per gene using genomic center - set window size.
# run this script for all populations.

# this script should parallelize looping through populations using foreach

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet_nokin_SHAPEIT5noref_anc_lowconf/low_conf_anc_annot_iHS_inputs")
library(data.table)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)

output_dir <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet_nokin_SHAPEIT5noref_anc_lowconf/low_conf_anc_annot_iHS_inputs"
files <- list.files(path = output_dir, pattern = "\\_iHS$",full.names = TRUE)

# set up temporary data frame to load all final contents
all_pops_df <- data.frame()

# read each file within the range
for (i in 1:length(files)) {
  currentFile <- fread(files[i], col.names = c("index", "rsID", "freq", "iHH_0", "iHH_1", "iHS", "std_iHS"),
                       stringsAsFactors = TRUE)
  currentFileName <- gsub(".*/([A-Za-z]+)_.*", "\\1", files[i])
  currentFile$pop <- c(currentFileName)

  # append the current file to all_pops_df
  all_pops_df <- rbind(all_pops_df, currentFile)
}

all_pops_df <- na.omit(all_pops_df)

######################################################################################
#####        Format all_pops_df with the fake coordinates from rsIDs             #####
######################################################################################
# Extract the chromosome number
# This regex captures digits following 'chr' and before 'pos'
all_pops_df$chr <- as.numeric(str_extract(all_pops_df$rsID, "(?<=chr)\\d+"))
all_pops_df$chr <- as.numeric(all_pops_df$chr)

# Extract the position number
# This regex captures digits following 'pos'
all_pops_df$pos <- as.numeric(str_extract(all_pops_df$rsID, "(?<=pos)\\d+"))
all_pops_df$pos <- as.numeric(all_pops_df$pos)

all_pops_df <- all_pops_df %>% arrange(pop, index)
all_pops_df <- all_pops_df[, c(1,9,10,3:8)]

# check structure
head(all_pops_df)
str(all_pops_df)

#########################################################################################
#####    Get dataframe of genes and their transcription start and end coordinates   #####
#########################################################################################
gff3 <- read.delim("/project/lbarreiro/USERS/bridget/references/GRCh37_ENSEMBL_annotations/gff3/Homo_sapiens.GRCh37.87.chr.gff3.genes_only", header=FALSE, stringsAsFactors=TRUE)
colnames(gff3) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gff3 <- gff3[,c("chr", "source", "type", "start", "end", "strand", "attributes")] #remove irrelevant columns

# extract gene ID and name from each row of gff3

gene_IDs <- gsub(".*ID=gene:(.+);Name.*", "\\1", gff3$attributes)
gene_names <- gsub(".*Name=(.+);biotype.*", "\\1", gff3$attributes)
gff3$gene_IDs <- gene_IDs
gff3$names <- gene_names

# filter out MT, X, Y chromosomes
gff3 <- gff3 %>% filter(!(chr== "MT" | chr == "X" | chr =="Y")) %>% droplevels()

#########################################################################################
#####             Establishing 100kb window from center of each gene                #####
#########################################################################################
gff3$length <- gff3$end - gff3$start # take a look at gene length

# for window size of 100 kb
window_size <- 100000 #100kb

# set lower and upper window boundaries for where iHS would be averaged
gff3$genomic_center <- round((gff3$start + gff3$end) / 2)
gff3$window_lower <- gff3$genomic_center - (window_size / 2)
gff3$window_upper <- gff3$genomic_center + (window_size / 2)

# check structure
head(gff3)
str(gff3)

#########################################################################################
#####                 Compute iHS per gene across populations                       #####
#########################################################################################

# function for calculating average |iHS| within window
compute_ihs <- function(x, pop_df){ #notice pop_df is in there as an argument. this is important
  chrom <- as.numeric(as.character(x[1]))
  lower_bound <- as.numeric(as.character(x[12])) # lower bound position of the window
  upper_bound <- as.numeric(as.character(x[13])) # upper bound position of the window
  snps_for_each_gene <- pop_df[pop_df$chr == chrom & pop_df$pos >= lower_bound & pop_df$pos <= upper_bound, ]
  ihs_for_gene <- mean(abs(snps_for_each_gene$std_iHS))
  return(ihs_for_gene)
}


# set number of cores -----------------------------------------------------
workers <- 6
cl <- parallel::makeCluster(workers)
# register the cluster for using foreach
doParallel::registerDoParallel(cl)

# loop through populations and apply function within for-loop
unique_pops <- unique(all_pops_df$pop)
gff3_ihs <- foreach (i=1:length(unique_pops), .combine = cbind, .packages=c('data.table', 'dplyr', 'foreach', 'doParallel')) %dopar% {
  pop_i <- unique_pops[i]
  pop_df <- all_pops_df %>% filter(pop == pop_i)
  print(unique(pop_df$pop))
  as.vector(apply(gff3, 1, compute_ihs, pop_df))
}

save.image(file="iHS_per_gene_all_pops_lowconf_anc_annot_BCpipeline.RData")

# gff3_ihs has the per-population results
gff3_results <- cbind(gff3, gff3_ihs)
names(gff3_results)[14:47] <- unique_pops

save(gff3_results, file = "new_iHS_per_gene_all_pops_lowconf_anc_annot_BCpipeline_gff3_df.RData")
```

## Gene set enrichment 

**Step 1) Test the enrichment of VIPs (`viral_effects_file_april2023.txt`) for iHS scores > 99% empirical percentile compared to background.**

`VIP_enrichment_iHS.R` → `VIP_enrichments.tsv`

```r
##############################################
######            READ ME           ##########
##############################################
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(ggpattern)

###########################################################################
#####       Reshape last step's gff3 df and reshape to long form      #####
###########################################################################
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet_nokin_SHAPEIT5noref_anc_lowconf/low_conf_anc_annot_iHS_inputs/new_iHS_per_gene_all_pops_lowconf_anc_annot_BCpipeline_gff3_df.RData")
gff3_results <- gff3_results %>% gather(pop, ihs, 14:47)

# combine gff3_results with metadata
metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
gff3_results <- merge(gff3_results, metadata, by.x = "pop", by.y = "pop_code", all.x = TRUE)
unique(gff3_results$region)

################################################
#####       Compare VIPs and non-VIPs      #####
################################################
VIP_list <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze1/metadata/from_providers/viral_effects_file_april2023.txt", header = FALSE)
gff3_results$type <- ifelse(gff3_results$gene_IDs %in% VIP_list$V1, "VIP", "non-VIP")


run_ihs_enrichments_VIP <- function(df_long) {
  library(data.table)
  df <- as.data.table(copy(df_long))
  df <- df[!is.na(ihs)]

  percentile_labels <- c("50", "70", "90", "95", "99", "99.5")
  percentile_cutoffs <- c(0.50, 0.70, 0.90, 0.95, 0.99, 0.995)

  results <- list()

  for (p in unique(df$pop)) {
    df_pop <- df[pop == p]
    ihs_cutoffs <- quantile(df_pop$ihs, probs = percentile_cutoffs, na.rm = TRUE)
    names(ihs_cutoffs) <- percentile_labels

    for (lbl in names(ihs_cutoffs)) {
      cut <- ihs_cutoffs[[lbl]]

      g1 <- df_pop[type == "VIP"]
      g0 <- df_pop[type == "non-VIP"]

      g1_top <- sum(g1$ihs > cut)
      g0_top <- sum(g0$ihs > cut)
      g1_total <- nrow(g1)
      g0_total <- nrow(g0)

      test <- chisq.test(
        matrix(c(g1_top, g1_total - g1_top, g0_top, g0_total - g0_top), nrow = 2),
        correct = FALSE
      )

      results[[paste(p, lbl, sep = "_")]] <- data.table(
        pop = p,
        Percentile = lbl,
        Proportion_VIP = g1_top / g1_total,
        Proportion_nonVIP = g0_top / g0_total,
        Enrichment_ratio = (g1_top / g1_total) / (g0_top / g0_total),
        n_VIP = g1_total,
        n_nonVIP = g0_total,
        count_VIP_top = g1_top,
        count_VIP_notop = g1_total - g1_top,
        count_nonVIP_top = g0_top,
        count_nonVIP_notop = g0_total - g0_top,
        chisq_pval = test$p.value,
        chisq_stat = unname(test$statistic)
      )
    }
  }

  vip_enrichment <- rbindlist(results)
  vip_enrichment$Percentile <- factor(vip_enrichment$Percentile, levels = percentile_labels, ordered = TRUE)

  return(vip_enrichment)
}


VIP_enrichment_df <- run_ihs_enrichments_VIP(gff3_results)

# Save results
fwrite(VIP_enrichment_df, "VIP_enrichments.tsv", sep = "\t")
save.image("VIP_enrichment_iHS.RData")
```

**Step 2) Test the enrichment of metabolic genes and subcategories of metabolic genes (e.g. carbohydrate, lipid, etc.) for iHS scores > 99% empirical percentile compared to background. Refer to the code page `0_annotation_curation` for documentation on how these gene sets were curated.**

`enrichment_iHS.R` → `ihs_enrichments.tsv`

```r
##########################
####      READ ME     ####
##########################

# conducting enrichment on metabolism-related gene sets on iHS scores

library(data.table)
library(dplyr)
library(ggplot2)

## ------------------------------------------------------------------
##  1. Import iHS df and metabolic subcat reference
## ------------------------------------------------------------------
setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS")
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/new_iHS_per_gene_all_pops_lowconf_anc_annot_BCpipeline_gff3_df.RData")
subcat_annot <- readRDS("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/fst/ihs/data/sub_ensg.rds")

## ------------------------------------------------------------------
##  2. Data cleaning
## ------------------------------------------------------------------

# Thin down unnecessary columns from iHS df
ihs <- gff3_results %>% select(!c(2,3,6,7))

# Build long-form subcategory table using hgnc_symbol (gene name)
subcat_dt <- rbindlist(lapply(names(subcat_annot), function(subcat_name) {
  data.table(
    hgnc_symbol = subcat_annot[[subcat_name]]$hgnc_symbol,
    subcat = subcat_name
  )
}))

# Identify genes that belong to multiple subcategories (ambiguous mappings) -- there's 414 of them
multi_subcat_genes <- subcat_dt[, .N, by = hgnc_symbol][N > 1]$hgnc_symbol

# Merge subcategory info iHS results table
# We're going to keep all the genes for now but need to create a many-to-many table
# for example, ADH1A is going to appear 6 times per population because it's associated with 6 subcats
ihs_multi <- merge(
  ihs,
  subcat_dt,
  by.x = "names",
  by.y = "hgnc_symbol",
  all.x = TRUE,
  allow.cartesian = TRUE
)

# Create metabolic vs. non-metabolic column
ihs_multi$metabolic <- ifelse(!is.na(ihs_multi$subcat), "metabolic", "non-metabolic")

# Save table of genes and annotations
write.table(ihs_multi, file = "all_pops_ihs_values_per_gene_with_annot.txt")

## ------------------------------------------------------------------
##  3. Build function for testing the 3 situations for enrichment in each population
## ------------------------------------------------------------------

run_ihs_enrichments_full_range <- function(df, pop_colname) {
  df_pop <- as.data.table(copy(df))
  df_pop[, ihs := get(pop_colname)]
  df_pop <- df_pop[!is.na(ihs)]

  # Define percentile thresholds
  percentile_labels <- c("50", "70", "90", "95", "99", "99.5")
  percentile_cutoffs <- c(0.50, 0.70, 0.90, 0.95, 0.99, 0.995)
  ihs_cutoffs <- quantile(df_pop$ihs, probs = percentile_cutoffs, na.rm = TRUE)
  pct_vec <- setNames(ihs_cutoffs, percentile_labels)

  results <- list()

  ## ---------------------------------------------------
  ## Situation 1: Metabolic vs Non-Metabolic
  ## ---------------------------------------------------
  df1 <- unique(df_pop, by = "names")
  group1 <- df1[metabolic == "metabolic"]
  group0 <- df1[metabolic == "non-metabolic"]

  for (lbl in names(pct_vec)) {
    cut <- pct_vec[[lbl]]
    g1_top <- sum(group1$ihs > cut)
    g0_top <- sum(group0$ihs > cut)
    g1_total <- nrow(group1)
    g0_total <- nrow(group0)

    results[[paste0("meta_vs_non_", lbl)]] <- data.table(
      pop = pop_colname,
      test_type = "vs_nonmetabolic",
      test_group = "metabolic",
      background = "non-metabolic",
      Percentile = lbl,
      Proportion_group = g1_top / g1_total,
      Proportion_background = g0_top / g0_total,
      Enrichment_ratio = (g1_top / g1_total) / (g0_top / g0_total),
      n_group = g1_total,
      n_background = g0_total,
      count_group_top = g1_top,
      count_group_notop = g1_total - g1_top,
      count_background_top = g0_top,
      count_background_notop = g0_total - g0_top
    )
  }

  ## ---------------------------------------------------
  ## Situations 2 & 3: Subcategory comparisons
  ## ---------------------------------------------------
  subcats <- unique(na.omit(df_pop$subcat))
  ambiguous_genes <- df_pop[!is.na(subcat), .N, by = names][N > 1]$names

  for (sc in subcats) {
    # Situation 2: Subcat vs All Other Genes
    df2 <- copy(df_pop)
    df2[, in_group := subcat == sc]
    df2[, background := subcat != sc | is.na(subcat)]
    df2 <- df2[!(names %in% ambiguous_genes & in_group == TRUE & background == TRUE)]

    for (lbl in names(pct_vec)) {
      cut <- pct_vec[[lbl]]
      g1 <- df2[in_group == TRUE]
      g0 <- df2[background == TRUE]
      g1_top <- sum(g1$ihs > cut)
      g0_top <- sum(g0$ihs > cut)
      g1_total <- nrow(g1)
      g0_total <- nrow(g0)

      results[[paste0("subcat_vs_all_", sc, "_", lbl)]] <- data.table(
        pop = pop_colname,
        test_type = "subcat_vs_all",
        test_group = sc,
        background = "all_other_genes",
        Percentile = lbl,
        Proportion_group = g1_top / g1_total,
        Proportion_background = g0_top / g0_total,
        Enrichment_ratio = (g1_top / g1_total) / (g0_top / g0_total),
        n_group = g1_total,
        n_background = g0_total,
        count_group_top = g1_top,
        count_group_notop = g1_total - g1_top,
        count_background_top = g0_top,
        count_background_notop = g0_total - g0_top
      )
    }

    # Situation 3: Subcat vs Other Metabolic Genes
    df3 <- df_pop[!is.na(subcat)]
    df3[, in_group := subcat == sc]
    df3[, background := subcat != sc]
    df3 <- df3[!(names %in% ambiguous_genes & in_group == TRUE & background == TRUE)]

    for (lbl in names(pct_vec)) {
      cut <- pct_vec[[lbl]]
      g1 <- df3[in_group == TRUE]
      g0 <- df3[background == TRUE]
      g1_top <- sum(g1$ihs > cut)
      g0_top <- sum(g0$ihs > cut)
      g1_total <- nrow(g1)
      g0_total <- nrow(g0)

      results[[paste0("subcat_vs_metab_", sc, "_", lbl)]] <- data.table(
        pop = pop_colname,
        test_type = "subcat_vs_other_metabolic",
        test_group = sc,
        background = "other_metabolic",
        Percentile = lbl,
        Proportion_group = g1_top / g1_total,
        Proportion_background = g0_top / g0_total,
        Enrichment_ratio = (g1_top / g1_total) / (g0_top / g0_total),
        n_group = g1_total,
        n_background = g0_total,
        count_group_top = g1_top,
        count_group_notop = g1_total - g1_top,
        count_background_top = g0_top,
        count_background_notop = g0_total - g0_top
      )
    }
  }

  return(rbindlist(results))
}


## ------------------------------------------------------------------
##  4. Implementation
## ------------------------------------------------------------------

# Define the columns representing each population
ihs_pop_columns <- setdiff(colnames(ihs_multi), c(
  "names", "chr", "start", "end", "gene_IDs", "length",
  "genomic_center", "window_lower", "window_upper",
  "metabolic", "subcat"
))

# Run the enrichment for each population
ihs_enrichments <- rbindlist(
  lapply(ihs_pop_columns, function(pop) {
    message("Running iHS enrichment for: ", pop)
    run_ihs_enrichments_full_range(ihs_multi, pop)
  })
)

# Save results
fwrite(ihs_enrichments, "ihs_enrichments.tsv", sep = "\t")

save.image("enrichment_iHS.RData")
```

**Step 3) Test the enrichment of endocrine genes and subcategories of endocrine genes (e.g. estrogen signaling, etc.) for iHS scores > 99% empirical percentile compared to background. Refer to the code page `0_geneset_annotations` for documentation on how these gene sets were curated.**

`endocrine_enrichment_iHS.R` → `ihs_endocrine_enrichments.tsv`

```r
## ------------------------------------------------------------------
##  README
## ------------------------------------------------------------------

# Calculating fold enrichment ratios for endocrine-related terms,
# in similar manner as the metabolic enrichment terms

library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(data.table)

## ------------------------------------------------------------------
##  1. Import iHS df and endocrine genes reference
## ------------------------------------------------------------------
setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS")
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/new_iHS_per_gene_all_pops_lowconf_anc_annot_BCpipeline_gff3_df.RData")

endocrine_genes <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/curate_genes/endocrine/all_endocrine_genes.txt", header=FALSE)
colnames(endocrine_genes) <- c("kegg_pathway_id", "kegg_pathway", "ncbi_gene_id")

endocrine_genes$ncbi_gene_id <- as.character(endocrine_genes$ncbi_gene_id)

# Use mapIds to fetch gene symbols
endocrine_genes$names <- mapIds(org.Hs.eg.db,
                                keys = endocrine_genes$ncbi_gene_id,
                                column = "SYMBOL",
                                keytype = "ENTREZID",
                                multiVals = "first")


## ------------------------------------------------------------------
##  2. Clean up endocrine gene list and merge endocrine status with iHS table
## ------------------------------------------------------------------
endocrine_genes <- as_tibble(endocrine_genes)

# Create endocrine subcategory-to-gene mapping table
subcat_dt_endocrine <- endocrine_genes %>%
  dplyr::select(hgnc_symbol = names, subcat = kegg_pathway) %>%
  dplyr::filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  data.table::as.data.table()

# Identify genes that appear in multiple endocrine subcategories
multi_subcat_genes_endocrine <- subcat_dt_endocrine[, .N, by = hgnc_symbol][N > 1]$hgnc_symbol

# Start from gff3_results
gff3_results <- gff3_results %>%
  dplyr::select(!c(2,3,6,7)) %>%
  as.data.table()

# Merge with subcategory table
ihs_multi_endocrine <- merge(
  gff3_results,
  subcat_dt_endocrine,
  by.x = "names",
  by.y = "hgnc_symbol",
  all.x = TRUE,
  allow.cartesian = TRUE
)

# Add binary annotation column
ihs_multi_endocrine$endocrine <- ifelse(!is.na(ihs_multi_endocrine$subcat), "endocrine", "non-endocrine")

# Save annotated table
fwrite(ihs_multi_endocrine, "all_pops_ihs_values_per_gene_with_endocrine_annot.txt", sep = "\t")


## ------------------------------------------------------------------
##  3. Build function for testing the 3 situations for enrichment in each population
## ------------------------------------------------------------------

run_ihs_enrichments_endocrine <- function(df, pop_colname) {
  df_pop <- as.data.table(copy(df))
  df_pop[, ihs := get(pop_colname)]
  df_pop <- df_pop[!is.na(ihs)]

  percentile_labels <- c("50", "70", "90", "95", "99", "99.5")
  percentile_cutoffs <- c(0.50, 0.70, 0.90, 0.95, 0.99, 0.995)
  ihs_cutoffs <- quantile(df_pop$ihs, probs = percentile_cutoffs, na.rm = TRUE)
  pct_vec <- setNames(ihs_cutoffs, percentile_labels)

  results <- list()

  ## ---------------------------------------------------
  ## Situation 1: Endocrine vs Non-Endocrine
  ## ---------------------------------------------------
  df1 <- unique(df_pop, by = "names")
  group1 <- df1[endocrine == "endocrine"]
  group0 <- df1[endocrine == "non-endocrine"]

  for (lbl in names(pct_vec)) {
    cut <- pct_vec[[lbl]]
    g1_top <- sum(group1$ihs > cut)
    g0_top <- sum(group0$ihs > cut)
    g1_total <- nrow(group1)
    g0_total <- nrow(group0)

    results[[paste0("endo_vs_non_", lbl)]] <- data.table(
      pop = pop_colname,
      test_type = "endocrine_vs_nonendocrine",
      test_group = "endocrine",
      background = "non-endocrine",
      Percentile = lbl,
      Proportion_group = g1_top / g1_total,
      Proportion_background = g0_top / g0_total,
      Enrichment_ratio = (g1_top / g1_total) / (g0_top / g0_total),
      n_group = g1_total,
      n_background = g0_total,
      count_group_top = g1_top,
      count_group_notop = g1_total - g1_top,
      count_background_top = g0_top,
      count_background_notop = g0_total - g0_top
    )
  }

  ## ---------------------------------------------------
  ## Situations 2 & 3: Subcategory comparisons
  ## ---------------------------------------------------
  subcats <- unique(na.omit(df_pop$subcat))
  ambiguous_genes <- df_pop[!is.na(subcat), .N, by = names][N > 1]$names

  for (sc in subcats) {
    # Situation 2: Subcat vs All Genes
    df2 <- copy(df_pop)
    df2[, in_group := subcat == sc]
    df2[, background := subcat != sc | is.na(subcat)]
    df2 <- df2[!(names %in% ambiguous_genes & in_group == TRUE & background == TRUE)]

    for (lbl in names(pct_vec)) {
      cut <- pct_vec[[lbl]]
      g1 <- df2[in_group == TRUE]
      g0 <- df2[background == TRUE]
      g1_top <- sum(g1$ihs > cut)
      g0_top <- sum(g0$ihs > cut)
      g1_total <- nrow(g1)
      g0_total <- nrow(g0)

      results[[paste0("subcat_vs_all_", sc, "_", lbl)]] <- data.table(
        pop = pop_colname,
        test_type = "subcat_vs_all",
        test_group = sc,
        background = "all_other_genes",
        Percentile = lbl,
        Proportion_group = g1_top / g1_total,
        Proportion_background = g0_top / g0_total,
        Enrichment_ratio = (g1_top / g1_total) / (g0_top / g0_total),
        n_group = g1_total,
        n_background = g0_total,
        count_group_top = g1_top,
        count_group_notop = g1_total - g1_top,
        count_background_top = g0_top,
        count_background_notop = g0_total - g0_top
      )
    }

    # Situation 3: Subcat vs Other Endocrine Genes
    df3 <- df_pop[!is.na(subcat)]
    df3[, in_group := subcat == sc]
    df3[, background := subcat != sc]
    df3 <- df3[!(names %in% ambiguous_genes & in_group == TRUE & background == TRUE)]

    for (lbl in names(pct_vec)) {
      cut <- pct_vec[[lbl]]
      g1 <- df3[in_group == TRUE]
      g0 <- df3[background == TRUE]
      g1_top <- sum(g1$ihs > cut)
      g0_top <- sum(g0$ihs > cut)
      g1_total <- nrow(g1)
      g0_total <- nrow(g0)

      results[[paste0("subcat_vs_endo_", sc, "_", lbl)]] <- data.table(
        pop = pop_colname,
        test_type = "subcat_vs_other_endocrine",
        test_group = sc,
        background = "other_endocrine",
        Percentile = lbl,
        Proportion_group = g1_top / g1_total,
        Proportion_background = g0_top / g0_total,
        Enrichment_ratio = (g1_top / g1_total) / (g0_top / g0_total),
        n_group = g1_total,
        n_background = g0_total,
        count_group_top = g1_top,
        count_group_notop = g1_total - g1_top,
        count_background_top = g0_top,
        count_background_notop = g0_total - g0_top
      )
    }
  }

  return(rbindlist(results))
}

## ------------------------------------------------------------------
##  4. Implementation
## ------------------------------------------------------------------

ihs_pop_columns <- setdiff(colnames(ihs_multi_endocrine), c(
  "names", "chr", "start", "end", "gene_IDs", "length",
  "genomic_center", "window_lower", "window_upper",
  "endocrine", "subcat"
))

ihs_endocrine_enrichments <- rbindlist(
  lapply(ihs_pop_columns, function(pop) {
    message("Running iHS enrichment for: ", pop)
    run_ihs_enrichments_endocrine(ihs_multi_endocrine, pop)
  })
)

# Save results
fwrite(ihs_endocrine_enrichments, "ihs_endocrine_enrichments.tsv", sep = "\t")
              
save.image("endocrine_enrichment_iHS.RData")    
```

**Step 4) Plot the enrichment of extreme iHS in VIPs, metabolic, and endocrine genes to create Figure 3 and Supplementary Table 3.**

`VIP_metabolic_endocrine_enrichment_comps.R` → `all_subcat_level_enrichments_95th_percentile_for_supplement.csv`

```r
## ------------------------------------------------------------------
##  README
## ------------------------------------------------------------------

# Consolidating the enrichments for VIP, metabolism, and endocrine gene sets
# for the manuscript.

# Ran chi-sq test for each population + enrichment category, then estimated
# ORs and CIs using the contigency tables

# Currently the significance levels are determined by chisq_pval < 0.01


library(dplyr)
library(ggplot2)
library(data.table)

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS")


## ------------------------------------------------------------------
##  Calculate chi-sq statistic if not available yet
## ------------------------------------------------------------------

# Calculate p-value for enrichment / depletion
enrichments <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/ihs_enrichments.tsv")
endocrine_enrichments <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/ihs_endocrine_enrichments.tsv")

# Chi-sq test for metabolism and endocrine
chisq <- enrichments %>%
  rowwise() %>%
  mutate(
    chisq_result = list(chisq.test(
      matrix(c(count_group_top, count_group_notop,
               count_background_top, count_background_notop),
             nrow = 2,
             byrow = TRUE),
      correct = FALSE
    )),
    chisq_pval = chisq_result$p.value,
    chisq_stat = chisq_result$statistic
  ) %>%
  ungroup() %>%
  dplyr::select(-chisq_result)

endocrine_chisq <- endocrine_enrichments %>%
  rowwise() %>%
  mutate(
    chisq_result = list(chisq.test(
      matrix(c(count_group_top, count_group_notop,
               count_background_top, count_background_notop),
             nrow = 2,
             byrow = TRUE),
      correct = FALSE
    )),
    chisq_pval = chisq_result$p.value,
    chisq_stat = chisq_result$statistic
  ) %>%
  ungroup() %>%
  dplyr::select(-chisq_result)

# VIP already has chi-sq calculated
vip_enrichments <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/VIP_enrichments.tsv")

# Reformat VIP table to match the other two
vip_chisq <- vip_enrichments %>%
  dplyr::mutate(
    test_type = "vip_vs_nonvip",
    test_group = "VIP",
    background = "non-VIP",
    Proportion_group = Proportion_VIP,
    Proportion_background = Proportion_nonVIP,
    n_group = n_VIP,
    n_background = n_nonVIP,
    count_group_top = count_VIP_top,
    count_group_notop = count_VIP_notop,
    count_background_top = count_nonVIP_top,
    count_background_notop = count_nonVIP_notop
  ) %>%
  dplyr::select(
    pop, test_type, test_group, background, Percentile,
    Proportion_group, Proportion_background, Enrichment_ratio,
    n_group, n_background,
    count_group_top, count_group_notop,
    count_background_top, count_background_notop,
    chisq_pval, chisq_stat
  )


# Bind 3 dfs together for plotting
all_chisq_enrichments <- rbindlist(list(
  chisq,
  endocrine_chisq,
  vip_chisq), use.names = TRUE, fill = TRUE)

# add in significance flag for plotting downstream
all_chisq_enrichments$significant <- all_chisq_enrichments$chisq_pval < 0.01


#####------------------------------------
###    PLOTTING PANEL A HIGH-LEVEL COMPS
#####------------------------------------
metadata <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")

# Merge
df_meta_full <- merge(all_chisq_enrichments, metadata, by.x = "pop", by.y = "pop_code")

# Define region order for clean plotting
pop_region_mapping <- df_meta_full %>%
  distinct(pop, region) %>%
  mutate(region = factor(region, levels = c("South Asia", "Southeast Asia", "East Asia",
                                            "Oceania", "West Central Africa", "East Africa",
                                            "Southern Africa", "East Central Africa"))) %>% arrange(region)

df_meta_full$pop <- factor(df_meta_full$pop, levels = pop_region_mapping$pop)
df_meta_full$region <- factor(df_meta_full$region, levels = levels(pop_region_mapping$region))
df_meta_full$Percentile <- factor(df_meta_full$Percentile, levels = c("50", "70", "90", "95", "99", "99.5"))


# Subset for 1 pair of populations for 6 regions with contrast (drop East Asia and Oceania cos there's no non-AGR)
# And only high-level comps
pops_to_plot <- c("AET", "PHP", "ADV", "SHL", "BCH", "CBU", "TWA", "KIG", "SAN", "NAM", "BAK", "FAN")
comps_to_plot <- c("VIP", "metabolic", "endocrine")
df_meta_subset <- df_meta_full %>% filter(pop %in% pops_to_plot) %>% filter(test_group %in% comps_to_plot) %>% filter(Percentile != 99.5)


plot_A <- ggplot(df_meta_subset,
       aes(x = Percentile,
           y = Enrichment_ratio,
           group = interaction(pop, test_group),
           color = subsistence,
           linetype = test_group)) +
  geom_line(size = 0.6) +
  geom_point(aes(fill = significant), size = 2, shape = 21, stroke = 0.3, show.legend = FALSE) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  facet_wrap(~ region, scales = "fixed", nrow = 1) +
  scale_color_manual(values = c(
    "Hunter-gatherer" = "#C11C84",
    "Agriculturalist" = "#FDAE61",
    "Pastoralist" = "#40E0D0")) +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
  scale_linetype_manual(values = c(
    "VIP" = "solid",
    "metabolic" = "longdash",
    "endocrine" = "dotted")) +
  labs(x = "iHS Percentile",
       y = "Enrichment Ratio",
       color = "Subsistence",
       linetype = "Test Group") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    strip.text = element_text(size = 8, face = "bold")
  ) +
  guides(color = guide_legend(override.aes = list(shape = NA)))

ggsave("Fig2A_broad_enrichment_curbes_by_region.pdf", plot = plot_A, width = 11, height = 3.5, units = "in", device = cairo_pdf)
#####---------------------------------------------
###    PLOTTING PANEL B SUBGROUP METABOLIC COMPS
#####---------------------------------------------

# I want to try plotting estimated odds ratios and confidence intervals
# Add odds ratio and CI columns
df_meta_full[, `:=` (
  OR = (count_group_top * count_background_notop) /
    (count_group_notop * count_background_top),

  logOR = log((count_group_top * count_background_notop) /
                (count_group_notop * count_background_top)),

  SE_logOR = sqrt(1 / count_group_top +
                    1 / count_group_notop +
                    1 / count_background_top +
                    1 / count_background_notop)
)]

df_meta_full[, `:=` (
  CI_lower = exp(logOR - 1.96 * SE_logOR),
  CI_upper = exp(logOR + 1.96 * SE_logOR)
)]


# Plot 4 subcategories for all populations, dropping those with n < 15
# Merge sample size info into main data
sample_numbers <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/samples_per_pop_baseline_filters_site_0.1_cohort_missing_nohet_SHAPEIT5phased_ancstate_nokin.txt")
df_meta_full <- merge(df_meta_full, sample_numbers, by.x = "pop", by.y = "population")

# Define test groups and correct test types
subcat_groups_metabolic <- c("Carbohydrate_metabolism", "Lipid_metabolism")
subcat_groups_endocrine <- c("estrogen_signaling_pathway")


# Filter to get the relevant rows for plotting
subplot <- df_meta_full[
  n >= 15 &
    (
      (test_group %in% subcat_groups_metabolic & test_type == "subcat_vs_other_metabolic") |
        (test_group %in% subcat_groups_endocrine & test_type == "subcat_vs_other_endocrine")
    ) &
    Percentile == "95"
]

# Define custom display labels and order
test_group_order <- c("Carbohydrate_metabolism", "Lipid_metabolism", "estrogen_signaling_pathway")
test_group_labels <- c(
  "Carbohydrate_metabolism" = "Carbohydrate Metabolism",
  "Lipid_metabolism" = "Lipid Metabolism",
  "estrogen_signaling_pathway" = "Estrogen Signaling"
)

# Apply factor levels to enforce facet order
subplot[, test_group := factor(test_group, levels = test_group_order)]

# Define subsistence order
subsistence_order <- c("Agriculturalist", "Hunter-gatherer", "Pastoralist")
subplot[, subsistence := factor(subsistence, levels = subsistence_order)]
setorder(subplot, subsistence, pop)
subplot[, pop := factor(pop, levels = unique(pop))]

# Define fill colour for circle
subplot[, fill_color := ifelse(significant, as.character(subsistence), "white")]

subplot$p_stars <- cut(subplot$chisq_pval,
                                 breaks = c(-Inf, 0.0005, 0.005, 0.05, Inf),
                                 labels = c("***", "**", "*", ""),
                                 right = FALSE
)

# Add significance flag for opacity
subplot <- subplot %>%
  mutate(SigFlag = chisq_pval < 0.05)

# Finally plot lol
plot_B <- ggplot(subplot, aes(x = OR, y = pop, color = subsistence)) +
  # Error bars with alpha by significance
  geom_errorbarh(
    aes(xmin = CI_lower, xmax = CI_upper, alpha = SigFlag),
    height = 0.5, linewidth = 0.3
  ) +
  # Points with alpha by significance
  geom_point(
    aes(alpha = SigFlag, fill = subsistence),
    size = 1.7, shape = 21, stroke = 0.7
  ) +
  # Significance asterisks
  geom_text(
    aes(label = p_stars),
    nudge_y = 0.3, size = 2.5, color = "black"
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  facet_wrap(~ test_group, scales = "free_x", nrow = 1,
             labeller = labeller(test_group = test_group_labels)) +
  scale_color_manual(values = c(
    "Hunter-gatherer" = "#C11C84",
    "Agriculturalist" = "#FDAE61",
    "Pastoralist" = "#40E0D0"
  )) +
  scale_fill_manual(values = c(
    "Hunter-gatherer" = "#C11C84",
    "Agriculturalist" = "#FDAE61",
    "Pastoralist" = "#40E0D0"
  )) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
  labs(x = "Odds Ratio", y = NULL, color = "Subsistence", fill = "Subsistence") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.05))) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 9),
    strip.background = element_rect(fill = "grey80", color = "grey90"),
    strip.text = element_text(color = "black", face = "bold",
                              margin = margin(t = 0.5, b = 0.5, l = 4, r = 4)),
    panel.spacing = unit(0.8, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.8),
    legend.position = "top"
  )

ggsave("Fig2B_subcat_ORs_by_pop.pdf", plot = plot_B, width = 9, height = 5, units = "in", device = cairo_pdf)

# Get ORs and CI ranges for paper

# Filter for carbohydrate metabolism rows
carb_ci <- subplot[test_group == "Carbohydrate_metabolism", .(pop, OR, CI_lower, CI_upper)]

range_lower <- range(carb_ci$CI_lower, na.rm = TRUE)
range_upper <- range(carb_ci$CI_upper, na.rm = TRUE)
range_or <- range(carb_ci$OR, na.rm = TRUE)

list(
  OR_range = range_or,
  CI_lower_range = range_lower,
  CI_upper_range = range_upper
)


# Filter for estrogen signalling rows
est_ci <- subplot[test_group == "estrogen_signaling_pathway", .(pop, OR, CI_lower, CI_upper)]

# Exclude some pops
range_or <- range(est_ci[!pop %in% c("NAM", "TWA", "SAN", "VAN"), OR], na.rm = TRUE)
range_ci_lower <- range(est_ci[!pop %in% c("NAM", "TWA", "SAN", "VAN"), CI_lower], na.rm = TRUE)
range_ci_upper <- range(est_ci[!pop %in% c("NAM", "TWA", "SAN", "VAN"), CI_upper], na.rm = TRUE)

list(
  OR_range = range_or,
  CI_lower_range = range_ci_lower,
  CI_upper_range = range_ci_upper
)

#####---------------------------------------------
###    OTHER STUFF
#####---------------------------------------------

# Filter df to get just significant subcats with at least 5 genes
filtered_significant_subcats <- df_meta_full[
  chisq_pval < 0.01 &
    count_group_top >= 5 &
    count_background_top >= 5 &
    Percentile == 95 &
    !(test_group %in% c("VIP", "metabolic", "endocrine"))
]

# Get table of all chi-sq enrichments at the 95th percentile for supplement
chisq_enrichments_95th_for_supp <- df_meta_full %>% filter(Percentile == 95.0) %>%
  filter(test_type == "subcat_vs_other_metabolic" | test_type == "subcat_vs_other_endocrine")

write.csv(chisq_enrichments_95th_for_supp, "all_subcat_level_enrichments_95th_percentile_for_supplement.csv",
          row.names = FALSE, quote = TRUE)

save.image("VIP_metabolic_endocrine_enrichment_comps.RData")
```

**Step 5) Plot the proportion of VIPs with a non-synonymous SNP that has iHS > 0.95% out of all other VIPs to create Supplementary Figure 6.**

`ihs_shared_enrich.Rmd {ihs nonsyn gene}`  

```r
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(plyranges)
library(rtracklayer)
library(tidyverse)

wd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/fst/ihs/"

load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/lowconf_anc_per_gene_iHS_inc_VIP_and_percentile_info_by_pop.RData")

genes <- unique(percentiles[,c("chr", "start", "end", "gene_IDs", "names", "window_lower", "window_upper")])
genes <- genes[!duplicated(genes$names),] # Some names are duplicated due to diff Ensembl ID, just take first one
gene_list <- split(genes, ~ chr)
gene_bod_list <- lapply(gene_list, function(x) makeGRangesFromDataFrame(x[,1:4], keep.extra.columns = T, ignore.strand = T, start.field = "start", end.field = "end"))
gene_win_list <- lapply(gene_list, function(x) makeGRangesFromDataFrame(x[,c(1,4:6)], keep.extra.columns = T, ignore.strand = T, start.field = "window_lower", end.field = "window_upper"))

anno_dir <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing/annotation/"

gene_bod_snps <- list()
gene_win_snps <- list()

# Get list of nonsynonymous SNPs within gene body or window regions
for (chr in names(gene_bod_list)) {
  annos <- fread(paste0(anno_dir, "chr", chr, "_anno.hg19_multianno_nogenotypes_filtered.txt"), select = c(1:5,9))
  colnames(annos) <- c("chr", "start", "end", "ref", "alt", "vartype")
  nonsyn <- annos %>% filter(vartype == "nonsynonymous SNV") %>% makeGRangesFromDataFrame(ignore.strand = T, keep.extra.columns = T)
  gene_bod_snps[[`chr`]] <- join_overlap_inner(nonsyn, gene_bod_list[[`chr`]])
  gene_win_snps[[`chr`]] <- join_overlap_inner(nonsyn, gene_win_list[[`chr`]])
}

# Filter out SNPs that are not in the nohet_nokin VCF
# First get all unique sites we extracted
nonsyn <- suppressWarnings(unlist(as(gene_win_snps, "GRangesList")))
nonsyn <- reduce_ranges(nonsyn) %>% as.data.frame() # 285,644
fwrite(nonsyn[1:2], paste0(wd, "ihs_baseline_nonsyn_targets.txt"), quote = F, col.names = F, sep = "\t")
nonsyn <- makeGRangesFromDataFrame(nonsyn)

# Extract nonsyn SNPs that exist in the VCF used
# VCF=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet_nokin_SHAPEIT5noref_anc_lowconf/freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes_anc_inclc_clean.vcf.gz
# bcftools query -T ihs_baseline_nonsyn_targets.txt -f '%CHROM\t%POS\t%REF\t%ALT\n' $VCF > ihs_nhnk_nonsyn_sites.txt
keep_nonsyn <- fread(paste0(wd, "ihs_nhnk_nonsyn_sites.txt"), select = 1:2)
colnames(keep_nonsyn) <- c("chr", "pos")
keep_nonsyn <- makeGRangesFromDataFrame(keep_nonsyn, keep.extra.columns = T, start.field = "pos", end.field = "pos")
keep_nonsyn <- split(keep_nonsyn, ~ seqnames)

# Keep only the nonsyn SNPs in the nohet nokin VCF
tmp <- lapply(names(gene_bod_snps), function(x) join_overlap_inner(gene_bod_snps[[`x`]], keep_nonsyn[[`x`]]))
names(tmp) <- names(gene_bod_snps)
gene_bod_snps <- tmp
gene_bod_snps <- suppressWarnings(unlist(as(gene_bod_snps, "GRangesList")))

tmp <- lapply(names(gene_win_snps), function(x) join_overlap_inner(gene_win_snps[[`x`]], keep_nonsyn[[`x`]]))
names(tmp) <- names(gene_win_snps)
gene_win_snps <- tmp
gene_win_snps <- suppressWarnings(unlist(as(gene_win_snps, "GRangesList")))

# Count the number of nonsyn SNPs in the gene body and windows and merge w/ ihs data
genes <- left_join(genes, plyr::count(gene_bod_snps$names), by = c("names" = "x"))
genes <- left_join(genes, plyr::count(gene_win_snps$names), by = c("names" = "x"))
colnames(genes)[8:9] <- c("n_nonsyn_gene", "n_nonsyn_window")
genes <- setnafill(genes, fill = 0, cols = 8:9)
saveRDS(genes, paste0(wd, "data/ihs_genes_num_nonsyn.rds"))

percentiles <- percentiles %>% inner_join(genes[,c(4,8,9)], by = "gene_IDs")
saveRDS(percentiles, paste0(wd, "data/ihs_percentiles_num_nonsyn.rds"))

ihs_files <- list.files(wd, pattern = "ihs.txt", full.names = T)
nonsyn_ihs <- lapply(ihs_files, function(x) makeGRangesFromDataFrame(na.omit(fread(x, col.names = c("chr", "pos", "id", "var_ihs"))), keep.extra.columns = T, ignore.strand = T, start.field = "pos", end.field = "pos"))
nonsyn_ihs <- GRangesList(nonsyn_ihs)

saveRDS(nonsyn_ihs, paste0(wd, "data/ihs_nonsyn.rds")) # Only 60,586 unique nonsyn SNPs with iHS

# Annotate vars as VIP/non-VIP
genes <- unique(percentiles[,c("chr", "start", "end", "type", "gene_IDs", "names", "window_lower", "window_upper")])
genes <- genes[!duplicated(genes$names),] # Some names are duplicated due to diff Ensembl ID, just take first one

# Combine gene info and percentiles w/ variant ihs
p95 <- percentiles %>% filter(within_pop_percentile > 95)
p95$popcode <- p95$pop
p95_list <- makeGRangesListFromDataFrame(p95, "pop", keep.extra.columns = T)

p99 <- percentiles %>% filter(within_pop_percentile > 99)
p99$popcode <- p99$pop
p99_list <- makeGRangesListFromDataFrame(p99, "pop", keep.extra.columns = T)

p95_wlist <- makeGRangesListFromDataFrame(p95[,c(-5,-6)], "pop", start.field = "window_lower", end.field = "window_upper", keep.extra.columns = T)
p99_wlist <- makeGRangesListFromDataFrame(p99[,c(-5,-6)], "pop", start.field = "window_lower", end.field = "window_upper", keep.extra.columns = T)

names(nonsyn_ihs) <- names(p95_list)

p95_olap <- lapply(1:34, function(x) join_overlap_intersect(p95_list[[x]], nonsyn_ihs[[x]]))
p95_olap <- suppressWarnings(unlist(as(p95_olap, "GRangesList"))) %>% as.data.frame()
p99_olap <- lapply(1:34, function(x) join_overlap_intersect(p99_list[[x]], nonsyn_ihs[[x]]))
p99_olap <- suppressWarnings(unlist(as(p99_olap, "GRangesList"))) %>% as.data.frame()

p95_wolap <- lapply(1:34, function(x) join_overlap_intersect(p95_wlist[[x]], nonsyn_ihs[[x]]))
p95_wolap <- suppressWarnings(unlist(as(p95_wolap, "GRangesList"))) %>% as.data.frame()
p99_wolap <- lapply(1:34, function(x) join_overlap_intersect(p99_wlist[[x]], nonsyn_ihs[[x]]))
p99_wolap <- suppressWarnings(unlist(as(p99_wolap, "GRangesList"))) %>% as.data.frame()

# Unique genes that are VIP in top 95/99 iHS
genes_95 <- p95_olap %>% dplyr::select(names, type) %>% unique() # 3661 unique genes
print(sum(genes_95 == "VIP") / nrow(genes_95)) # 0.26
print(sum(p95_olap$type == "VIP") / nrow(p95_olap)) # 0.26

genes_99 <- p99_olap %>% dplyr::select(names, type) %>% unique() # 827 unique genes
print(sum(genes_99 == "VIP") / nrow(genes_99)) # 0.30
print(sum(p99_olap$type == "VIP") / nrow(p99_olap)) # 0.27

wgenes_95 <- p95_wolap %>% dplyr::select(names, type) %>% unique() # 6112 unique genes
print(sum(wgenes_95 == "VIP") / nrow(wgenes_95)) # 0.23
print(sum(p95_wolap$type == "VIP") / nrow(p95_wolap)) # 0.22

wgenes_99 <- p99_wolap %>% dplyr::select(names, type) %>% unique() # 1345 unique genes
print(sum(wgenes_99 == "VIP") / nrow(wgenes_99)) # 0.26
print(sum(p99_wolap$type == "VIP") / nrow(p99_wolap)) # 0.20

# Proportion VIPs that exceed gene-wise iHS 95/99 percentile score (min per-population iHS)
pop_ihs_ecdfs <- lapply(names(p95_list), function(x) ecdf(percentiles[percentiles$pop == x,]$ihs))
names(pop_ihs_ecdfs) <- names(p95_list)
ihs_95 <- sapply(pop_ihs_ecdfs, function(x) quantile(x, 0.95))
names(ihs_95) <- names(p95_list)
ihs_99 <- sapply(pop_ihs_ecdfs, function(x) quantile(x, 0.99))
names(ihs_99) <- names(p95_list)

nonsyn_extreme_95 <- lapply(names(p95_list), function(x) p95_olap %>% filter(popcode == x) %>% mutate(extreme = ifelse(abs(var_ihs) >  ihs_95[[x]], T, F))) %>% bind_rows()
nonsyn_extreme_99 <- lapply(names(p95_list), function(x) p99_olap %>% filter(popcode == x) %>% mutate(extreme = ifelse(abs(var_ihs) >  ihs_99[[x]], T, F))) %>% bind_rows()
nonsyn_wextreme_95 <- lapply(names(p95_list), function(x) p95_wolap %>% filter(popcode == x) %>% mutate(extreme = ifelse(abs(var_ihs) >  ihs_95[[x]], T, F))) %>% bind_rows()
nonsyn_wextreme_99 <- lapply(names(p95_list), function(x) p99_wolap %>% filter(popcode == x) %>% mutate(extreme = ifelse(abs(var_ihs) >  ihs_99[[x]], T, F))) %>% bind_rows()

print(nrow(nonsyn_extreme_95[nonsyn_extreme_95$type == "VIP",]) / nrow(nonsyn_extreme_95)) # 0.27
print(nrow(nonsyn_extreme_99[nonsyn_extreme_99$type == "VIP",]) / nrow(nonsyn_extreme_99)) # 0.26
print(nrow(nonsyn_wextreme_95[nonsyn_wextreme_95$type == "VIP",]) / nrow(nonsyn_wextreme_95)) # 0.22
print(nrow(nonsyn_wextreme_99[nonsyn_wextreme_99$type == "VIP",]) / nrow(nonsyn_wextreme_99)) # 0.21

# Collapse extreme to gene-level (gene is extreme if at least 1 nonsyn w/ extreme score)
exgene_95 <- nonsyn_extreme_95 %>% group_by(popcode, names) %>% summarise(gene_extreme = ifelse(sum(extreme > 0), T, F), type = type[1])
exgene_95 <- mutate(exgene_95, label = paste(ifelse(gene_extreme == 1, "Extreme", "Normal"), type))
exgene_99 <- nonsyn_extreme_99 %>% group_by(popcode, names) %>% summarise(gene_extreme = ifelse(sum(extreme > 0), T, F), type = type[1])
exgene_99 <- mutate(exgene_99, label = paste(ifelse(gene_extreme == 1, "Extreme", "Normal"), type))
exwind_95 <- nonsyn_wextreme_95 %>% group_by(popcode, names) %>% summarise(gene_extreme = ifelse(sum(extreme > 0), T, F), type = type[1])
exwind_95 <- mutate(exwind_95, label = paste(ifelse(gene_extreme == 1, "Extreme", "Normal"), type))
exwind_99 <- nonsyn_wextreme_99 %>% group_by(popcode, names) %>% summarise(gene_extreme = ifelse(sum(extreme > 0), T, F), type = type[1])
exwind_99 <- mutate(exwind_99, label = paste(ifelse(gene_extreme == 1, "Extreme", "Normal"), type))

pdf(paste0(wd, "figures/ihs_extreme_nonsyn_genes.pdf"), height = 8, width = 12)
exgene_95 %>% group_by(popcode, label) %>% summarise(n = n(), .groups = "keep") %>%
ggplot(aes(x = popcode, y = n, fill = label)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), size = 3, position = position_fill(0.5)) +
  scale_fill_manual(values = c("lightblue", "maroon", "steelblue", "deeppink4")) +
  labs(x = "Population", y = "Proportion of genes", fill = "Gene type", title = "iHS 95 percentile genes that contain nonsynonymous variant (gene body)") +
  theme_classic()

exgene_99 %>% group_by(popcode, label) %>% summarise(n = n(), .groups = "keep") %>%
ggplot(aes(x = popcode, y = n, fill = label)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), size = 3, position = position_fill(0.5)) +
  scale_fill_manual(values = c("lightblue", "maroon", "steelblue", "deeppink4")) +
  labs(x = "Population", y = "Proportion of genes", fill = "Gene type", title = "iHS 99 percentile genes that contain nonsynonymous variant (gene body)") +
  theme_classic()

exwind_95 %>% group_by(popcode, label) %>% summarise(n = n(), .groups = "keep") %>%
ggplot(aes(x = popcode, y = n, fill = label)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), size = 3, position = position_fill(0.5)) +
  scale_fill_manual(values = c("darkorchid", "goldenrod", "purple4", "orange4")) +
  labs(x = "Population", y = "Proportion of genes", fill = "Gene type", title = "iHS 95 percentile genes that contain nonsynonymous variant (100kb_window)") +
  theme_classic()

exwind_99 %>% group_by(popcode, label) %>% summarise(n = n(), .groups = "keep") %>%
ggplot(aes(x = popcode, y = n, fill = label)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), size = 3, position = position_fill(0.5)) +
  scale_fill_manual(values = c("darkorchid", "goldenrod", "purple4", "orange4")) +
  labs(x = "Population", y = "Proportion of genes", fill = "Gene type", title = "iHS 99 percentile genes that contain nonsynonymous variant (100kb window)") +
  theme_classic()
dev.off()
```

## Enrichment of regulatory features

**Step 1) Average absolute iHS scores in 100kb non-overlapping windows across the genome.**

`new_iHS_per_window_100kb.R` → `new_iHS_per_window_100kb_df.RData`

```r
#########################################################################################
#####                                   READ ME                                     #####
#########################################################################################

# computation of standardised iHS values per 100kb window for all populations.

# this script should parallelize looping through populations using foreach

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet_nokin_SHAPEIT5noref_anc_lowconf/low_conf_anc_annot_iHS_inputs")
library(data.table)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)

output_dir <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet_nokin_SHAPEIT5noref_anc_lowconf/low_conf_anc_annot_iHS_inputs"
files <- list.files(path = output_dir, pattern = "\\_iHS$",full.names = TRUE)

# set up temporary data frame to load all final contents
all_pops_df <- data.frame()

# read each file within the range
for (i in 1:length(files)) {
  currentFile <- fread(files[i], col.names = c("index", "rsID", "freq", "iHH_0", "iHH_1", "iHS", "std_iHS"),
                       stringsAsFactors = TRUE)
  currentFileName <- gsub(".*/([A-Za-z]+)_.*", "\\1", files[i])
  currentFile$pop <- c(currentFileName)

  # append the current file to all_pops_df
  all_pops_df <- rbind(all_pops_df, currentFile)
}

all_pops_df <- na.omit(all_pops_df)


######################################################################################
#####        Format all_pops_df with the fake coordinates from rsIDs             #####
######################################################################################
# Extract the chromosome number
# This regex captures digits following 'chr' and before 'pos'
all_pops_df$chr <- as.numeric(str_extract(all_pops_df$rsID, "(?<=chr)\\d+"))
all_pops_df$chr <- as.numeric(all_pops_df$chr)

# Extract the position number
# This regex captures digits following 'pos'
all_pops_df$pos <- as.numeric(str_extract(all_pops_df$rsID, "(?<=pos)\\d+"))
all_pops_df$pos <- as.numeric(all_pops_df$pos)

all_pops_df <- all_pops_df %>% arrange(pop, index)
all_pops_df <- all_pops_df[, c(1,9,10,3:8)]

# check structure
head(all_pops_df)
str(all_pops_df)

#######################################################################
#####    Get dataframe of chromosome lengths for build 37/hg19    #####
#######################################################################

chr_lengths <- read.csv("/project/lbarreiro/USERS/bridget/references/chromosome_lengths_GRCh37.csv", header=FALSE)
chr_lengths <- chr_lengths[1:22,1:2]
chr_lengths$V2 <- as.numeric(gsub(",", "", chr_lengths$V2))
chr_lengths <- setNames(chr_lengths$V2, paste0("chr", chr_lengths$V1)) # convert to named vector

# Function to create intervals
create_intervals <- function(chr, length, interval_size = 100000) {
  starts <- seq(1, length, by = interval_size)
  ends <- pmin(starts + interval_size - 1, length)
  data.frame(
    chromosome = chr,
    start = starts,
    end = ends
  )
}

# Generate intervals for all chromosomes
genome_intervals <- do.call(rbind, lapply(names(chr_lengths), function(chr) {
  create_intervals(chr, chr_lengths[chr])
}))

genome_intervals$chromosome <- gsub("chr", "", genome_intervals$chromosome)

write.csv(genome_intervals, "GRCh37_100kb_intervals.csv", row.names = FALSE)
save(genome_intervals, file= "GRCh37_100kb_intervals.RData")

######################################################################################
#####             Compute iHS per interval across populations                    #####
######################################################################################

# function for calculating average |iHS| within the interval using all variants
compute_ihs <- function(x, pop_df){ #notice pop_df is in there as an argument. this is important
  chrom <- as.numeric(as.character(x[1]))
  lower_bound <- as.numeric(as.character(x[2])) # lower bound position of the interval
  upper_bound <- as.numeric(as.character(x[3])) # upper bound position of the interval
  snps_for_each_int <- pop_df[pop_df$chr == chrom & pop_df$pos >= lower_bound & pop_df$pos <= upper_bound, ]
  ihs_for_interval <- mean(abs(snps_for_each_int$std_iHS))
  return(ihs_for_interval)
}

# function to count the number of SNPs within each interval
count_snps <- function(x, pop_df){
  chrom <- as.numeric(as.character(x[1]))
  lower_bound <- as.numeric(as.character(x[2])) # lower bound position of the interval
  upper_bound <- as.numeric(as.character(x[3])) # upper bound position of the interval
  snps_for_each_int <- pop_df[pop_df$chr == chrom & pop_df$pos >= lower_bound & pop_df$pos <= upper_bound, ]
  snp_count <- nrow(snps_for_each_int) # Count the number of SNPs
  return(snp_count) # Return the count
}

# set number of cores -----------------------------------------------------
workers <- 6
cl <- parallel::makeCluster(workers)
# register the cluster for using foreach
doParallel::registerDoParallel(cl)

# loop through populations and apply function within for-loop
unique_pops <- unique(all_pops_df$pop)

# calculate the average iHS score for each interval for each population
ihs_100kb_interval <- foreach (i=1:length(unique_pops), .combine = cbind, .packages=c('data.table', 'dplyr', 'foreach', 'doParallel')) %dopar% {
  pop_i <- unique_pops[i]
  pop_df <- all_pops_df %>% filter(pop == pop_i)
  print(unique(pop_df$pop))
  as.vector(apply(genome_intervals, 1, compute_ihs, pop_df))
}

# calculate the number of SNPs in each interval for each population
nsnps_100kb_interval <- foreach (i=1:length(unique_pops), .combine = cbind, .packages=c('data.table', 'dplyr', 'foreach', 'doParallel')) %dopar% {
  pop_i <- unique_pops[i]
  pop_df <- all_pops_df %>% filter(pop == pop_i)
  print(unique(pop_df$pop))
  as.vector(apply(genome_intervals, 1, count_snps, pop_df))
}

save.image(file="/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet_nokin_SHAPEIT5noref_anc_lowconf/low_conf_anc_annot_iHS_inputs/new_iHS_per_window_100kb.RData")
save(ihs_100kb_interval, file="new_iHS_per_window_100kb_df.RData")
save(nsnps_100kb_interval, file="new_iHS_per_window_100kb_nsnps.df.RData")
```

**Step 2) Use logistic regression to test the enrichment of immune, brain, and metabolic tissue (adipose-tissue, liver, and pancreas) cCREs separately for windowed iHS scores > 99% empirical percentile compared to background. Aggregate ORs across populations with at least 15 individuals using Fisher's combined method. cCRE annotations are included in `${tissue}_cCREs_GRCh37.tsv`.**

`enrich_cCREs_100kb.R ${tissue}` → `${tissue}_logistic_regression_ORs_CIs_pvals_all_pops_100kb.txt`, `${tissue}_iHS_100kb_with_cCRE_annot_all_pops.txt`, and `${tissue}_cCREs_100kb_fishers_meta.rds`

```r
#####################################
#####         READ ME           #####
#####################################

# annotating 100kb iHS windows with putative cCRE for a given tissue regulatory elements
# run this script for all populations.

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/ihs/ccres/")
library(data.table)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(scales)
args <- commandArgs(T)
tissue <- args[1]

# import the putative cCREs & 100kb stats
cCREs <- read.delim(paste0("data/", tissue, "_cCREs_GRCh37.tsv"))
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/functional_annotations/new_iHS_per_window_100kb_df.RData")
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/functional_annotations/new_iHS_per_window_100kb_nsnps.df.RData")
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/functional_annotations/GRCh37_100kb_intervals.RData")
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/unique_pops.RData")

# recode iHS window value as NaN if window has less than 5 snps
ihs_100kb_interval[nsnps_100kb_interval < 5] <- NaN
cCREs$chr <- gsub("^chr", "", cCREs$chr)

colnames(ihs_100kb_interval) <- unique_pops

########################################################################
#####      first annotate each 100kb window with the cCREs        ######
########################################################################
genome_intervals_gr <- GRanges(
  seqnames = genome_intervals$chromosome,
  ranges = IRanges(start = genome_intervals$start, end = genome_intervals$end)
)

cCREs_gr <- GRanges(
  seqnames = cCREs$chr,
  ranges = IRanges(start = cCREs$start, end = cCREs$end),
  cCRE = cCREs$cCRE
)

# Get unique classes of cCREs
unique_classes <- unique(cCREs$cCRE)

# Initialize new columns in genome_intervals for each class
for (c in unique_classes) {
  genome_intervals[[`c`]] <- 0
}

# Check for overlaps and update the columns
overlaps <- findOverlaps(genome_intervals_gr, cCREs_gr)

# For each overlap, set the corresponding class column to 1
for (i in seq_along(overlaps)) {
  genome_row <- queryHits(overlaps)[i]
  row <- subjectHits(overlaps)[i]
  cCRE <- cCREs$cCRE[row]
  genome_intervals[genome_row, cCRE] <- 1
}

###############################################################
####      Run Logistic Regression for all populations      ####
###############################################################

ihs_100kb_interval <- cbind(genome_intervals, ihs_100kb_interval)
colnames(ihs_100kb_interval) <- gsub("-", "_", colnames(ihs_100kb_interval))

# Initialize dataframes to store results for each regulatory element
results_pELS <- data.frame(Population = character(), Odds_Ratio = numeric(), Lower_CI = numeric(), Upper_CI = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
results_dELS <- data.frame(Population = character(), Odds_Ratio = numeric(), Lower_CI = numeric(), Upper_CI = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
results_PLS <- data.frame(Population = character(), Odds_Ratio = numeric(), Lower_CI = numeric(), Upper_CI = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
results_H3K4me3 <- data.frame(Population = character(), Odds_Ratio = numeric(), Lower_CI = numeric(), Upper_CI = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

# Function to run logistic regression and store results, including confidence intervals
run_logistic_regression <- function(element_col, results_df) {
  for (pop in unique_pops) {
    # Compute the threshold for high iHS (99th percentile) for each population
    threshold <- quantile(ihs_100kb_interval[[`pop`]], 0.99, na.rm = TRUE)

    # Create a binary high iHS column specifically for the current population
    high_iHS_col <- paste0(pop, "_high_iHS")
    ihs_100kb_interval[[`high_iHS_col`]] <- ifelse(ihs_100kb_interval[[`pop`]] > threshold, 1, 0)

    # Fit the logistic regression model using the binary high iHS column
    model <- glm(formula = as.formula(paste(high_iHS_col, "~", element_col)), data = ihs_100kb_interval, family = binomial)

    # Extract coefficients and standard errors
    log_odds <- coef(model)[element_col]
    standard_error <- summary(model)$coefficients[element_col, "Std. Error"]
    odds_ratio <- exp(log_odds)

    # Calculate 95% CI for the odds ratio
    conf_int <- exp(confint(model)[element_col, ])
    lower_ci <- conf_int[1]
    upper_ci <- conf_int[2]

    # Extract p-value
    p_value <- summary(model)$coefficients[element_col, "Pr(>|z|)"]

    # Append results to the provided results dataframe
    results_df <- rbind(results_df, data.frame(
      Population = pop,
      Log_Odds_Ratio = log_odds,
      Standard_Error = standard_error,
      Odds_Ratio = odds_ratio,
      Lower_CI = lower_ci,
      Upper_CI = upper_ci,
      P_Value = p_value,
      row.names = NULL
    ))
  }
  return(results_df)
}

# Run logistic regression for each regulatory element
results_pELS <- run_logistic_regression("pELS", results_pELS)
results_dELS <- run_logistic_regression("dELS", results_dELS)
results_PLS <- run_logistic_regression("PLS", results_PLS)
results_H3K4me3 <- run_logistic_regression("CA_H3K4me3", results_H3K4me3)

# Add element type column to each dataframe
results_pELS$element_type <- "pELS"
results_dELS$element_type <- "dELS"
results_PLS$element_type <- "PLS"
results_H3K4me3$element_type <- "CA_H3K4me3"

# Combine results into one dataframe
combined_results <- rbind(results_pELS, results_dELS, results_PLS, results_H3K4me3)

# Write the combined results to a file
write.table(combined_results, file = paste0("data/", tissue, "_logistic_regression_ORs_CIs_pvals_all_pops_100kb.txt"), row.names = FALSE)

# Save ihs_100kb_interval df that has annotations for presence of cCREs
write.table(ihs_100kb_interval, file = paste0("data/", tissue, "_iHS_100kb_with_cCRE_annot_all_pops.txt"), row.names = FALSE)


################################
####      Meta p-value      ####
################################

# Only test populations with at least 15 individuals
metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
combined_results <- merge(combined_results, metadata, by.x = "Population", by.y = "pop_code")
sample_numbers <- read.csv("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/samples_per_pop_baseline_filters_site_0.1_cohort_missing_nohet_SHAPEIT5phased_ancstate_nokin.txt", sep="")

combined_results <- merge(combined_results, sample_numbers, by.x = "Population", by.y = "population")

filtered_combined_results <- combined_results %>% filter(n>=15)

# Fisher's combined method uses the sum of the natural log of the p-vals from different tests, then you compare this combined
# statistic against a chi-squared distribution with 2k degrees of freedom
fishers_combined_pval <- filtered_combined_results %>%
  group_by(element_type) %>%
  summarise(Combined_P_Value = pchisq(-2 * sum(log(P_Value)), df = 2 * n(), lower.tail = FALSE))

saveRDS(fishers_combined_pval, paste0("data/", tissue, "_cCREs_100kb_fishers_meta.rds"))
```

**Step 3) Use logistic regression including immune, brain, and metabolic tissue cCREs as covariates in the same model to test enrichment for windowed iHS scores > 99% empirical percentile compared to background. Aggregate ORs across populations with at least 15 individuals using Fisher's combined method.**

`enrich_cCREs_100kb_multiregression.R` → `multiregression_tissue-specific.txt` and `fe_100kb_all_tissues.txt`

```r
#####################################
#####         READ ME           #####
#####################################

# annotating 100kb iHS windows with putative cCRE regulatory elements
# & running logistic regression model with combinations or all cCREs as covariates

library(data.table)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(scales)
library(metafor)

load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/unique_pops.RData")

# Import 100kb iHS cCRE annots for each tissue
# iHS annots are the same per population anyways so drop those cols after first
# Select only cols for PLS, pELS, dELS
immune_ccres <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/functional_annotations/iHS_100kb_with_i-cCRE_annot_all_pops.txt",
                     select = c(1,2,3,8,4,7,13:46))
brain_ccres <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/functional_annotations/iHS_100kb_with_b-cCRE_annot_all_pops.txt",                    select = c(1,2,6,7,5))
adipose_ccres <- fread("data/adipose-tissue_iHS_100kb_with_cCRE_annot_all_pops.txt",
                      select = c(1,2,4,5,6))
pancreas_ccres <- fread("data/pancreas_iHS_100kb_with_cCRE_annot_all_pops.txt",
                       select = c(1,2,7,4,9))
liver_ccres <- fread("data/liver_iHS_100kb_with_cCRE_annot_all_pops.txt",
                    select = c(1,2,6,7,5))

# Rename elements so they're different cols
colnames(immune_ccres)[4:6] <- paste0("immune_", colnames(immune_ccres)[4:7])
colnames(brain_ccres)[3:5] <- paste0("brain_", colnames(brain_ccres)[3:5])
colnames(adipose_ccres)[3:5] <- paste0("adipose_", colnames(adipose_ccres)[3:5])
colnames(pancreas_ccres)[3:5] <- paste0("pancreas_", colnames(pancreas_ccres)[3:5])
colnames(liver_ccres)[3:5] <- paste0("liver_", colnames(liver_ccres)[3:5])

# Merge all annots based on windows
ihs_100kb <- left_join(immune_ccres, brain_ccres, by = c("chromosome", "start")) %>%
       left_join(adipose_ccres, by = c("chromosome", "start")) %>%
       left_join(pancreas_ccres, by = c("chromosome", "start")) %>%
       left_join(liver_ccres, by = c("chromosome", "start"))
fwrite(ihs_100kb, "data/ihs_100kb_ccres_all_tissues.txt", sep = "\t")
ihs_100kb <- fread("data/ihs_100kb_ccres_all_tissues.txt")

###########################################################################
#####       Analyse correlation between reg element & iHS scores      #####
###########################################################################

# Function to run logistic regression with designated covariates
run_logistic_regression_with_covariates <- function(element_col, covariates) {
  results_df <- data.frame(
    Population = character(),
    Predictor = character(),      # To store the name of the predictor (main or covariate)
    Log_Odds_Ratio = numeric(),   # Log-transformed coefficient
    Odds_Ratio = numeric(),       # Exponentiated odds ratio
    Lower_CI = numeric(),         # Lower bound of confidence interval
    Upper_CI = numeric(),         # Upper bound of confidence interval
    Std_Error = numeric(),        # Standard error of the coefficient
    P_Value = numeric(),          # P-value of the predictor
    stringsAsFactors = FALSE
  )

  for (pop in unique_pops) {
    # Compute the threshold for high iHS (99th percentile) for each population
    threshold <- quantile(ihs_100kb[[pop]], 0.99, na.rm = TRUE)

    # Create a binary high iHS column specifically for the current population
    high_iHS_col <- paste0(pop, "_high_iHS")
    ihs_100kb[[high_iHS_col]] <- ifelse(ihs_100kb[[pop]] > threshold, 1, 0)

    # Construct the formula dynamically, including covariates
    formula <- as.formula(
      paste(high_iHS_col, "~", paste(c(element_col, covariates), collapse = " + "))
    )

    # Fit the logistic regression model
    model <- glm(formula = formula, data = ihs_100kb, family = binomial)

    # Extract coefficients, standard errors, and p-values for all predictors
    coef_summary <- summary(model)$coefficients
    conf_intervals <- exp(confint(model))

    # Loop over all predictors (main predictor and covariates)
    for (predictor in c(element_col, covariates)) {
      log_odds_ratio <- coef_summary[predictor, "Estimate"]
      std_error <- coef_summary[predictor, "Std. Error"]
      p_value <- coef_summary[predictor, "Pr(>|z|)"]

      # Calculate odds ratio and 95% confidence intervals
      odds_ratio <- exp(log_odds_ratio)
      lower_ci <- conf_intervals[predictor, 1]
      upper_ci <- conf_intervals[predictor, 2]

      # Append results to the dataframe
      results_df <- rbind(results_df, data.frame(
        Population = pop,
        Predictor = predictor,
        Log_Odds_Ratio = log_odds_ratio,
        Odds_Ratio = odds_ratio,
        Lower_CI = lower_ci,
        Upper_CI = upper_ci,
        Std_Error = std_error,
        P_Value = p_value,
        row.names = NULL
      ))
    }
  }
  # Add the predictors used to a column
  results_df$formula <- paste(c(element_col, covariates), collapse = " + ")
  return(results_df)
}

# 4: Put everything in the same model per-tissue and for all tissues -> meta analysis
tissues <- c("immune", "brain", "adipose", "pancreas", "liver", "immune|brain|adipose|pancreas|liver")

res4 <- lapply(tissues, function(x) {
              tmp <- grep(x, colnames(ihs_100kb), value = T)
              return(run_logistic_regression_with_covariates(tmp[1], tmp[-1]))
               })
combined_results4 <- bind_rows(res4)
combined_results4 <- merge(combined_results4, metadata, by.x = Population", by.y = "pop_code")
combined_results4 <- merge(combined_results4, sample_numbers, by.x = Population", by.y = "population")
fwrite(combined_results4, "data/multiregression_tissue-specific.txt", sep = "\t")
combined_results4 <- fread("data/multiregression_tissue-specific.txt")

###################################################################
#####      Fixed Effect Meta-Analysis of the Single Model     #####
###################################################################

run_fixed_effects <- function(predictor, data, string) {
        meta <- data %>%
                dplyr::filter(Predictor == predictor) %>%
                dplyr::select(Population, Log_Odds_Ratio, Std_Error) %>%
                dplyr::mutate(Variance = Std_Error^2)
        meta_res <- rma(data = meta, method = "FE",
                        yi = Log_Odds_Ratio, vi = Variance)
        meta_res$formula <- string
        return(meta_res)
}

# Calculate fixed effects for single-tissue models and the all tissues model
all_ccres <- grep("immune|brain|adipose|pancreas|liver", colnames(ihs_100kb), value = T)
all_formula <- paste(all_ccres, collapse = " + ")

fe_single_tissue <- lapply(all_ccres, function(x)
                           run_fixed_effects(predictor = x, string = "single",
                                             data = dplyr::filter(combined_results4, formula != all_formula)))
fe_all_tissues <- lapply(all_ccres, function(x)
                         run_fixed_effects(predictor = x, string = "all",
                                            data = dplyr::filter(combined_results4, formula == all_formula)))
save(fe_single_tissue, fe_all_tissues, file = "data/fe_test.RData")

fixed_effect_data <- rbind(bind_rows(fe_single_tissue), bind_rows(fe_all_tissues))

fixed_effect_data <- fixed_effect_data %>%
  dplyr::mutate(

fixed_effect_data <- fixed_effect_data %>%
  dplyr::mutate(
    Odds_Ratio = exp(Estimate),
    CI_Lower_OR = exp(CI_Lower),
    CI_Upper_OR = exp(CI_Upper))


# Write the fixed_effect_data results to a file
#write.table(fixed_effect_data, file = "data/fe_100kb_all_tissues.txt", row.names = FALSE)
```

**Step 4) Run meta-analysis of single-cCRE models within subsistence groups rather than all populations to create Supplementary Table 4.**

`hg-ag_fe_models.R` → `supp_tab_4_hg_ag_meta_pairwise_pval.csv`

```r
library(ggplot2)
library(dplyr)
library(data.table)
library(metafor)

# Supplementary Table 4: Generate subsistence-based meta OR and p-vals for metabolic tissues
tissues <- c("adipose-tissue", "liver", "pancreas")
tissue_logres <- lapply(tissues, function(x) fread(paste0("./data/", x, "_logistic_regression_ORs_CIs_pvals_all_pops_100kb.txt")))
tissue_logres <- lapply(1:3, function(x) mutate(tissue_logres[[x]], tissue = tissues[x]))
combined_results <- bind_rows(tissue_logres)

metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
sample_numbers <- read.csv("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/samples_per_pop_baseline_filters_site_0.1_cohort_missing_nohet_SHAPEIT5phased_ancstate_nokin.txt", sep="")

combined_results <- merge(combined_results, metadata, by.x = "Population", by.y = "pop_code")
combined_results <- merge(combined_results, sample_numbers, by.x = "Population", by.y = "population")
combined_results <- combined_results %>% filter(n>=15)

# Run FE meta analysis within each tissue and subsistence group
run_meta_by_group_single <- function(df, tis, predictor, group) {
  sub_df <- filter(df, tissue == tis & element_type == predictor & subsistence == group)

  if (nrow(sub_df) < 2) {
    warning(paste("Not enough data for", predictor, "in", group))
    return(NULL)
  }

  meta_result <- rma(
    yi = Log_Odds_Ratio,
    sei = Standard_Error,
    method = "FE",
    data = sub_df
  )

  data.frame(
    Tissue = tis,
    Predictor = predictor,
    Group = group,
    Log_OR = meta_result$b,
    SE = meta_result$se,
    CI_Lower = meta_result$ci.lb,
    CI_Upper = meta_result$ci.ub,
    OR = exp(meta_result$b),
    OR_CI_Lower = exp(meta_result$ci.lb),
    OR_CI_Upper = exp(meta_result$ci.ub),
    P_Value = meta_result$pval
  )
}

predictors_of_interest <- c("pELS", "dELS", "PLS")
subsistence_groups <- unique(metadata$subsistence)
test_df <- expand.grid(tissues, predictors_of_interest, subsistence_groups)

meta_results <- apply(test_df, 1, function(x) run_meta_by_group_single(combined_results, x[1], x[2], x[3]))
meta_results <- bind_rows(meta_results)

# Compute p-values via Z-test for pairwise comparisons
calculate_group_diff_pval <- function(summary_df, tissue, predictor, group1, group2) {
  est1 <- summary_df %>% filter(Tissue == tissue, Predictor == predictor, Group == group1) %>% pull(Log_OR)
  se1 <- summary_df %>% filter(Tissue == tissue, Predictor == predictor, Group == group1) %>% pull(SE)
  est2 <- summary_df %>% filter(Tissue == tissue, Predictor == predictor, Group == group2) %>% pull(Log_OR)
  se2 <- summary_df %>% filter(Tissue == tissue, Predictor == predictor, Group == group2) %>% pull(SE)

  if (length(est1) == 0 | length(est2) == 0) {
    return(NULL)
  }

  z <- (est1 - est2) / sqrt(se1^2 + se2^2)
  p_val <- 2 * pnorm(-abs(z))
  group_name <- paste0(group1, "_vs_", group2)

  data.frame(
    Tissue = tissue,
    Predictor = predictor,
    Group = group_name,
    Z = z,
    Diff_P_Value = p_val
  )
}

test_df2 <- expand.grid(tissues, predictors_of_interest, subsistence_groups, subsistence_groups) %>%
  filter(Var3 != Var4) %>%
  filter(Var3 == "Pastoralist" | (Var3 == "Agriculturalist" & Var4 == "Hunter-gatherer"))

pairwise_comparisons <- apply(test_df2, 1, function(x) calculate_group_diff_pval(meta_results, x[1], x[2], x[3], x[4]))
pairwise_comparisons <- bind_rows(pairwise_comparisons)

# Combine tables w/ immune results
i_dir <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/functional_annotations"
i_meta <- fread(paste0(i_dir, "/meta_results_by_subsistence_single_predictor.txt"))
i_pair <- fread(paste0(i_dir, "/meta_by_subsistence_single_predictor_pairwise_pvals.txt"))

i_meta$Tissue <- "immune"
meta_results <- rbind(meta_results, i_meta) %>% arrange(Tissue)

i_pair$Tissue <- "immune"
i_pair$Group <- case_when(
  i_pair$Group2 == "Hunter-gatherer" ~ "Agriculturalist_vs_Hunter-gatherer",
  i_pair$Group1 == "Agriculturalist" ~ "Pastoralist_vs_Agriculturalist",
  T ~ "Pastoralist_vs_Hunter-gatherer"
)
pairwise_comparisons <- rbind(pairwise_comparisons, select(i_pair, -Group1, -Group2)) %>% arrange(Group, Tissue, Predictor)
pairwise_comparisons$Significant <- pairwise_comparisons$Diff_P_Value < 0.05

fwrite(meta_results, "data/meta_results_by_subsistence_single_predictor.txt")
fwrite(pairwise_comparisons, "data/meta_by_subsistence_single_predictor_pairwise_pvals.txt")

# For supplementary table only keep hg vs ag comparisons
supp_tab_4 <- filter(pairwise_comparisons, Group == "Agriculturalist_vs_Hunter-gatherer")
fwrite(supp_tab_4, "figures/supp_tab_4_hg_ag_meta_pairwise_pval.csv")
```

**Step 5) Perform logistic regression enrichment and fixed-effects meta analysis across populations including all cell type-specific immune cCRE covariates in a single model. Cell type assignments are included as part of the immune cCRE data `i-cCREs_GRCh37.tsv`.**

`immune_cCREs_100kb_CTspecific_condensed_model.R` → `immune_cCREs_100kb_CTspecific_condensed_model.RData`

```r
#####################################
#####         READ ME           #####
#####################################

# Fixed effect meta-analysis on cell-type-specific immune cCRE 100kb logistic regression results
# This version models all cell types together without considering the *present vs. *specific scenarios as v1 did

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/functional_annotations")
library(data.table)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(metafor)
library(stringr)

############################################
#####      Set up data structures      #####
############################################

# import the putative i-cCREs & 100kb stats
immune_cCREs <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/functional_annotations/i-cCREs_GRCh37.tsv")
load("new_iHS_per_window_100kb_df.RData")
load("new_iHS_per_window_100kb_nsnps.df.RData")
load("GRCh37_100kb_intervals.RData")
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/unique_pops.RData")

# recode iHS window value as NaN if window has less than 5 snps
ihs_100kb_interval[nsnps_100kb_interval < 5] <- NaN

immune_cCREs$chr <- gsub("^chr", "", immune_cCREs$chr)
immune_cCREs$class <- gsub("-", "_", immune_cCREs$class)
immune_cCREs <- immune_cCREs %>% filter(class == "pELS" | class == "dELS" | class == "PLS" | class == "CA_H3K4me3")
immune_cCREs <- immune_cCREs %>% filter(chr != "X" & chr != "Y")
colnames(ihs_100kb_interval) <- unique_pops
ihs_100kb_interval <- cbind(genome_intervals, ihs_100kb_interval)

# Import metadata for plotting
metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
sample_numbers <- read.csv("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/samples_per_pop_baseline_filters_site_0.1_cohort_missing_nohet_SHAPEIT5phased_ancstate_nokin.txt", sep="")

# Columns to check in immune_cCREs
columns_to_check <- c("Calderon.Corces_Bcells", "Calderon.Corces_CD4_Tcells", "Calderon.Corces_CD8_Tcells",
                      "Calderon.Corces_Erythroblasts", "Calderon.Corces_Myeloid", "Calderon.Corces_NK",
                      "Calderon.Corces_Progenitors", "ENCODE_Bcells", "ENCODE_CD4_Tcells", "ENCODE_CD8_Tcells",
                      "ENCODE_Myeloid", "ENCODE_NK", "ENCODE_Progenitors")



######################################################################################
#####    Annotate each 100kb window with 28 columns: (CellType)_(cCRE_class)      ####
######################################################################################

cell_types <- c("Bcells", "CD4_Tcells", "CD8_Tcells","Erythroblasts", "Myeloid", "NK", "Progenitors")

ccre_classes <- c("pELS", "dELS", "PLS")

# Convert genome_intervals to GRanges once
genome_intervals_gr <- GRanges(
  seqnames = genome_intervals$chromosome,
  ranges   = IRanges(start = genome_intervals$start, end = genome_intervals$end)
)

for (ct in cell_types) {

  # We know each cell type should have a Calderon column
  col_calderon <- paste0("Calderon.Corces_", ct)
  has_calderon <- col_calderon %in% colnames(immune_cCREs)  # Should be TRUE if the data has that column

  # We also *might* have an ENCODE column
  col_encode <- paste0("ENCODE_", ct)
  has_encode <- col_encode %in% colnames(immune_cCREs)      # FALSE for Erythroblasts

  for (class_i in ccre_classes) {

    new_col <- paste0(ct, "_", class_i)
    genome_intervals[[new_col]] <- 0  # initialize

    if (has_calderon && has_encode) {
      # If we have both columns, filter on either column being TRUE
      cCRE_subset <- immune_cCREs %>%
        filter(
          class == class_i &
            (
              (.data[[col_calderon]] == TRUE | .data[[col_calderon]] == "TRUE") |
                (.data[[col_encode]]   == TRUE | .data[[col_encode]]   == "TRUE")
            )
        )
    } else if (has_calderon && !has_encode) {
      # If we only have Calderon, filter just on that
      cCRE_subset <- immune_cCREs %>%
        filter(
          class == class_i &
            (.data[[col_calderon]] == TRUE | .data[[col_calderon]] == "TRUE")
        )
    } else if (!has_calderon && has_encode) {
      # If, for some reason, only ENCODE was there (unlikely for Erythroblasts),
      cCRE_subset <- immune_cCREs %>%
        filter(
          class == class_i &
            (.data[[col_encode]] == TRUE | .data[[col_encode]] == "TRUE")
        )
    } else {
      # If neither column exists, skip this combination

      # (unlikely but could happen for Erythroblasts if there’s no Calderon col either)
      next
    }

    # Convert that subset to GRanges and find overlaps
    cCRE_subset_gr <- GRanges(
      seqnames = cCRE_subset$chr,
      ranges   = IRanges(start = cCRE_subset$start, end = cCRE_subset$end)
    )

    ov <- findOverlaps(genome_intervals_gr, cCRE_subset_gr, ignore.strand = TRUE)
    if (length(ov) > 0) {
      genome_intervals[[new_col]][queryHits(ov)] <- 1
    }
  }
}


# Double-check you have 21 new columns, one for each celltype_class combination
predictor_cols <- setdiff(
  colnames(genome_intervals),
  c("chromosome","start","end")  # your original columns
)

length(predictor_cols)

##########################################
#####      Running singular model    #####
##########################################
ihs_100kb_interval <- cbind(genome_intervals, ihs_100kb_interval)
ihs_100kb_interval <- ihs_100kb_interval[ , -c(25:27)]

run_logistic_regression_all_predictors <- function(results_df) {

  for (pop in unique_pops) {

    # Compute the 99th percentile threshold for that pop
    threshold <- quantile(ihs_100kb_interval[[pop]], 0.99, na.rm = TRUE)

    # Create a binary indicator of "high iHS" (1 if above 99th percentile)
    high_iHS_col <- paste0(pop, "_high_iHS")
    ihs_100kb_interval[[high_iHS_col]] <- ifelse(ihs_100kb_interval[[pop]] > threshold, 1, 0)

    # Construct formula: high_iHS_col ~ all 21 columns
    formula_str <- paste(high_iHS_col, "~", paste(predictor_cols, collapse=" + "))
    formula_obj <- as.formula(formula_str)

    # Fit logistic regression
    model <- glm(formula = formula_obj, data = ihs_100kb_interval, family = binomial)

    # Grab coefficients
    coef_summary <- summary(model)$coefficients

    # Try to get profile-likelihood CIs; fallback to Wald if error
    conf_intervals <- tryCatch({
      exp(confint(model))
    }, error = function(e) {
      message("Using Wald CIs for population ", pop, ": ", e$message)
      wald_lower <- exp(coef_summary[, "Estimate"] - 1.96 * coef_summary[, "Std. Error"])
      wald_upper <- exp(coef_summary[, "Estimate"] + 1.96 * coef_summary[, "Std. Error"])
      mat <- cbind(wald_lower, wald_upper)
      rownames(mat) <- rownames(coef_summary)
      mat
    })

    # Loop over each predictor
    for (predictor in predictor_cols) {
      if (predictor %in% rownames(coef_summary)) {
        log_or   <- coef_summary[predictor, "Estimate"]
        std_err  <- coef_summary[predictor, "Std. Error"]
        p_val    <- coef_summary[predictor, "Pr(>|z|)"]

        odds_ratio <- exp(log_or)
        lower_ci   <- conf_intervals[predictor, 1]
        upper_ci   <- conf_intervals[predictor, 2]

        results_df <- rbind(results_df, data.frame(
          Population    = pop,
          Predictor     = predictor,
          Log_Odds_Ratio= log_or,
          Odds_Ratio    = odds_ratio,
          Lower_CI      = lower_ci,
          Upper_CI      = upper_ci,
          Std_Error     = std_err,
          P_Value       = p_val,
          row.names     = NULL
        ))
      }
    }
  }
  return(results_df)
}

# Initialize df for storing results
results_all_predictors <- data.frame(
  Population     = character(),
  Predictor      = character(),
  Log_Odds_Ratio = numeric(),
  Odds_Ratio     = numeric(),
  Lower_CI       = numeric(),
  Upper_CI       = numeric(),
  Std_Error      = numeric(),
  P_Value        = numeric(),
  stringsAsFactors=FALSE
)

results_all_predictors <- run_logistic_regression_all_predictors(results_all_predictors)

##############################################
#####      Fixed effect meta-analysis    #####
##############################################

meta_summary <- data.frame(
  Predictor   = character(),
  FE_Estimate = numeric(),
  FE_CI_Lower = numeric(),
  FE_CI_Upper = numeric(),
  FE_P_Value  = numeric(),
  stringsAsFactors = FALSE
)

unique_predictors <- unique(results_all_predictors$Predictor)

for (pred in unique_predictors) {
  # Subset to that predictor across all populations
  subset_df <- subset(results_all_predictors, Predictor == pred,
                      select = c(Population, Log_Odds_Ratio, Std_Error))

  # Compute variance as square of Std_Error
  subset_df$Variance <- subset_df$Std_Error^2

  # Meta-analyze with fixed effects
  ma_result <- rma(
    yi     = Log_Odds_Ratio,
    vi     = Variance,
    data   = subset_df,
    method = "FE"
  )

  # Extract overall log(OR) estimate and p-value
  fe_estimate <- as.numeric(ma_result$b)
  fe_ci_lower <- as.numeric(ma_result$ci.lb)
  fe_ci_upper <- as.numeric(ma_result$ci.ub)
  fe_p_value  <- as.numeric(ma_result$pval)

  meta_summary <- rbind(meta_summary, data.frame(
    Predictor   = pred,
    FE_Estimate = fe_estimate,
    FE_CI_Lower = fe_ci_lower,
    FE_CI_Upper = fe_ci_upper,
    FE_P_Value  = fe_p_value,
    stringsAsFactors = FALSE
  ))
}

meta_summary <- meta_summary %>%
  mutate(
    Odds_Ratio  = exp(FE_Estimate),
    CI_Lower_OR = exp(FE_CI_Lower),
    CI_Upper_OR = exp(FE_CI_Upper)
  )


plot_df <- meta_summary %>%
  mutate(
    Odds_Ratio  = exp(FE_Estimate),
    CI_Lower_OR = exp(FE_CI_Lower),
    CI_Upper_OR = exp(FE_CI_Upper),
    SigFlag     = FE_P_Value < 0.05
  ) %>%
  mutate(
    cCRE_class = case_when(
      grepl("_pELS$",   Predictor) ~ "pELS",
      grepl("_dELS$",   Predictor) ~ "dELS",
      grepl("_PLS$",    Predictor) ~ "PLS",
      grepl("_CA_H3K4me3$", Predictor) ~ "CA_H3K4me3"
    ),
    cell_type = case_when(
      grepl("^Bcells_",       Predictor) ~ "Bcells",
      grepl("^CD4_Tcells_",   Predictor) ~ "CD4_Tcells",
      grepl("^CD8_Tcells_",   Predictor) ~ "CD8_Tcells",
      grepl("^Erythroblasts_",Predictor) ~ "Erythroblasts",
      grepl("^Myeloid_",      Predictor) ~ "Myeloid",
      grepl("^NK_",           Predictor) ~ "NK",
      grepl("^Progenitors_",  Predictor) ~ "Progenitors"
    )
  )

####################################################################
####    Manuscript Fig 5C: condensed cell-type specific model   ####
####################################################################
plot_df$p_stars <- cut(plot_df$FE_P_Value,
                       breaks = c(-Inf, 0.0005, 0.005, 0.05, Inf),
                       labels = c("***", "**", "*", ""),
                       right = FALSE)

desired_order <- c(
  "Progenitors", "Myeloid", "Erythroblasts",
  "NK", "CD8_Tcells", "CD4_Tcells", "Bcells"
)

plot_df <- plot_df %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = desired_order))

save.image("immune_cCREs_100kb_CTspecific_condensed_model_v2.RData")
```

**Step 6) Assess enrichment of eQTLs in windows thresholded at increasingly extreme iHS percentiles. Condition `${cond}` refers to the treatment of the cells (which only differs for the Harrison data: CTL, LPS, GARD; Natri is therefore it's own condition). These eQTLs are included in `${cond}/result_original.txt` for Harrison, or `Natri2022_S1_Permutation_significant_eQTLs.tsv` and `Natri2022_S3_eQTL_nominal.tsv` for Natri.**

`eQTL_enrichment_all_conditions.R` → Harrison: `${cond}_eQTL_enrichment.RData` and `Harrison_subset_for_5D.RData`, Natri: `enrichment_of_perm_sig_eQTLs.RData` and `Natri_subset_for_5D.RData`

```r
###############################################
####   Enrichment of Harrison‑et‑al. eQTL  ####
####       for GARD, CTL and LPS           ####
###############################################

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

## ------------------------------------------------------------------
##  USER paths
## ------------------------------------------------------------------
base_dir <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/eQTL_integration/Harrison"
ihs_dir  <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/inputs"
meta_file <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt"
setwd(base_dir)

## ------------------------------------------------------------------
##  helper: prepare merged eQTL results table for each condition, keeping unique SNPs only (the one with lowest FDR)
## ------------------------------------------------------------------
prepare_eqtl_table <- function(cond_dir) {
  eqtl_raw <- fread(file.path(cond_dir, "raw_results/result_original.txt"), drop = 1)
  ref      <- fread("SNP_positions.txt", drop = 1)

  eqtl_tbl <- merge(eqtl_raw, ref, by.x = "snps", by.y = "snp", all.x = TRUE)
  eqtl_tbl[, chr := as.numeric(sub("^chr", "", chr))]

  eqtl_tbl <- eqtl_tbl[order(FDR)][, .SD[1], by = snps]   # keep lowest FDR per SNP
  eqtl_tbl[, sig_eQTL := (FDR < 0.05)]
  eqtl_tbl
}

## ------------------------------------------------------------------
##  helper: enrichment for one population
## ------------------------------------------------------------------
process_population <- function(file_path, eqtl_tbl) {

  ihs <- fread(file_path,
               col.names = c("index","rsID","freq",
                             "iHH_0","iHH_1","iHS","std_iHS"))
  pop <- sub(".*/([A-Za-z]+)_.*", "\\1", file_path)

  ihs[, chr := as.numeric(str_extract(rsID, "(?<=chr)\\d+"))]
  ihs[, pos := as.numeric(str_extract(rsID, "(?<=pos)\\d+"))]

  common <- merge(ihs, eqtl_tbl, by = c("chr","pos"), all = FALSE)
  common[, abs_iHS := abs(std_iHS)]
  common  <- common[!is.na(abs_iHS)]

  pct_vec <- c(
    "50"   = as.numeric(quantile(common$abs_iHS, 0.50)),
    "70"   = as.numeric(quantile(common$abs_iHS, 0.70)),
    "90"   = as.numeric(quantile(common$abs_iHS, 0.90)),
    "95"   = as.numeric(quantile(common$abs_iHS, 0.95)),
    "99"   = as.numeric(quantile(common$abs_iHS, 0.99)),
    "99.5" = as.numeric(quantile(common$abs_iHS, 0.995))
  )

  enr_cut <- function(cut) {
    prop_sig  <- sum(common$abs_iHS > cut &  common$sig_eQTL) / sum(common$sig_eQTL)
    prop_nsig <- sum(common$abs_iHS > cut & !common$sig_eQTL) / sum(!common$sig_eQTL)
    c(prop_sig, prop_nsig, prop_sig / prop_nsig)
  }

  rbindlist(lapply(names(pct_vec), function(lbl) {
    v <- enr_cut(pct_vec[[lbl]])
    data.table(Population = pop,
               Percentile = lbl,
               Proportion_sig     = v[1],
               Proportion_non_sig = v[2],
               Enrichment_ratio   = v[3])
  }))
}

## ------------------------------------------------------------------
##  Run pipeline for each condition
## ------------------------------------------------------------------
conditions <- c("CTL", "LPS", "GARD")
for (cond in conditions) {

  message("Processing ", cond, " ...")

  # build unique‑SNP eQTL table
  cond_dir  <- file.path(base_dir, cond)
  eqtl_tbl  <- prepare_eqtl_table(cond_dir)

  # enrichment for all populations
  pop_files <- list.files(ihs_dir, pattern = "_iHS$", full.names = TRUE)
  all_populations_enrichment <- rbindlist(
    lapply(pop_files, process_population, eqtl_tbl = eqtl_tbl)
  )

  # save TSV
  fwrite(all_populations_enrichment,
         file = paste0("all_populations_enrichment_", cond, ".tsv"),
         sep  = "\t")

  # plot
  plot_enrichment(all_populations_enrichment,
                  cond,
                  paste0("enrichment_plot_", cond, ".png"))

  save.image(paste0(cond, "_eQTL_enrichment.RData"))
}

## ------------------------------------------------------------------
##  Extract TWA and KIG curves for manuscript
## ------------------------------------------------------------------
harrison_subset <- df_meta %>% filter(Population == "KIG" | Population == "TWA")
save(harrison_subset, file = "Harrison_subset_for_5D.RData")

###############################################
####     Enrichment of Natri et al. eQTL   ####
###############################################

## ------------------------------------------------------------------
##  USER paths + import files
## ------------------------------------------------------------------
base_dir <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/eQTL_integration/Natri"
ihs_dir  <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/inputs"
meta_file <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt"
setwd(base_dir)

# Import significant eQTLs
eqtl_tbl <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/eQTL_integration/Natri/Natri2022_S1_Permutation_significant_eQTLs.tsv", header=FALSE)
colnames(eqtl_tbl) <- c("gene_ID", "chr", "gene_pos", "n_variants_tested", "dist_between_gene_and_variant", "rsID", "snp_pos", "slope", "nominal_pval", "adjusted_pval")


# Import nominal eQTLs
nominal_tbl <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/eQTL_integration/Natri/Natri2022_S3_eQTL_nominal.tsv", header=FALSE)
colnames(nominal_tbl) <- c("gene_ID", "chr", "gene_pos", "n_variants_tested", "dist_between_gene_and_variant", "rsID", "snp_pos", "nominal_pval", "slope")
nominal_tbl_min <- nominal_tbl[order(nominal_pval)][, .SD[1], by = rsID] # retain entry with lowest nominal_pval per rsID (this is arbitrary tbh you just need to retain the unique rsIDs)

# Create table with the variants that are not eQTLs
non_eqtl_tbl <- nominal_tbl_min[!rsID %in% eqtl_tbl$rsID]

# Change the colnames for pos
setnames(eqtl_tbl, old = "snp_pos", new = "pos")
setnames(non_eqtl_tbl, old = "snp_pos", new = "pos")

## ------------------------------------------------------------------
##  Enrichment test
## ------------------------------------------------------------------
# Enrichment function for one population
process_population <- function(file_path, eqtl_tbl, non_eqtl_tbl) {
  ihs <- fread(file_path,
               col.names = c("index", "rsID", "freq", "iHH_0", "iHH_1", "iHS", "std_iHS"))
  pop <- sub(".*/([A-Za-z]+)_.*", "\\1", file_path)

  # Extract chr and pos from rsID
  ihs[, chr := as.numeric(str_extract(rsID, "(?<=chr)\\d+"))]
  ihs[, pos := as.numeric(str_extract(rsID, "(?<=pos)\\d+"))]

  # Merge significant and non-significant eQTLs with iHS
  sig_merged    <- as.data.table(merge(eqtl_tbl, ihs, by = c("chr", "pos"), all = FALSE))
  nonsig_merged <- as.data.table(merge(non_eqtl_tbl, ihs, by = c("chr", "pos"), all = FALSE))

  sig_merged[, abs_iHS := abs(std_iHS)]
  nonsig_merged[, abs_iHS := abs(std_iHS)]

  sig_merged <- sig_merged[!is.na(abs_iHS)]
  nonsig_merged <- nonsig_merged[!is.na(abs_iHS)]

  # Combined thresholds
  percentile_labels <- c("50", "70", "90", "95", "99", "99.5")
  percentile_cutoffs <- c(0.50, 0.70, 0.90, 0.95, 0.99, 0.995)

  pct_vec <- setNames(
    quantile(c(sig_merged$abs_iHS, nonsig_merged$abs_iHS), probs = percentile_cutoffs),
    percentile_labels
  )

  enr_cut <- function(cut) {
    prop_sig  <- mean(sig_merged$abs_iHS > cut)
    prop_nsig <- mean(nonsig_merged$abs_iHS > cut)
    c(prop_sig, prop_nsig, prop_sig / prop_nsig)
  }

  rbindlist(lapply(names(pct_vec), function(lbl) {
    v <- enr_cut(pct_vec[[lbl]])
    data.table(Population = pop,
               Percentile = lbl,
               Proportion_sig     = v[1],
               Proportion_non_sig = v[2],
               Enrichment_ratio   = v[3])
  }))
}


# Run for all the pops!
pop_files <- list.files(ihs_dir, pattern = "_iHS$", full.names = TRUE)

all_enrichment <- rbindlist(
  lapply(pop_files, process_population,
         eqtl_tbl = eqtl_tbl,
         non_eqtl_tbl = non_eqtl_tbl)
)

# Load metadata
metadata <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")

# Merge enrichment with metadata
df_meta <- merge(all_enrichment, metadata, by.x = "Population", by.y = "pop_code")

# Order populations by region
pop_region_mapping <- df_meta %>%
  distinct(Population, region) %>%
  mutate(region = factor(region, levels = c("South Asia","Southeast Asia",
                                            "East Asia","Oceania",
                                            "West Central Africa","East Africa",
                                            "Southern Africa","East Central Africa"))) %>%
  arrange(region)

df_meta$Population <- factor(df_meta$Population,
                             levels = pop_region_mapping$Population)
df_meta$region <- factor(df_meta$region,
                         levels = levels(pop_region_mapping$region))

# Set Percentile as ordered factor
df_meta$Percentile <- factor(df_meta$Percentile,
                             levels = c("50", "70", "90", "95", "99", "99.5"))


# Save results
fwrite(all_enrichment, "enrichment_sig_vs_nonsig_eQTLs.tsv", sep = "\t")

save.image("enrichment_of_perm_sig_eQTLs.RData")


## ------------------------------------------------------------------
##  Extract KWI, MTW, SMB curves for manuscript
## ------------------------------------------------------------------

natri_subset <- df_meta %>% filter(Population == "KWI" | Population == "MTW" | Population == "SMB")
save(natri_subset, file = "Natri_subset_for_5D.RData")
```

**Step 7) Get the distribution of MAFs for eQTLs and non-eQTLs for each eQTL dataset (Harrison: CTL, LPS, and GARD; and Natri). Also generate chi-squared enrichment p-values to create Supplementary Table 5.**

`eqtl_maf_bin_ihs.R ${cond}` → `eqtl_pop_maf_dist_${cond}.tsv` and `Harrison_Natri_chisq.tsv`

```r
library(dplyr)
library(data.table)
library(stringr)

args <- commandArgs(T)
data <- args[1]
load(data)

if (grepl("Natri", data)) {
  eqtl_tbl$sig_eQTL <- T
  non_eqtl_tbl$sig_eQTL <- F
  eqtl_tbl <- rbind(eqtl_tbl[, -10], non_eqtl_tbl)
  cond <- "natri"
}

pop_files <- list.files("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/inputs", pattern = "iHS", full.names = T)

# Get MAF for all eqtls
process_population <- function(file_path, eqtl_tbl, ihs_quantiles) {
  ihs <- fread(file_path,
               select = c(2, 3, 7),
               col.names = c("rsID", "freq", "std_iHS"))
  pop <- sub(".*/([A-Za-z]+)_.*", "\\1", file_path)

  ihs[, chr := as.numeric(str_extract(rsID, "(?<=chr)\\d+"))]
  ihs[, pos := as.numeric(str_extract(rsID, "(?<=pos)\\d+"))]

  common <- merge(ihs, eqtl_tbl, by = c("chr","pos"), all = FALSE)
  common[, abs_iHS := abs(std_iHS)]
  common <- common[!is.na(abs_iHS)]
  common$maf <- ifelse(common$freq > 0.5, 1 - common$freq, common$freq)

  # Summaries at diff cutoffs
  cuts <- quantile(common$abs_iHS, c(0.90, 0.95, 0.99))
  all_eqtl <- summary(common$maf[common$sig_eQTL])
  all_non_eqtl <- summary(common$maf[!common$sig_eQTL])
  ihs_0.90_eqtl <- summary(common$maf[common$sig_eQTL &
                                      common$abs_iHS > cuts[1]])
  ihs_0.90_non_eqtl <- summary(common$maf[!common$sig_eQTL &
                                          common$abs_iHS > cuts[1]])
  ihs_0.95_eqtl <- summary(common$maf[common$sig_eQTL &
                                      common$abs_iHS > cuts[2]])
  ihs_0.95_non_eqtl <- summary(common$maf[!common$sig_eQTL &
                                          common$abs_iHS > cuts[2]])
  ihs_0.99_eqtl <- summary(common$maf[common$sig_eQTL &
                                      common$abs_iHS > cuts[3]])
  ihs_0.99_non_eqtl <- summary(common$maf[!common$sig_eQTL &
                                          common$abs_iHS > cuts[3]])
  sums <- rbind(all_eqtl, all_non_eqtl, ihs_0.90_eqtl, ihs_0.90_non_eqtl, ihs_0.95_eqtl, ihs_0.95_non_eqtl, ihs_0.99_eqtl, ihs_0.99_non_eqtl)
  sums <- as.data.frame(sums)
  sums$group <- rownames(sums)
  sums$pop <- pop
  return(sums)
}

all_populations_enrichment <- rbindlist(
  lapply(pop_files, function(x) process_population(x, eqtl_tbl, qs))
)

# save TSV
fwrite(all_populations_enrichment,
       file = paste0("eqtl_pop_maf_dist_", cond, ".tsv"),
       sep  = "\t")

# Get chi-squared enrichments 
enr_cut <- function(common, cut) {
  ihs_sig <- common$abs_iHS > cut
  chisq <- chisq.test(table(ihs_sig, common$sig_eQTL), correct = F)
  print(cut)
  print(chisq$observed)
  print(chisq$stres)

  return(c(chisq$p.value, chisq$statistic))
}

process_population <- function(file_path, eqtl_tbl, ihs_quantiles) {
  ihs <- fread(file_path,
               select = c(2, 7),
               col.names = c("rsID", "std_iHS"))
  pop <- sub(".*/([A-Za-z]+)_.*", "\\1", file_path)
  print(pop)

  ihs[, chr := as.numeric(str_extract(rsID, "(?<=chr)\\d+"))]
  ihs[, pos := as.numeric(str_extract(rsID, "(?<=pos)\\d+"))]

  common <- merge(ihs, eqtl_tbl, by = c("chr","pos"), all = FALSE)
  common[, abs_iHS := abs(std_iHS)]
  common <- common[!is.na(abs_iHS)]

  pct_vec <- as.numeric(quantile(common$abs_iHS, ihs_quantiles))
  names(pct_vec) <- as.character(ihs_quantiles * 100)

  rbindlist(lapply(names(pct_vec), function(lbl) {
    v <- enr_cut(common, pct_vec[[lbl]])
    data.table(Population = pop,
               Percentile = lbl,
               chisq_p = v[1],
               chisq_stat = v[2])
  }))
}

qs <- c(0.50, 0.70, 0.90, 0.95, 0.99, 0.995)
all_populations_enrichment <- rbindlist(
  lapply(pop_files, function(x) process_population(x, eqtl_tbl, qs))
)

# save TSV
fwrite(all_populations_enrichment,
       file = paste0("all_pops_chisq_", cond, ".tsv"),

# Run in command line after getting all condition enrichments to subset to Supplementary Table 5
# head -n 1 all_pops_chisq_CTL.tsv > Harrison_Natri_chisq.tsv
# grep -E "KIG|TWA" all_pops_chisq_CTL.tsv >> Harrison_Natri_chisq.tsv
# grep -E "KWI|MTW|SMB" all_pops_chisq_natri.tsv >> Harrison_Natri_chisq.tsv
```

**Step 8) Plot the distribution of MAF for eQTLs and non-eQTLs only for the populations where the eQTLs were discovered to create Supplementary Figure 7.**

`ihs_shared_enrich.Rmd {eqtl maf dist}`

```r
library(data.table)
library(ggplot2)
library(tidyverse)

wd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/fst/ihs/"

# Supplementary Figure 7: eQTL MAF distributions
metadata <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
cond_files <- list.files(paste0(wd, "data"), pattern = "eqtl_pop_maf_dist", full.names = T)
conds <- c("CTL", "GARD", "LPS", "Natri")
eqtl_list <- lapply(cond_files, fread)
eqtl_list <- lapply(1:4, function(x) mutate(eqtl_list[[x]], condition = conds[x]))
df_meta <- bind_rows(eqtl_list)

df_meta$eqtl <- ifelse(grepl( "non", df_meta$group), "non-eQTL", "eQTL")
df_meta$cut <- case_when(
  grepl("all", df_meta$group) ~ 0,
  grepl("0.90", df_meta$group) ~ 0.90,
  grepl("0.95", df_meta$group) ~ 0.95,
  grepl("0.99", df_meta$group) ~ 0.99,
)
colnames(df_meta)[c(2,5)] <- c("q25", "q75")
df_meta <- left_join(df_meta, metadata, by = c("pop" = "pop_code"))

pdf(paste0(wd, "figures/supp_fig7_eqtl_maf_dist_subset.pdf"), height = 3, width = 5)
filter(df_meta, cut == 0 &
         (condition != "Natri" & pop %in% c("TWA", "KIG") |
          condition == "Natri" & pop %in% c("KWI", "MTW", "SMB"))) %>%
  ggplot(aes(x = pop, y = Median, color = eqtl, group = eqtl)) +
    geom_point(position = position_dodge(width = 0.8)) +
    geom_errorbar(aes(ymin = q25, ymax = q75), position = position_dodge(width = 0.8)) +
    scale_color_manual(values = c("tomato", "deepskyblue")) +
    facet_wrap(~condition, nrow = 1, scales = "free_x") +
    labs(y = "Median MAF", x = "", color = "") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          strip.background = element_blank())
dev.off()
```

**Step 9) Assemble all enrichment results from this section to create Figure 4.**

`mccre_plots.Rmd`

```r
library(ggplot2)
library(dplyr)
library(data.table)
library(khroma)
library(ggragged)
library(ggh4x)
library(cowplot)

muted <- color("muted")
colors <- c("dark grey", muted(4))

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/ihs/metabolic_ccres")

# Figure 4a: Single covariate logistic regression plots
metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
sample_numbers <- read.csv("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/samples_per_pop_baseline_filters_site_0.1_cohort_missing_nohet_SHAPEIT5phased_ancstate_nokin.txt", sep="")
metadata <- left_join(metadata, sample_numbers, by = c("pop_code" = "population"))
metadata$subsistence <- factor(metadata$subsistence, levels = c("Hunter-gatherer", "Agriculturalist", "Pastoralist"))

t1 <- fread("data/adipose-tissue_logistic_regression_ORs_CIs_pvals_all_pops_100kb.txt")
t1$tissue <- "adipose"
t2 <- fread("data/pancreas_logistic_regression_ORs_CIs_pvals_all_pops_100kb.txt")
t2$tissue <- "pancreas"
t3 <- fread("data/liver_logistic_regression_ORs_CIs_pvals_all_pops_100kb.txt")
t3$tissue <- "liver"
t4 <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/functional_annotations/immune_cCREs_logistic_regression_outputs_all_pops_100kb.txt")
t4$tissue <- "immune"

logistic_enrich_tissue <- rbind(t1, t2, t3, t4[,c(1,4:9)])
logistic_enrich_tissue <- left_join(metadata, logistic_enrich_tissue, by = c("pop_code" = "Population"))
logistic_enrich_tissue$significant <- ifelse(logistic_enrich_tissue$P_Value < 0.05, "1", "0")
logistic_enrich_tissue$sig_lab <- case_when(
  logistic_enrich_tissue$P_Value < 0.0005 ~ "***",
  logistic_enrich_tissue$P_Value < 0.005 ~ "**",
  logistic_enrich_tissue$P_Value < 0.05 ~ "*",
  T ~ "")
logistic_enrich_tissue$pop_code <- factor(logistic_enrich_tissue$pop_code, levels = metadata$pop_code[order(metadata$subsistence)])
logistic_enrich_tissue <- filter(logistic_enrich_tissue, element_type != "CA_H3K4me3" & n >= 15)

p4a <- ggplot(logistic_enrich_tissue,
              aes(x = Odds_Ratio, y = pop_code,
                  color = tissue, alpha = significant, shape = subsistence)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 2.5, ymax = 17.5),
            fill = "grey91", color = "grey91") +
  geom_point(fill = "white", size = 1) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.4) +
  geom_text(aes(label = sig_lab), nudge_y = 0.1, hjust = 0.5, size = 2, color = "black") +
  scale_color_manual(values = colors[2:5]) +
  scale_alpha_manual(values = c(0.3, 1)) +
  scale_shape_manual(values = c(25, 21, 22)) +
  labs(color = "", shape = "", x = "Odds ratio") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  facet_nested(tissue ~ element_type) +
  guides(alpha = "none", color = "none", fill = "none",
         shape = guide_legend(override.aes = list(fill = "black"))) +
  theme_classic(base_size = 7) +
  theme(legend.position = "top",
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 5),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_blank(),
        panel.spacing = unit(1, "lines"))
p4a
saveRDS(p4a, "figures/fig_4a.rds")

# Figure 4b: Multiregression logistic regression plots
load("data/fe_test.RData")

names <- strsplit("immune_pELS + immune_dELS + immune_PLS + brain_pELS + brain_dELS + brain_PLS + adipose_pELS + adipose_dELS + adipose_PLS + pancreas_pELS + pancreas_dELS + pancreas_PLS + liver_pELS + liver_dELS + liver_PLS", " [+] ") %>% unlist()

fe_df_all <- lapply(fe_all_tissues, function(x) data.frame(
  Estimate = x$b,
  CI_Lower = x$ci.lb,
  CI_Upper = x$ci.ub,
  P_Value = x$pval
))
fe_df_all <- bind_rows(fe_df_all)
fe_df_all$Element <- names
fe_df_all$formula <- "all"

fe_df_sin <- lapply(fe_single_tissue, function(x) data.frame(
  Estimate = x$b,
  CI_Lower = x$ci.lb,
  CI_Upper = x$ci.ub,
  P_Value = x$pval
))
fe_df_sin <- bind_rows(fe_df_sin)
fe_df_sin$Element <- names
fe_df_sin$formula <- "single"

fe_df <- rbind(fe_df_all, fe_df_sin)
fe_df$tissue <- lapply(fe_df$Element, function(x) strsplit(x, "_")[[1]][1]) %>% unlist()
fe_df$ccre <- lapply(fe_df$Element, function(x) strsplit(x, "_")[[1]][2]) %>% unlist()
fe_df$significant <- ifelse(fe_df$P_Value < 0.05, "1", "0")
fe_df <- fe_df %>%
  dplyr::mutate(
    Odds_Ratio = exp(Estimate),
    CI_Lower_OR = exp(CI_Lower),
    CI_Upper_OR = exp(CI_Upper))

fwrite(fe_df, "data/all_tissues_fe.txt")

# Plot for figure 4b
fe_df <- fread("data/all_tissues_fe.txt")
fe_df$tissue <- factor(fe_df$tissue, levels = c("brain", "immune", "adipose", "pancreas", "liver"))
fe_df$sig_lab <- factor(case_when(fe_df$P_Value < 0.0005 ~ "***", fe_df$P_Value < 0.005 ~ "**", fe_df$P_Value < 0.05 ~ "*", T ~ ""))
fe_df$ccre <- factor(fe_df$ccre)

p4b <- filter(fe_df, formula == "all") %>%
  ggplot(aes(x = Odds_Ratio, y = ccre, color = tissue, group = tissue, alpha = as.factor(significant), label = sig_lab)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = CI_Lower_OR, xmax = CI_Upper_OR), height = 0.4, position = position_dodge(width = 0.6)) +
  geom_text(position = position_dodge(width = 0.6), vjust = -0.025, size = 3, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_colour_manual(values = colors) +
  scale_alpha_manual(values = c(0.3, 1)) +
  labs(x = "Odds ratio", y = "", color = "") +
  guides(alpha = "none") +
  theme_classic(base_size = 7) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold"),
        panel.spacing = unit(0.5, "lines"))
p4b
saveRDS(p4b, "figures/fig_4b.rds")

# Figure 4c: Celltype-specific immune cCRE enrichments plot
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/functional_annotations/immune_cCREs_100kb_CTspecific_condensed_model.RData")

fwrite(plot_df, "data/plot_4c.txt")
plot_df <- fread("data/plot_4c.txt")
plot_df$cell_type <- case_match(plot_df$cell_type,
                                "Bcells" ~ "B cells",
                                "CD4_Tcells" ~ "CD4 T cells",
                                "CD8_Tcells" ~ "CD8 T cells",
                                "Myeloid" ~ "Myeloid cells",
                                "NK" ~ "NK cells",
                                "Erythroblasts" ~ "Erythroblasts",
                                "Progenitors" ~ "Progenitors")
plot_df$cCRE_class <- factor(plot_df$cCRE_class, levels = c("PLS", "pELS", "dELS"))

p4c <- plot_df %>%
  ggplot(aes(x = Odds_Ratio, y = cell_type, alpha = SigFlag, label = p_stars)) +
  geom_point(size = 1, color = "#2179B4") +
  geom_errorbarh(aes(xmin = CI_Lower_OR, xmax = CI_Upper_OR), height = 0.2, color = "#2179B4") +
  geom_text(nudge_y = 0.1, hjust = 0.5, size = 3, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  facet_wrap(~ cCRE_class, ncol = 3) +
  labs(x = "Odds ratio", y = "Cell type") +
  scale_alpha_manual(values = c(0.3, 1), guide = "none") +
  theme_classic(base_size = 7) +
  theme(
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    #panel.spacing = unit(0.8, "lines"),
    #panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.8)
  )
p4c
saveRDS(p4c, "figures/fig_4c.rds")

# Figure 4d: eQTL enrichments plot
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/eQTL_integration/Natri/Natri_subset_for_5D.RData")
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/eQTL_integration/Harrison/Harrison_subset_for_5D.RData")
fig_5d <- rbind(natri_subset, harrison_subset, fill = TRUE)

# Recode "CTL" and NA to unstimulated
fig_5d$condition <- ifelse(is.na(fig_5d$condition) | fig_5d$condition == "CTL",
                           "unstimulated", fig_5d$condition)
fig_5d$condition <- factor(fig_5d$condition, levels = c("unstimulated", "GARD", "LPS"))

# Drop 99.5 percentile to be consistent with other figures
fig_5d <- fig_5d %>%
  arrange(country, Population) %>%
  mutate(
    country = factor(country, levels = c("Indonesia", "Uganda")),
    Population = factor(Population, levels = c("KWI", "MTW", "SMB", "TWA", "KIG")),
    subsistence = factor(subsistence, levels = c("Hunter-gatherer", "Agriculturalist"))
  )

p4d <- filter(fig_5d, Percentile != "99.5" & condition == "unstimulated") %>%
  ggplot(aes(x = Percentile, y = Enrichment_ratio,
             color = subsistence, fill = subsistence, shape = subsistence, group = Population)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_line(show.legend = F) +
  geom_point(na.rm = TRUE, size = 2) +
  facet_ragged_rows(rows = vars(country), cols = vars(Population), axes = "all_x") +
  scale_discrete_manual(c("color", "fill"), values = c("#C21F84", "#FDAE61")) +
  scale_shape_manual(values = c(25, 21)) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2),
    limits = c(min(fig_5d$Enrichment_ratio, na.rm = TRUE),
               max(fig_5d$Enrichment_ratio, na.rm = TRUE))) +
  labs(x = "iHS percentile",
    y = "Enrichment ratio",
    fill = "", color = "") +
  guides(shape = "none", fill = "none",
         color = guide_legend(override.aes = list(fill = c("#C21F84", "#FDAE61"), shape = c(25, 21)))) +
  theme_classic(base_size = 7) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        legend.position = "top")
p4d
saveRDS(p4d, "figures/fig_4d.rds")

# Arrange all plots
p4a <- readRDS("figures/fig_4a.rds")
p4b <- readRDS("figures/fig_4b.rds")
p4c <- readRDS("figures/fig_4c.rds")
p4d <- readRDS("figures/fig_4d.rds")

ccre_figures <- plot_grid(
    plot_grid(p4a, labels = c("a"), label_size = 8),
    plot_grid(p4b, p4c, p4d, labels = c("b", "c", "d"),
              rel_widths = c(1, 0.8), nrow = 3, label_size = 8),
  rel_widths = c(1, 1), ncol = 2)

ccre_figures
save_plot("figures/ccre_figures.pdf", ccre_figures, base_width = 7, base_height = 9)
```
