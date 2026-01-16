# VIP and metabolic enrichment analysis

Working directory: `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/fst/ihs` 

## iHS and Tajima’s enrichments for VIP and metabolic genes

---

```r
# library(biomaRt)
library(data.table)
library(ggplot2)
library(ggragged)
library(tidyverse)

wd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/fst/ihs/"
#load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/freeze2_iHS_per_gene_all_pops_BCpipeline.RData")
metadata <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
metadata$region <- factor(metadata$region, levels = c("West Central Africa", "East Central Africa", "East Africa", "Southern Africa", "South Asia", "East Asia", "Southeast Asia", "Oceania"))
metadata$subsistence <- factor(metadata$subsistence, levels = c("Hunter-gatherer", "Agriculturalist", "Pastoralist"))
VIP_list <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze1/metadata/from_providers/viral_effects_file_april2023.txt", header = FALSE)
# metabolic_genes <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/curate_genes/KEGG_09100_Metabolism/09100_METABOLISM_ALLGENES.txt", header = F)
# metabolic_genes <- sapply(unique(metabolic_genes$V2),  function(x) strsplit(x, ":")[[1]][2])
# ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# metabolic_ensg <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id'),
#       filters = 'entrezgene_id',
#       values = metabolic_genes,
#       mart = ensembl)
# saveRDS(metabolic_ensg, paste0(wd, "metabolic_ensg.rds")) # Doesn't map 1:1, 1679 unique Entrez/KEGG IDs vs. 1858 unique Ensembl IDs
metabolic_ensg <- readRDS(paste0(wd, "metabolic_ensg.rds"))

# Format per-population iHS results and calculate percentile per-gene
# gff3_results <- cbind(gff3, gff3_ihs)
# names(gff3_results)[14:47] <- metadata$pop_code
# gff3_results <- gff3_results %>% gather(pop, ihs, 14:47)
# gff3_results <- merge(gff3_results, metadata, by.x = "pop", by.y = "pop_code", all.x = TRUE)
# gff3_results$type <- ifelse(gff3_results$gene_IDs %in% VIP_list$V1, "VIP", "non-VIP")
# gff3_results$metabolic <- ifelse(gff3_results$gene_IDs %in% metabolic_ensg$ensembl_gene_id, "metabolic", "non-metabolic")
# gff3_results <- gff3_results %>%
#   group_by(pop) %>%
#   mutate(within_pop_percentile = if_else(is.na(ihs),
#                                          as.numeric(NA),
#                                          rank(ihs, ties.method = "average") / sum(!is.na(ihs)) * 100))
# saveRDS(gff3_results, paste0(wd, "ihs_gff3_result.rds"))

gff3_results <- readRDS(paste0(wd, "ihs_gff3_result.rds"))

# Format percentiles to plot proportion of genes that in a category (VIP, metabolic) at certain cutoffs
prop_pct <- function(percentiles, test, category, cutoff) {
  # Get the smaller tail of dist
  if(cutoff < 50) {
    tmp <- percentiles %>% filter(within_pop_percentile <= cutoff) %>% group_by(pop, get(category)) %>% count()
  } else {
    tmp <- percentiles %>% filter(within_pop_percentile >= cutoff) %>% group_by(pop, get(category)) %>% count()
  }
  colnames(tmp)[3] <- "n_category"
  tmp <- percentiles %>% filter(!is.na(get(test))) %>% group_by(pop, get(category)) %>% count() %>% left_join(tmp, by = c("pop", "get(category)"))
  colnames(tmp) <- c("pop", category, "n", "n_category")
  tmp[[as.character(cutoff)]] <- tmp$n_category / tmp$n
  return(tmp[,c(1,2,5)])
}

# Convert proportion of percentiles dataframe to enrichment ratio
pct_to_er <- function(percentiles, category, enrich_str) {
  tmp_list <- lapply(unique(percentiles$pop), function(x) data.frame(
    pop = rep(x, length(unique(percentiles$percentile))),
    percentiles = unique(percentiles$percentile),
    enrich_ratio = percentiles$proportion[percentiles$pop == x & percentiles[[`category`]] == enrich_str] / percentiles$proportion[percentiles$pop == x & percentiles[[`category`]] != enrich_str]))
  tmp <- bind_rows(tmp_list)
  tmp <- tmp %>% left_join(unique(percentiles[,c(1,3,4)]), by = "pop")
  return(tmp)
}

cutoffs <- c(50, 95, 98, 99, 99.5)

# Get enrichment ratios for VIPs and metabolic genes
prop_vip_list <- lapply(cutoffs, prop_pct, percentiles = gff3_results, test = "ihs", category = "type")
prop_vip_list <- purrr::reduce(prop_vip_list, left_join, by = c("pop", "type")) %>% left_join(metadata[,c(1,3,5)], by = c("pop" = "pop_code")) %>% pivot_longer(cols = as.character(cutoffs), names_to = "percentile", values_to = "proportion")
prop_vip <- pct_to_er(prop_vip_list, "type", "VIP")

prop_metabolic_list <- lapply(cutoffs, prop_pct, percentiles = gff3_results, test = "ihs", category = "metabolic")
prop_metabolic_list <- purrr::reduce(prop_metabolic_list, left_join, by = c("pop", "metabolic")) %>% left_join(metadata[,c(1,3,5)], by = c("pop" = "pop_code")) %>% pivot_longer(cols = as.character(cutoffs), names_to = "percentile", values_to = "proportion")
prop_metabolic <- pct_to_er(prop_metabolic_list, "metabolic", "metabolic")

pdf(paste0(wd, "ihs_vip_metabolic_enrich.pdf"), height = 10, width = 10)
# Plot enriched proportions
prop_vip %>% ggplot(aes(x = percentiles, y = enrich_ratio, color = subsistence, group = 1)) +
  geom_hline(yintercept = 1, color = "black", linetype = 2) +
  geom_line() +
  geom_point(na.rm = TRUE) +
  facet_ragged_rows(vars(region), vars(pop), axes = "all_x") +
  theme_classic() +
  labs(title = "VIP genes iHS enrichment", x = "Percentile", y = "Enrichment Ratio", color = "Subsistence") +
  theme(strip.text.y.right = element_text(angle = 0),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))

prop_metabolic %>% ggplot(aes(x = percentiles, y = enrich_ratio, color = subsistence, group = 1)) +
  geom_hline(yintercept = 1, color = "black", linetype = 2) +
  geom_line() +
  geom_point(na.rm = TRUE) +
  facet_ragged_rows(vars(region), vars(pop), axes = "all_x") +
  theme_classic() +
  labs(title = "Metabolic genes iHS enrichment", x = "Percentile", y = "Enrichment Ratio", color = "Subsistence") +
  theme(strip.text.y.right = element_text(angle = 0),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
dev.off()

# Process Tajima's gene-wise data
# tajima_files <- list.files("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/fst/thetas/", pattern = "gene_quantile", full.names = T)
# tajimas_list <- lapply(tajima_files, fread)
# names(tajimas_list) <- metadata$pop_code
# tajimas <- bind_rows(tajimas_list, .id = "pop")
# # gene_data <- readRDS("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/gs_data/hg19_genes.rds")
#
# tajimas <- tajimas %>% left_join(unique(gff3_results[,c("gene_IDs", "names")]), by = c("gene_name" = "names")) %>% na.omit()
# tajimas$type <- ifelse(tajimas$gene_IDs %in% VIP_list$V1, "VIP", "non-VIP")
# tajimas$metabolic <- ifelse(tajimas$gene_IDs %in% metabolic_ensg$ensembl_gene_id, "metabolic", "non-metabolic")
# colnames(tajimas) <- c("pop", "names", "tajima", "within_pop_percentile", "gene_IDs", "type", "metabolic")
# tajimas$within_pop_percentile <- tajimas$within_pop_percentile * 100
# saveRDS(tajimas, paste0(wd, "tajimas.rds"))
tajimas <- readRDS(paste0(wd, "tajimas.rds"))

cutoffs <- c(0.05, 1, 5, 50, 95, 99, 99.5)

# Get enrichment ratios for VIPs and metabolic genes
tajima_vip_list <- lapply(cutoffs, prop_pct, percentiles = tajimas, test = "tajima", category = "type")
tajima_vip_list <- purrr::reduce(tajima_vip_list, left_join, by = c("pop", "type")) %>% left_join(metadata[,c(1,3,5)], by = c("pop" = "pop_code")) %>% pivot_longer(cols = as.character(cutoffs), names_to = "percentile", values_to = "proportion")
tajima_vip <- pct_to_er(tajima_vip_list, "type", "VIP")

tajima_metabolic_list <- lapply(cutoffs, prop_pct, percentiles = tajimas, test = "tajima", category = "metabolic")
tajima_metabolic_list <- purrr::reduce(tajima_metabolic_list, left_join, by = c("pop", "metabolic")) %>% left_join(metadata[,c(1,3,5)], by = c("pop" = "pop_code")) %>% pivot_longer(cols = as.character(cutoffs), names_to = "percentile", values_to = "proportion")
tajima_metabolic <- pct_to_er(tajima_metabolic_list, "metabolic", "metabolic")

pdf(paste0(wd, "tajima_vip_metabolic_enrich.pdf"), height = 10, width = 10)
# Plot enriched proportions
tajima_vip %>% ggplot(aes(x = percentiles, y = log2(enrich_ratio), color = subsistence, group = 1)) +
  geom_hline(yintercept = 0, color = "black", linetype = 2) +
  geom_line() +
  geom_point(na.rm = TRUE) +
  facet_ragged_rows(vars(region), vars(pop), axes = "all_x") +
  theme_classic() +
  labs(title = "VIP genes Tajima's enrichment", x = "Percentile", y = expression(paste("log"[2], " Enrichment Ratio")), color = "Subsistence") +
  theme(strip.text.y.right = element_text(angle = 0),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))

tajima_metabolic %>% ggplot(aes(x = percentiles, y = log2(enrich_ratio), color = subsistence, group = 1)) +
  geom_hline(yintercept = 0, color = "black", linetype = 2) +
  geom_line() +
  geom_point(na.rm = TRUE) +
  facet_ragged_rows(vars(region), vars(pop), axes = "all_x") +
  theme_classic() +
  labs(title = "Metabolic genes Tajima's enrichment", x = "Percentile", y = expression(paste("log"[2], " Enrichment Ratio")), color = "Subsistence") +
  theme(strip.text.y.right = element_text(angle = 0),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
dev.off()

```

## Bridget’s Code

---

```r
##############################################
######            READ ME           ##########
##############################################

# start with outputs from freeze2_iHS_per_gene_all_pops_BCpipeline.R 

# the gff3_ihs matrix has the standardized iHS value per gene for all the populations

# this script includes some cursory analysis of VIP enrichment

library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/freeze2_iHS_per_gene_all_pops_BCpipeline.RData")

#########################################################################################
#####       Bind results from the last step's gff3 df and reshape to long form      #####
#########################################################################################

# clean up objects from before
gc()
rm(currentFile)

# gff3_ihs has the per-population results
gff3_results <- cbind(gff3, gff3_ihs)

# rename columns 14 - 47
names(gff3_results)[14:47] <- unique_pops
gff3_results <- gff3_results %>% gather(pop, ihs, 14:47)

# combine gff3_results with metadata 
metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
gff3_results <- merge(gff3_results, metadata, by.x = "pop", by.y = "pop_code", all.x = TRUE)
unique(gff3_results$region)

# gff3_results <- na.omit(gff3_results) 

##################################
#####       Exploratory      #####
##################################

# genome-wide distributions of iHS values per population, separated by regions-------------
gff3_results %>% filter(region == "Southeast Asia" | region == "East Asia") %>% ggplot(aes(x = pop, y = ihs)) + 
  ggdist::stat_halfeye(aes(fill = pop, fill = after_scale(colorspace::lighten(fill, .5))),
                       position = position_nudge(x = -.3), interval_size = 1) +
  guides(fill="none") + ggtitle("Asian Populations") + theme(plot.title = element_text(hjust = 0.5))

gff3_results %>% filter(region == "West Central Africa" | region == "East Africa" | region == "Southern Africa" | region == "East Central Africa") %>% ggplot(aes(x = pop, y = ihs)) + 
  ggdist::stat_halfeye(aes(fill = pop, fill = after_scale(colorspace::lighten(fill, .5))),
                       position = position_nudge(x = -.3), interval_size = 1) +
  guides(fill="none") + ggtitle("African Populations") + theme(plot.title = element_text(hjust = 0.5))

gff3_results %>% filter(region == "Oceania") %>% ggplot(aes(x = pop, y = ihs)) + 
  ggdist::stat_halfeye(aes(fill = pop, fill = after_scale(colorspace::lighten(fill, .5))),
                       position = position_nudge(x = -.3), interval_size = 1) +
  guides(fill="none") + ggtitle("Oceanian Populations") + theme(plot.title = element_text(hjust = 0.5))

################################################
#####       Compare VIPs and non-VIPs      #####
################################################
VIP_list <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze1/metadata/from_providers/viral_effects_file_april2023.txt", header = FALSE)
gff3_results$type <- ifelse(gff3_results$gene_IDs %in% VIP_list$V1, "VIP", "non-VIP")

# box plot showing per-population distributions of the per-gene average iHS for VIPs vs. non-VIPS
ggplot(gff3_results, aes(x = pop, y = ihs, fill = type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_minimal() +
  labs(title = "iHS Scores by Population for VIPs vs. non-VIPs (Freeze 2)",
       x = "Population",
       y = "iHS",
       fill = "Type") +
  theme(plot.title = element_text(hjust = 0.5))

mean_iHS_comp <- gff3_results %>% group_by(pop, type) %>% summarise(mean_iHS=mean(ihs, na.rm = TRUE))

# histogram showing the mean genome-wide iHS scores for VIPs vs. non-VIPS per population
ggplot(mean_iHS_comp, aes(x = pop, y = mean_iHS, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Mean iHS for VIPs vs.non-VIPs by Population",
       x = "Population",
       y = "Mean iHS",
       fill = "Type") +
  theme(plot.title = element_text(hjust = 0.5))

# counts of genes that have iHS value (non-NAs) per population
n_genes_with_iHS <- gff3_results %>%
  filter(!is.na(ihs)) %>%
  group_by(pop, type) %>%
  count()

# find the percentile of each iHS value (gene-level) within each population
percentiles <- gff3_results %>%
  group_by(pop) %>%
  mutate(within_pop_percentile = if_else(is.na(ihs), 
                                         as.numeric(NA), 
                                         rank(ihs, ties.method = "average") / sum(!is.na(ihs)) * 100))

# save dataframe with percentiles
save(percentiles, file="freeze2_per_gene_iHS_inc_VIP_and_percentile_info_by_pop.RData")

###############################################################################
#### Inspecting outliers: filter for 95th percentile and above iHS scores #####
###############################################################################

percentile_95 <- percentiles %>% filter(within_pop_percentile >= 95)

ggplot(percentile_95, aes(x = pop, y = ihs, fill = type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_minimal() +
  labs(title = "Boxplot of iHS Scores above 95th percentile by Population for VIPs vs. non-VIPs",
       x = "Population",
       y = "iHS",
       fill = "Type") +
  theme(plot.title = element_text(hjust = 0.5))

# raw counts of VIPs vs non VIPs in the 95th percentile
prop_VIP_in_percentile_95 <- percentile_95 %>% group_by(pop, type) %>% count()
colnames(prop_VIP_in_percentile_95) <- c("pop", "type", "n_in_95th_percentile")

# add column of total genes in the 95th percentile
prop_VIP_in_percentile_95 <- prop_VIP_in_percentile_95 %>%
  group_by(pop) %>%
  mutate(total_n_in_95th_percentile = sum(n_in_95th_percentile, na.rm = TRUE))

# calculate prop of VIPs vs non VIPS in the 95th percentile
prop_VIP_in_percentile_95$prop_in_95th_percentile <- prop_VIP_in_percentile_95$n_in_95th_percentile / prop_VIP_in_percentile_95$total_n_in_95th_percentile

# add in counts of VIPs vs non VIPS that have iHS values
prop_VIP_in_percentile_95$n_with_iHS <- n_genes_with_iHS$n

# calculate prop of VIPs in the 95th percentile out of all VIPs with iHS values
prop_VIP_in_percentile_95$prop_in_95th_percentile_among_type_with_iHS <- prop_VIP_in_percentile_95$n_in_95th_percentile / prop_VIP_in_percentile_95$n_with_iHS

#############################################################################
##### plot the proportions of VIPs vs. non-VIPs in the 95th percentile  #####
#############################################################################

# first merge with metadata
prop_VIP_in_percentile_95 <- merge(prop_VIP_in_percentile_95, metadata, by.x = "pop", by.y = "pop_code", all.x = TRUE)

prop_VIP_in_percentile_95 %>%
  filter(region %in% c("Southeast Asia", "East Asia", "South Asia")) %>% ggplot(aes(x = type, y = prop_in_95th_percentile, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.03), width = 0.7) + # Adjust bar width and distance between bars here
  facet_wrap(~ pop, scales = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(y = "Proportion in 95th Percentile", fill = "Type") +
  ggtitle("Proportion of VIPs vs. non-VIPs among genes in the 95th percentile (Asian populations)")

prop_VIP_in_percentile_95 %>%
  filter(region == "West Central Africa" | region == "East Africa" | region == "Southern Africa" | region == "East Central Africa") %>% ggplot(aes(x = type, y = prop_in_95th_percentile, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.03), width = 0.7) + # Adjust bar width and distance between bars here
  facet_wrap(~ pop, scales = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(y = "Proportion in 95th Percentile", fill = "Type") +
  ggtitle("Proportion of VIPs vs. non-VIPs among genes in the 95th percentile (African populations)")

prop_VIP_in_percentile_95 %>%
  filter(region == "Oceania") %>% ggplot(aes(x = type, y = prop_in_95th_percentile, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.03), width = 0.7) + # Adjust bar width and distance between bars here
  facet_wrap(~ pop, scales = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(y = "Proportion in 95th Percentile", fill = "Type") +
  ggtitle("Proportion of VIPs vs. non-VIPs among genes in the 95th percentile (Oceanian populations)")

################################################################################################################
##### plot the proportion of VIPs vs. non-VIPs out of all genes in their own category that have iHS values #####
################################################################################################################

prop_VIP_in_percentile_95 %>%
  filter(region %in% c("Southeast Asia", "East Asia", "South Asia")) %>% ggplot(aes(x = type, y = prop_in_95th_percentile_among_type_with_iHS, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.03), width = 0.7) + # Adjust bar width and distance between bars here
  facet_wrap(~ pop, scales = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(y = "Proportion", fill = "Type") +
  ggtitle("Proportion of VIPs and non-VIPs appearing in the 95th percentile among genes of the same type (Asian populations)")

prop_VIP_in_percentile_95 %>%
  filter(region == "West Central Africa" | region == "East Africa" | region == "Southern Africa" | region == "East Central Africa") %>% ggplot(aes(x = type, y = prop_in_95th_percentile_among_type_with_iHS, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.03), width = 0.7) + # Adjust bar width and distance between bars here
  facet_wrap(~ pop, scales = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(y = "Proportion", fill = "Type") +
  ggtitle("Proportion of VIPs and non-VIPs appearing in the 95th percentile among genes of the same type (African populations)")

prop_VIP_in_percentile_95 %>%
  filter(region == "Oceania") %>% ggplot(aes(x = type, y = prop_in_95th_percentile_among_type_with_iHS, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.03), width = 0.7) + # Adjust bar width and distance between bars here
  facet_wrap(~ pop, scales = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(y = "Proportion", fill = "Type") +
  ggtitle("Proportion of VIPs and non-VIPs appearing in the 95th percentile among genes of the same type (Oceanian populations)")

###############################################################################
#### Inspecting outliers: filter for 90th percentile and above iHS scores #####
###############################################################################

percentile_90 <- percentiles %>% filter(within_pop_percentile >= 90)

ggplot(percentile_90, aes(x = pop, y = ihs, fill = type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_minimal() +
  labs(title = "Boxplot of iHS Scores above 90th percentile by Population for VIPs vs. non-VIPs",
       x = "Population",
       y = "iHS",
       fill = "Type") +
  theme(plot.title = element_text(hjust = 0.5))

# raw counts of VIPs vs non VIPs in the 90th percentile
prop_VIP_in_percentile_90 <- percentile_90 %>% group_by(pop, type) %>% count()
colnames(prop_VIP_in_percentile_90) <- c("pop", "type", "n_in_90th_percentile")

# add column of total genes in the 90th percentile
prop_VIP_in_percentile_90 <- prop_VIP_in_percentile_90 %>%
  group_by(pop) %>%
  mutate(total_n_in_90th_percentile = sum(n_in_90th_percentile, na.rm = TRUE))

# calculate prop of VIPs vs non VIPS in the 90th percentile
prop_VIP_in_percentile_90$prop_in_90th_percentile <- prop_VIP_in_percentile_90$n_in_90th_percentile / prop_VIP_in_percentile_90$total_n_in_90th_percentile

# add in counts of VIPs vs non VIPS that have iHS values
prop_VIP_in_percentile_90$n_with_iHS <- n_genes_with_iHS$n

# calculate prop of VIPs in the 90th percentile out of all VIPs with iHS values
prop_VIP_in_percentile_90$prop_in_90th_percentile_among_type_with_iHS <- prop_VIP_in_percentile_90$n_in_90th_percentile / prop_VIP_in_percentile_90$n_with_iHS

###############################################################################
#### Inspecting outliers: filter for 97th percentile and above iHS scores #####
###############################################################################

percentile_97 <- percentiles %>% filter(within_pop_percentile >= 97)

ggplot(percentile_97, aes(x = pop, y = ihs, fill = type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_minimal() +
  labs(title = "Boxplot of iHS Scores above 97th percentile by Population for VIPs vs. non-VIPs",
       x = "Population",
       y = "iHS",
       fill = "Type") +
  theme(plot.title = element_text(hjust = 0.5))

# raw counts of VIPs vs non VIPs in the 97th percentile
prop_VIP_in_percentile_97 <- percentile_97 %>% group_by(pop, type) %>% count()
colnames(prop_VIP_in_percentile_97) <- c("pop", "type", "n_in_97th_percentile")

# add column of total genes in the 97th percentile
prop_VIP_in_percentile_97 <- prop_VIP_in_percentile_97 %>%
  group_by(pop) %>%
  mutate(total_n_in_97th_percentile = sum(n_in_97th_percentile, na.rm = TRUE))

# calculate prop of VIPs vs non VIPS in the 97th percentile
prop_VIP_in_percentile_97$prop_in_97th_percentile <- prop_VIP_in_percentile_97$n_in_97th_percentile / prop_VIP_in_percentile_97$total_n_in_97th_percentile

# add in counts of VIPs vs non VIPS that have iHS values
prop_VIP_in_percentile_97$n_with_iHS <- n_genes_with_iHS$n

# calculate prop of VIPs in the 97th percentile out of all VIPs with iHS values
prop_VIP_in_percentile_97$prop_in_97th_percentile_among_type_with_iHS <- prop_VIP_in_percentile_97$n_in_97th_percentile / prop_VIP_in_percentile_97$n_with_iHS

###############################################################################
#### Inspecting outliers: filter for 98th percentile and above iHS scores #####
###############################################################################
percentile_98 <- percentiles %>% filter(within_pop_percentile >= 98)

ggplot(percentile_98, aes(x = pop, y = ihs, fill = type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_minimal() +
  labs(title = "Boxplot of iHS Scores above 98th percentile by Population for VIPs vs. non-VIPs",
       x = "Population",
       y = "iHS",
       fill = "Type") +
  theme(plot.title = element_text(hjust = 0.5))

# raw counts of VIPs vs non VIPs in the 98th percentile
prop_VIP_in_percentile_98 <- percentile_98 %>% group_by(pop, type) %>% count()
colnames(prop_VIP_in_percentile_98) <- c("pop", "type", "n_in_98th_percentile")

# add column of total genes in the 98th percentile
prop_VIP_in_percentile_98 <- prop_VIP_in_percentile_98 %>%
  group_by(pop) %>%
  mutate(total_n_in_98th_percentile = sum(n_in_98th_percentile, na.rm = TRUE))

# calculate prop of VIPs vs non VIPS in the 98th percentile
prop_VIP_in_percentile_98$prop_in_98th_percentile <- prop_VIP_in_percentile_98$n_in_98th_percentile / prop_VIP_in_percentile_98$total_n_in_98th_percentile

# add in counts of VIPs vs non VIPS that have iHS values
prop_VIP_in_percentile_98$n_with_iHS <- n_genes_with_iHS$n

# calculate prop of VIPs in the 98th percentile out of all VIPs with iHS values
prop_VIP_in_percentile_98$prop_in_98th_percentile_among_type_with_iHS <- prop_VIP_in_percentile_98$n_in_98th_percentile / prop_VIP_in_percentile_98$n_with_iHS

###############################################################################
#### Inspecting outliers: filter for 99th percentile and above iHS scores #####
###############################################################################
percentile_99 <- percentiles %>% filter(within_pop_percentile >= 99)

ggplot(percentile_99, aes(x = pop, y = ihs, fill = type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_minimal() +
  labs(title = "Boxplot of iHS Scores above 99th percentile by Population for VIPs vs. non-VIPs",
       x = "Population",
       y = "iHS",
       fill = "Type") +
  theme(plot.title = element_text(hjust = 0.5))

# raw counts of VIPs vs non VIPs in the 99th percentile
prop_VIP_in_percentile_99 <- percentile_99 %>% group_by(pop, type) %>% count()
colnames(prop_VIP_in_percentile_99) <- c("pop", "type", "n_in_99th_percentile")

# add column of total genes in the 99th percentile
prop_VIP_in_percentile_99 <- prop_VIP_in_percentile_99 %>%
  group_by(pop) %>%
  mutate(total_n_in_99th_percentile = sum(n_in_99th_percentile, na.rm = TRUE))

# calculate prop of VIPs vs non VIPS in the 99th percentile
prop_VIP_in_percentile_99$prop_in_99th_percentile <- prop_VIP_in_percentile_99$n_in_99th_percentile / prop_VIP_in_percentile_99$total_n_in_99th_percentile

# add in counts of VIPs vs non VIPS that have iHS values
prop_VIP_in_percentile_99$n_with_iHS <- n_genes_with_iHS$n

# calculate prop of VIPs in the 99th percentile out of all VIPs with iHS values
prop_VIP_in_percentile_99$prop_in_99th_percentile_among_type_with_iHS <- prop_VIP_in_percentile_99$n_in_99th_percentile / prop_VIP_in_percentile_99$n_with_iHS

#################################################################################
#### Inspecting outliers: filter for 99.5th percentile and above iHS scores #####
#################################################################################

percentile_99.5 <- percentiles %>% filter(within_pop_percentile >= 99.5)

ggplot(percentile_99.5, aes(x = pop, y = ihs, fill = type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_minimal() +
  labs(title = "Boxplot of iHS Scores above 99.5th percentile by Population for VIPs vs. non-VIPs",
       x = "Population",
       y = "iHS",
       fill = "Type") +
  theme(plot.title = element_text(hjust = 0.5))

# raw counts of VIPs vs non VIPs in the 99.5th percentile
prop_VIP_in_percentile_99.5 <- percentile_99.5 %>% group_by(pop, type) %>% count()
colnames(prop_VIP_in_percentile_99.5) <- c("pop", "type", "n_in_99.5th_percentile")

# add column of total genes in the 99.5th percentile
prop_VIP_in_percentile_99.5 <- prop_VIP_in_percentile_99.5 %>%
  group_by(pop) %>%
  mutate(total_n_in_99.5th_percentile = sum(n_in_99.5th_percentile, na.rm = TRUE))

# calculate prop of VIPs vs non VIPS in the 99.5th percentile
prop_VIP_in_percentile_99.5$prop_in_99.5th_percentile <- prop_VIP_in_percentile_99.5$n_in_99.5th_percentile / prop_VIP_in_percentile_99.5$total_n_in_99.5th_percentile

# add in counts of VIPs vs non VIPS that have iHS values
prop_VIP_in_percentile_99.5$n_with_iHS <- n_genes_with_iHS$n

# calculate prop of VIPs in the 99.5th percentile out of all VIPs with iHS values
prop_VIP_in_percentile_99.5$prop_in_99.5th_percentile_among_type_with_iHS <- prop_VIP_in_percentile_99.5$n_in_99.5th_percentile / prop_VIP_in_percentile_99.5$n_with_iHS

###############################################################################
#### Inspecting outliers: filter for 50th percentile and above iHS scores #####
###############################################################################
percentile_50 <- percentiles %>% filter(within_pop_percentile >= 50)

ggplot(percentile_50, aes(x = pop, y = ihs, fill = type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_minimal() +
  labs(title = "Boxplot of iHS Scores above 50th percentile by Population for VIPs vs. non-VIPs",
       x = "Population",
       y = "iHS",
       fill = "Type") +
  theme(plot.title = element_text(hjust = 0.5))

# raw counts of VIPs vs non VIPs in the 50th percentile
prop_VIP_in_percentile_50 <- percentile_50 %>% group_by(pop, type) %>% count()
colnames(prop_VIP_in_percentile_50) <- c("pop", "type", "n_in_50th_percentile")

# add column of total genes in the 50th percentile
prop_VIP_in_percentile_50 <- prop_VIP_in_percentile_50 %>%
  group_by(pop) %>%
  mutate(total_n_in_50th_percentile = sum(n_in_50th_percentile, na.rm = TRUE))

# calculate prop of VIPs vs non VIPS in the 50th percentile
prop_VIP_in_percentile_50$prop_in_50th_percentile <- prop_VIP_in_percentile_50$n_in_50th_percentile / prop_VIP_in_percentile_50$total_n_in_50th_percentile

# add in counts of VIPs vs non VIPS that have iHS values
prop_VIP_in_percentile_50$n_with_iHS <- n_genes_with_iHS$n

# calculate prop of VIPs in the 50th percentile out of all VIPs with iHS values
prop_VIP_in_percentile_50$prop_in_50th_percentile_among_type_with_iHS <- prop_VIP_in_percentile_50$n_in_50th_percentile / prop_VIP_in_percentile_50$n_with_iHS

#######################################################################
####  Construct table of outliers at diff. percentile thresholds  #####
#######################################################################
prop_VIP_by_percentiles <- cbind(prop_VIP_in_percentile_50$prop_in_50th_percentile, 
                                 prop_VIP_in_percentile_50$prop_in_50th_percentile_among_type_with_iHS,
                                 prop_VIP_in_percentile_90$prop_in_90th_percentile, 
                                 prop_VIP_in_percentile_90$prop_in_90th_percentile_among_type_with_iHS,
                                 prop_VIP_in_percentile_95, 
                                 prop_VIP_in_percentile_97$prop_in_97th_percentile, 
                                 prop_VIP_in_percentile_97$prop_in_97th_percentile_among_type_with_iHS, 
                                 prop_VIP_in_percentile_98$prop_in_98th_percentile, 
                                 prop_VIP_in_percentile_98$prop_in_98th_percentile_among_type_with_iHS,
                                 prop_VIP_in_percentile_99$prop_in_99th_percentile, 
                                 prop_VIP_in_percentile_99$prop_in_99th_percentile_among_type_with_iHS,
                                 prop_VIP_in_percentile_99.5$prop_in_99.5th_percentile, 
                                 prop_VIP_in_percentile_99.5$prop_in_99.5th_percentile_among_type_with_iHS)

prop_VIP_by_percentiles <- prop_VIP_by_percentiles[,c(5,6,1:4,9,11,16:23)]
colnames(prop_VIP_by_percentiles) <- gsub(".*\\$", "", colnames(prop_VIP_by_percentiles))

# calculate enrichment ratios at various percentiles
prop_VIP_enrichment <- prop_VIP_by_percentiles %>%
  group_by(pop) %>%
  summarize(VIP_enrichment_50th = prop_in_50th_percentile_among_type_with_iHS[type == "VIP"] / prop_in_50th_percentile_among_type_with_iHS[type == "non-VIP"],
            VIP_enrichment_90th = prop_in_90th_percentile_among_type_with_iHS[type == "VIP"] / prop_in_90th_percentile_among_type_with_iHS[type == "non-VIP"],
            VIP_enrichment_95th = prop_in_95th_percentile_among_type_with_iHS[type == "VIP"] / prop_in_95th_percentile_among_type_with_iHS[type == "non-VIP"],
            VIP_enrichment_97th = prop_in_97th_percentile_among_type_with_iHS[type == "VIP"] / prop_in_97th_percentile_among_type_with_iHS[type == "non-VIP"],
            VIP_enrichment_98th = prop_in_98th_percentile_among_type_with_iHS[type == "VIP"] / prop_in_98th_percentile_among_type_with_iHS[type == "non-VIP"],
            VIP_enrichment_99th = prop_in_99th_percentile_among_type_with_iHS[type == "VIP"] / prop_in_99th_percentile_among_type_with_iHS[type == "non-VIP"],
            VIP_enrichment_99.5th = prop_in_99.5th_percentile_among_type_with_iHS[type == "VIP"] / prop_in_99.5th_percentile_among_type_with_iHS[type == "non-VIP"]) %>%
  ungroup()

# plot by population
prop_VIP_enrichment_long <- prop_VIP_enrichment %>%
  pivot_longer(
    cols = starts_with("VIP_enrichment"), # Select columns to gather
    names_to = "Percentile", # New column for the 'VIP_enrichment' labels
    values_to = "Value" # New column for the values
  )

prop_VIP_enrichment_long <- prop_VIP_enrichment_long %>%
  mutate(Percentile = str_replace(Percentile, "VIP_enrichment_", ""), # Remove 'VIP_enrichment_'
         Percentile = str_replace(Percentile, "th", "")) # Remove 'th'

prop_VIP_enrichment_long$Percentile <- as.numeric(prop_VIP_enrichment_long$Percentile)
prop_VIP_enrichment_long$Percentile <- factor(prop_VIP_enrichment_long$Percentile,
                                              levels = unique(prop_VIP_enrichment_long$Percentile),
                                              ordered = TRUE)

# add back metadata
prop_VIP_enrichment_long <- merge(prop_VIP_enrichment_long, metadata, by.x = "pop", by.y = "pop_code", all.x = "TRUE")

# order mapping of populations by region
pop_region_mapping <- prop_VIP_enrichment_long %>%
  distinct(pop, region) %>%
  mutate(region = factor(region, levels = c("South Asia","Southeast Asia", "East Asia", "Oceania", 
                                            "West Central Africa", "East Africa", 
                                            "Southern Africa", "East Central Africa"))) %>%
  arrange(region)

prop_VIP_enrichment_long$pop <- factor(prop_VIP_enrichment_long$pop, levels = pop_region_mapping$pop)

ggplot(prop_VIP_enrichment_long, aes(x = Percentile, y = Value, group = interaction(pop, subsistence), color = region, shape = subsistence)) +
  geom_line() +
  geom_point(na.rm = TRUE) + # Consider adding na.rm = TRUE to remove NA warning
  facet_wrap(~pop, scales = "free_y") +
  scale_x_discrete(breaks = unique(prop_VIP_enrichment_long$Percentile)) +
  theme_minimal() +
  labs(title = "VIP Enrichment by Population",
       x = "Percentile",
       y = "Enrichment Ratio") +
  scale_color_manual(values = c("South Asia" = "red", "Southeast Asia" = "deeppink", "East Asia" = "plum2", 
                                "Oceania" = "darkolivegreen2", "West Central Africa" = "darkslategray2",
                                "East Africa" = "blue4", "Southern Africa" = "cornflowerblue", 
                                "East Central Africa" = "slateblue3")) +
  scale_shape_manual(values = c("Hunter-gatherer" = 17, "Agriculturalist" = 1, "Pastoralist" = 3)) +
  scale_y_continuous(limits = c(min(prop_VIP_enrichment_long$Value), max(prop_VIP_enrichment_long$Value))) # min and max values
```

## Results

---

**Metabolic enrichments stratified by the KEGG subcategories (e.g. carbohydrate vs. amino acid vs. lipid metabolism enrichment):**

[tajima_metabolic_subcategory_enrich (1).pdf](tajima_metabolic_subcategory_enrich_(1).pdf)

[ihs_metabolic_subcategory_enrich (1).pdf](ihs_metabolic_subcategory_enrich_(1).pdf)

**Overall metabolic enrichment and VIP enrichment plots generated by Evan:**

[tajima_vip_metabolic_enrich.pdf](tajima_vip_metabolic_enrich%201.pdf)

[ihs_vip_metabolic_enrich (1).pdf](ihs_vip_metabolic_enrich_(1).pdf)