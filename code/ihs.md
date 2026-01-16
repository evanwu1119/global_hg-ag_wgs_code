# iHS

- Working directory: `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet_SHAPEIT5phased_ancstate_nokin/per_pop`
- Jump to the bottom of this page to navigate to the Akbari et al., and cCRE analyses

## Run iHS

---

**Below you're going to call every population in the poplist_filtered.txt file line by line**

**In terminal:**

```bash
cd /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet_SHAPEIT5phased_ancstate_nokin/per_pop

for file in *.impute.hap; do
    pop=${file%.impute.hap}
    export pop
    sbatch hapbin_ihs.sbatch
done
```

- **`for file in *.impute.hap; do ... done`**: This loops over all **`.impute.hap`** files.
- **`pop=${file%.impute.hap}`**: This extracts the prefix from the filename (e.g., **`AET`** from **`AET.impute.hap`**).
- **`export pop`**: This makes the **`pop`** variable available to the Slurm script.
- **`sbatch hapbin_ihs.sbatch`**: This submits the Slurm job.

**In `hapbin_ihs.sbatch`:**

```bash
#!/bin/bash
#SBATCH --partition=caslake
#SBATCH --account=pi-\lbarreiro
#SBATCH --time=10:00:00
#SBATCH --mem=90G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6

ihsbin --hap ${pop}.impute.hap --map ${pop}.map --out ${pop}_iHS --minmaf 0.05 --cutoff 0.1

# iHS calculated for all variants with a minor allele frequency greater than 5% and up to the point at which EHH drops below 0.1
```

## iHS analysis - BC pipeline

---

Step 1) Generate iHS stats per gene using 100kb windows centered on each gene with `freeze2_iHS_per_gene_all_pops_BCpipeline.R` and `freeze2_iHS_per_gene_all_pops_BCpipeline.sbatch`

```bash
#!/bin/bash
#SBATCH --partition=caslake
#SBATCH --account=pi-lbarreiro
#SBATCH --time=36:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=iHS_per_gene_all_pops

DIR=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS

cd ${DIR}

module load R/4.1.0

Rscript freeze2_iHS_per_gene_all_pops_BCpipeline.R
```

Resulting output is `freeze2_iHS_per_gene_all_pops_BCpipeline.RData`

- Within this .RData file, the `gff3_ihs` object has the per-population iHS results for each gene present in the gff3 gene list

## Exploratory iHS Analysis (including VIP enrichment) - BC pipeline

---

Use `freeze2_exploratory_analysis_iHS_per_gene_all_pops_BCpipeline.R`

Image is saved in `freeze2_exploratory_analysis_iHS_per_gene_all_pops_BCpipeline.RData`

![image.png](image.png)

- **Mean iHS scores by population for VIPs vs. non-VIPs (expand to view)**
    
    
    | **pop** | **type** | **mean_iHS** |
    | --- | --- | --- |
    | ADV | VIP | 0.8979815 |
    | ADV | non-VIP | 0.8300674 |
    | AET | VIP | 0.9020806 |
    | AET | non-VIP | 0.8302765 |
    | AGT | VIP | 0.8901982 |
    | AGT | non-VIP | 0.8405501 |
    | ATI | VIP | 0.8966457 |
    | ATI | non-VIP | 0.8375965 |
    | ATY | VIP | 0.8255931 |
    | ATY | non-VIP | 0.8017765 |
    | BAK | VIP | 0.8561564 |
    | BAK | non-VIP | 0.8099761 |
    | BCH | VIP | 0.8197378 |
    | BCH | non-VIP | 0.7836707 |
    | CBU | VIP | 0.8645471 |
    | CBU | non-VIP | 0.8069544 |
    | FAN | VIP | 0.8384088 |
    | FAN | non-VIP | 0.7948631 |
    | GAL | VIP | 0.8469315 |
    | GAL | non-VIP | 0.8028093 |
    | HMB | VIP | 0.8577092 |
    | HMB | non-VIP | 0.8039231 |
    | KIG | VIP | 0.8341722 |
    | KIG | non-VIP | 0.7925454 |
    | KOT | VIP | 0.8394823 |
    | KOT | non-VIP | 0.7972493 |
    | KWI | VIP | 0.8693245 |
    | KWI | non-VIP | 0.8256120 |
    | MNO | VIP | 0.8848648 |
    | MNO | non-VIP | 0.8191030 |
    | MTW | VIP | 0.8940705 |
    | MTW | non-VIP | 0.8329736 |
    | MWA | VIP | 0.8843143 |
    | MWA | non-VIP | 0.8155954 |
    | NAM | VIP | 0.8232502 |
    | NAM | non-VIP | 0.7777201 |
    | NZE | VIP | 0.8463063 |
    | NZE | non-VIP | 0.8045748 |
    | PAI | VIP | 0.8000008 |
    | PAI | non-VIP | 0.7860549 |
    | PHP | VIP | 0.8797662 |
    | PHP | non-VIP | 0.8172447 |
    | PNG | VIP | 0.8301695 |
    | PNG | non-VIP | 0.7895384 |
    | POL | VIP | 0.8679895 |
    | POL | non-VIP | 0.8222525 |
    | RAM | VIP | 0.8988901 |
    | RAM | non-VIP | 0.8356672 |
    | SAN | VIP | 0.8012317 |
    | SAN | non-VIP | 0.7644359 |
    | SCI | VIP | 0.8858259 |
    | SCI | non-VIP | 0.8270895 |
    | SHK | VIP | 0.8318676 |
    | SHK | non-VIP | 0.7902012 |
    | SHL | VIP | 0.8839827 |
    | SHL | non-VIP | 0.8180139 |
    | SLI | VIP | 0.8848965 |
    | SLI | non-VIP | 0.8152154 |
    | SMB | VIP | 0.8856888 |
    | SMB | non-VIP | 0.8302133 |
    | TEK | VIP | 0.8535571 |
    | TEK | non-VIP | 0.8029487 |
    | TSO | VIP | 0.8246203 |
    | TSO | non-VIP | 0.8035867 |
    | TWA | VIP | 0.8526046 |
    | TWA | non-VIP | 0.7986489 |
    | VAN | VIP | 0.8282964 |
    | VAN | non-VIP | 0.7947710 |

![image.png](image%201.png)

![image.png](image%202.png)

![image.png](image%203.png)

![image.png](image%204.png)

## GSEA analysis w/GO terms - BC pipeline

---

- Refer to `freeze2_iHS_GSEA_GO_all_pops_BCpipeline.R`  and                     
`freeze2_iHS_GSEA_GO_all_pops_BCpipeline.RData`

Note: for my GSEA pipeline, I’ve only used Gene Ontology terms. This is different than Evan’s annotation steps, which use GO, KEGG, and MSigDB terms. The per-gene statistics are also based on a 100kb window centered on each gene.

**Heatmap of significant terms, subsetted GO terms to second children**: 

[ORA.pdf](ORA.pdf)

[GSEA.pdf](GSEA.pdf)

## GSEA and ORA analysis w/GO, KEGG, and MSigDB terms - EW pipeline

---

- Refer to `freeze2_iHS_GSEA_GO_all_pops_v2_EWannot.R`
    - Enrichment and over-representation analysis results are stored in `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/freeze2_iHS_GSEA_GO_all_pops_v2_EWannot_outputs`
- Heatmaps are generated with `freeze2_iHS_GSEA_GO_all_pops_v2_EWannot_plots.R`
    
    [GSEA.pdf](GSEA%201.pdf)
    
    [ORA.pdf](ORA%201.pdf)
    

**Re-doing with all variants within a gene:**

- Refer to `freeze2_iHS_GSEA_GO_all_pops_v3_allvariants.R`
    - Enrichment and over-representation analysis results are stored in `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/freeze2_iHS_GSEA_GO_all_pops_v3_allvariants_outputs`

## iHS per 100kb window

---

- Refer to `freeze2_iHS_per_window_100kb.R`
    - Per-window statistics for each population are saved as `freeze2_iHS_per_window_100kb_df.RData`
        - use `unique_pops.RData` for the column names (populations)
        - use `GRCh37_100kb_intervals.RData` for the intervals
    - Per-window SNP count for each population is saved as
    `freeze2_iHS_per_window_100kb_count_snps_df.RData`
    - Per-population manhattan plots are saved in `/freeze2_iHS_per_window_100kb_manhattan_plots`

Example distribution of SNP counts:

![image.png](image%205.png)

[nsnps_per_interval_per_population.pdf](nsnps_per_interval_per_population.pdf)

**Note: x-axis limit is set to 500 but there are windows with > 1000 SNPs. For all analyses I’ve filtered out windows with 5 or less SNPs**

Example per-population manhattan plot:

![image.png](image%206.png)

Example of cross-pop comparisons at 99th percentile:

![image.png](image%207.png)

![image.png](image%208.png)

## **Regional HG + AGR comparisons**

---

![image.png](image%209.png)

![image.png](image%2010.png)

![image.png](image%2011.png)

![image.png](image%2012.png)

# Evan’s additional iHS analyses

---

Working directory: `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/fst/ihs` 

I used the data provided by Bridget for these analyses (`/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/freeze2_iHS_per_gene_all_pops_BCpipeline.RData`), as well as my own Tajima’s statistics (`tajimas.rds`) where applicable. All code is contained in `ihs_shared_enrich.Rmd`

**Step 1) Test the enrichment of metabolic genesets for extreme iHS values in each population**

[ihs_vip_metabolic_enrich.pdf](ihs_vip_metabolic_enrich.pdf)

[ihs_metabolic_subcategory_enrich.pdf](ihs_metabolic_subcategory_enrich.pdf)

```r
library(clusterProfiler)
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(plyranges)
library(rtracklayer)
library(tidyverse)
library(biomaRt)
library(ggridges)
library(ggragged)

load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/new_iHS_per_gene_all_pops_lowconf_anc_annot_BCpipeline_gff3_df.RData")
metadata <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
metadata$region <- factor(metadata$region, levels = c("West Central Africa", "East Central Africa", "East Africa", "Southern Africa", "South Asia", "East Asia", "Southeast Asia", "Oceania"))
metadata$subsistence <- factor(metadata$subsistence, levels = c("Hunter-gatherer", "Agriculturalist", "Pastoralist"))
metadata$nindiv <- c(12,28,17,18,10,23,37,24,30,24,31,25,30,11,23,26,16,40,28,7,20,8,17,18,57,14,34,13,34,28,25,8,41,141)
metadata$trust <- ifelse(metadata$nindiv < 15, "no", "yes")
VIP_list <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze1/metadata/from_providers/viral_effects_file_april2023.txt", header = FALSE)
# metabolic_genes <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/curate_genes/KEGG_09100_Metabolism/09100_METABOLISM_ALLGENES.txt", header = F)
# metabolic_genes <- sapply(unique(metabolic_genes$V2),  function(x) strsplit(x, ":")[[1]][2])
# ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# metabolic_ensg <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id'),
#       filters = 'entrezgene_id',
#       values = metabolic_genes,
#       mart = ensembl)
# saveRDS(metabolic_ensg, paste0(wd, "data/metabolic_ensg.rds")) # Doesn't map 1:1, 1679 unique Entrez/KEGG IDs vs. 1858 unique Ensembl IDs
metabolic_ensg <- readRDS(paste0(wd, "data/metabolic_ensg.rds"))

# Format per-population iHS results and calculate percentile per-gene
# gff3_results <- cbind(gff3, gff3_ihs)
# names(gff3_results)[14:47] <- metadata$pop_code
gff3_results <- gff3_results %>% gather(pop, ihs, 14:47)
gff3_results <- merge(gff3_results, metadata, by.x = "pop", by.y = "pop_code", all.x = TRUE)
gff3_results$type <- ifelse(gff3_results$gene_IDs %in% VIP_list$V1, "VIP", "non-VIP")
gff3_results$metabolic <- ifelse(gff3_results$gene_IDs %in% metabolic_ensg$ensembl_gene_id, "metabolic", "non-metabolic")
gff3_results <- gff3_results %>%
  group_by(pop) %>%
  mutate(within_pop_percentile = if_else(is.na(ihs),
                                         as.numeric(NA),
                                         rank(ihs, ties.method = "average") / sum(!is.na(ihs)) * 100))
saveRDS(gff3_results, paste0(wd, "data/ihs_gff3_result.rds"))
gff3_results <- readRDS(paste0(wd, "data/ihs_gff3_result.rds"))

# Format percentiles to plot proportion of genes that in a category (VIP, metabolic) at certain cutoffs, don't calculate if less than 5 genes at a cutoff
prop_pct <- function(percentiles, test, category, cutoff) {
  # Get the smaller tail of dist
  if(cutoff < 50) {
    tmp <- percentiles %>% filter(within_pop_percentile <= cutoff) %>% group_by(pop, get(category)) %>% count()
  } else {
    tmp <- percentiles %>% filter(within_pop_percentile >= cutoff) %>% group_by(pop, get(category)) %>% count()
  }
  colnames(tmp)[3] <- "n_category"
  tmp <- percentiles %>% filter(!is.na(get(test))) %>% group_by(pop, get(category)) %>% count() %>% left_join(tmp, by = c("pop", "get(category)"))
  colnames(tmp) <- c("pop", category, "n_total", "n_category")
  tmp$n_category[tmp$n_category < 5] <- NA
  tmp[[as.character(cutoff)]] <- tmp$n_category / tmp$n_total
  return(tmp[,c(1,2,5)])
  # colnames(tmp)[4] <- c(paste(cutoff, "n_category"))
  # return(tmp)
}

# Convert percentiles dataframe to enrichment ratio
pct_to_er <- function(percentiles, category, enrich_str) {
  tmp_list <- lapply(unique(percentiles$pop), function(x) data.frame(
    pop = rep(x, length(unique(percentiles$percentile))),
    percentiles = unique(percentiles$percentile),
    enrich_ratio = percentiles$proportion[percentiles$pop == x & percentiles[[`category`]] == enrich_str] / percentiles$proportion[percentiles$pop == x & percentiles[[`category`]] != enrich_str]))
  tmp <- bind_rows(tmp_list)
  tmp <- tmp %>% left_join(unique(percentiles[,c(1,3,4,5)]), by = "pop")
  return(tmp)
}

cutoffs <- c(50, 90, 95, 97, 98, 99)

prop_vip_list <- lapply(cutoffs, prop_pct, percentiles = gff3_results, test = "ihs", category = "type")
prop_vip_list <- purrr::reduce(prop_vip_list, left_join, by = c("pop", "type")) %>% left_join(metadata[,c(1,3,5,7)], by = c("pop" = "pop_code")) %>% pivot_longer(cols = as.character(cutoffs), names_to = "percentile", values_to = "proportion")
prop_vip <- pct_to_er(prop_vip_list, "type", "VIP")

prop_metabolic_list <- lapply(cutoffs, prop_pct, percentiles = gff3_results, test = "ihs", category = "metabolic")
prop_metabolic_list <- purrr::reduce(prop_metabolic_list, left_join, by = c("pop", "metabolic")) %>% left_join(metadata[,c(1,3,5,7)], by = c("pop" = "pop_code")) %>% pivot_longer(cols = as.character(cutoffs), names_to = "percentile", values_to = "proportion")
prop_metabolic <- pct_to_er(prop_metabolic_list, "metabolic", "metabolic")

pdf(paste0(wd, "figures/ihs_vip_metabolic_enrich.pdf"), height = 10, width = 10)
# Plot enriched proportions
prop_vip %>% ggplot(aes(x = percentiles, y = enrich_ratio, color = subsistence, group = 1)) +
  geom_hline(yintercept = 1, color = "black", linetype = 2) +
  geom_line() +
  geom_point(na.rm = TRUE) +
  scale_color_manual(values = c(
    "Hunter-gatherer" = "#FFB6C1",
    "Agriculturalist" = "#FDAE61",
    "Pastoralist" = "#40E0D0")) +
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
  scale_color_manual(values = c(
    "Hunter-gatherer" = "#FFB6C1",
    "Agriculturalist" = "#FDAE61",
    "Pastoralist" = "#40E0D0")) +
  facet_ragged_rows(vars(region), vars(pop), axes = "all_x") +
  theme_classic() +
  labs(title = "Metabolic genes iHS enrichment", x = "Percentile", y = "Enrichment Ratio", color = "Subsistence") +
  theme(strip.text.y.right = element_text(angle = 0),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
dev.off()

# Get genes from 11 subcategories of all metabolic genes
# metabolic_sub <- list.files(path = "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/curate_genes/KEGG_09100_Metabolism/", pattern = "ALLGENES", full.names = T) %>% grep(pattern = "09100_METABOLISM", value = T, invert = T)
# metabolic_sub_list <- lapply(metabolic_sub, function(x) sapply(unique(fread(x, header = F)$V2), function(y) strsplit(y, ":")[[1]][2]))
# names(metabolic_sub_list) <- str_sub(basename(metabolic_sub), 7, -14)
# ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# metabolic_sub_ensg <- lapply(metabolic_sub_list, function(x) getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id'), filters = 'entrezgene_id', values = x, mart = ensembl))
# saveRDS(metabolic_sub_ensg, paste0(wd, "data/metabolic_sub_ensg.rds"))
metabolic_sub_ensg <- readRDS(paste0(wd, "data/metabolic_sub_ensg.rds"))

# iHS enrichments
cutoffs <- c(50, 80, 90, 93, 94, 95)
ihs_metabolic_sub <- list()

# Calculate raw numbers
# for (m in names(metabolic_sub_ensg)) {
#   gff3_results[[`m`]] <-ifelse(gff3_results$gene_IDs %in% metabolic_sub_ensg[[m]]$ensembl_gene_id, "metabolic", "non-metabolic")
#
#   ihs_metabolic_list <- lapply(cutoffs, prop_pct, percentiles = gff3_results, test = "ihs", category = m)
#   ihs_metabolic_sub[[`m`]] <- purrr::reduce(ihs_metabolic_list, left_join, by = c("pop", m, "n_total"))
# }
# ihs_metabolic_sub_summary <- lapply(ihs_metabolic_sub, function(x) summary(x[x[,2] == "metabolic", c(4,6,8,10,12)]))

for (m in names(metabolic_sub_ensg)) {
  gff3_results[[`m`]] <-ifelse(gff3_results$gene_IDs %in% metabolic_sub_ensg[[m]]$ensembl_gene_id, "metabolic", "non-metabolic")

  ihs_metabolic_list <- lapply(cutoffs, prop_pct, percentiles = gff3_results, test = "ihs", category = m)
  ihs_metabolic_list <- purrr::reduce(ihs_metabolic_list, left_join, by = c("pop", m)) %>% left_join(metadata[,c(1,3,5,7)], by = c("pop" = "pop_code")) %>% pivot_longer(cols = as.character(cutoffs), names_to = "percentile", values_to = "proportion")
  ihs_metabolic_sub[[`m`]] <- pct_to_er(ihs_metabolic_list, m, "metabolic")
}

pdf(paste0(wd, "figures/ihs_metabolic_subcategory_enrich.pdf"), height = 10, width = 10)
for (m in names(ihs_metabolic_sub)) {
  # Plot enriched proportions
  p <- ihs_metabolic_sub[[m]] %>% ggplot(aes(x = percentiles, y = enrich_ratio, color = subsistence, group = 1)) +
    #geom_rect(aes(fill = trust), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, show.legend = F, inherit.aes = F) +
    geom_hline(yintercept = 1, color = "black", linetype = 2) +
    geom_line() +
    geom_point(na.rm = TRUE) +
    #scale_fill_manual(values = c("grey", "white")) +
    scale_color_manual(values = c(
      "Hunter-gatherer" = "#FFB6C1",
      "Agriculturalist" = "#FDAE61",
      "Pastoralist" = "#40E0D0")) +
    facet_ragged_rows(vars(region), vars(pop), axes = "all_x") +
    theme_classic() +
    labs(title = paste(gsub("_", " ", m), "genes iHS enrichment"), x = "Percentile", y = "Enrichment Ratio", color = "Subsistence") +
    theme(strip.text.y.right = element_text(angle = 0),
          strip.background = element_blank(),
          #panel.background = element_rect(color = back_color),
          axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 5))
  print(p)
}
dev.off()

```

**Step 2) Do the same but two-sided test for Tajima’s since we care about both positive and negative outliers**

[tajima_vip_metabolic_enrich.pdf](tajima_vip_metabolic_enrich.pdf)

[tajima_metabolic_subcategory_enrich.pdf](tajima_metabolic_subcategory_enrich.pdf)

```r
# Do for Tajima's
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
# saveRDS(tajimas, paste0(wd, "data/tajimas.rds"))
tajimas <- readRDS(paste0(wd, "data/tajimas.rds"))

cutoffs <- c(1, 3, 5, 10, 50, 90, 95, 97, 99)
tajima_vip_list <- lapply(cutoffs, prop_pct, percentiles = tajimas, test = "tajima", category = "type")
tajima_vip_list <- purrr::reduce(tajima_vip_list, left_join, by = c("pop", "type")) %>% left_join(metadata[,c(1,3,5,7)], by = c("pop" = "pop_code")) %>% pivot_longer(cols = as.character(cutoffs), names_to = "percentile", values_to = "proportion")
tajima_vip <- pct_to_er(tajima_vip_list, "type", "VIP")

tajima_metabolic_list <- lapply(cutoffs, prop_pct, percentiles = tajimas, test = "tajima", category = "metabolic")
tajima_metabolic_list <- purrr::reduce(tajima_metabolic_list, left_join, by = c("pop", "metabolic")) %>% left_join(metadata[,c(1,3,5,7)], by = c("pop" = "pop_code")) %>% pivot_longer(cols = as.character(cutoffs), names_to = "percentile", values_to = "proportion")
tajima_metabolic <- pct_to_er(tajima_metabolic_list, "metabolic", "metabolic")

pdf(paste0(wd, "figures/tajima_vip_metabolic_enrich.pdf"), height = 10, width = 10)
# Plot enriched proportions
tajima_vip %>% ggplot(aes(x = percentiles, y = log2(enrich_ratio), color = subsistence, group = 1)) +
  geom_hline(yintercept = 0, color = "black", linetype = 2) +
  geom_line() +
  geom_point(na.rm = TRUE) +
  scale_color_manual(values = c(
    "Hunter-gatherer" = "#FFB6C1",
    "Agriculturalist" = "#FDAE61",
    "Pastoralist" = "#40E0D0")) +
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
  scale_color_manual(values = c(
    "Hunter-gatherer" = "#FFB6C1",
    "Agriculturalist" = "#FDAE61",
    "Pastoralist" = "#40E0D0")) +
  facet_ragged_rows(vars(region), vars(pop), axes = "all_x") +
  theme_classic() +
  labs(title = "Metabolic genes Tajima's enrichment", x = "Percentile", y = expression(paste("log"[2], " Enrichment Ratio")), color = "Subsistence") +
  theme(strip.text.y.right = element_text(angle = 0),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
dev.off()

# Do metabolic subcategories for Tajima's
cutoffs <- c(5, 6, 7, 10, 20, 50, 80, 90, 93, 94, 95)
tajima_metabolic_sub <- list()

# Calculate raw numbers
# for (m in names(metabolic_sub_ensg)) {
#   tajimas[[`m`]] <-ifelse(tajimas$gene_IDs %in% metabolic_sub_ensg[[m]]$ensembl_gene_id, "metabolic", "non-metabolic")
#   tajima_metabolic_list <- lapply(cutoffs, prop_pct, percentiles = tajimas, test = "tajima", category = m)
#   tajima_metabolic_sub[[`m`]] <- purrr::reduce(tajima_metabolic_list, left_join, by = c("pop", m, "n_total"))
# }
# tajima_metabolic_sub_summary <- lapply(tajima_metabolic_sub, function(x) summary(x[x[,2] == "metabolic", c(4,6,8,10,12,14,16)]))
# save(ihs_metabolic_sub, ihs_metabolic_sub_summary, tajima_metabolic_sub, tajima_metabolic_sub_summary, file = paste0(wd, "data/metabolic_subcategory_raw_gene_cutoffs.RData"))

for (m in names(metabolic_sub_ensg)) {
  tajimas[[`m`]] <-ifelse(tajimas$gene_IDs %in% metabolic_sub_ensg[[m]]$ensembl_gene_id, "metabolic", "non-metabolic")

  tajima_metabolic_list <- lapply(cutoffs, prop_pct, percentiles = tajimas, test = "tajima", category = m)
  tajima_metabolic_list <- purrr::reduce(tajima_metabolic_list, left_join, by = c("pop", m)) %>% left_join(metadata[,c(1,3,5,7)], by = c("pop" = "pop_code")) %>% pivot_longer(cols = as.character(cutoffs), names_to = "percentile", values_to = "proportion")
  tajima_metabolic_sub[[`m`]] <- pct_to_er(tajima_metabolic_list, m, "metabolic")
  tajima_metabolic_sub[[`m`]]$percentiles <- factor(tajima_metabolic_sub[[`m`]]$percentiles, levels = cutoffs)
}

pdf(paste0(wd, "figures/tajima_metabolic_subcategory_enrich.pdf"), height = 10, width = 10)
for (m in names(tajima_metabolic_sub)) {
  # Plot enriched proportions
  p <- tajima_metabolic_sub[[m]] %>% ggplot(aes(x = percentiles, y = log2(enrich_ratio), color = subsistence, group = 1)) +
    #geom_rect(aes(fill = trust), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, show.legend = F, inherit.aes = F) +
    geom_hline(yintercept = 0, color = "black", linetype = 2) +
    geom_line() +
    geom_point(na.rm = TRUE) +
    #scale_fill_manual(values = c("grey", "white")) +
    scale_color_manual(values = c(
      "Hunter-gatherer" = "#FFB6C1",
      "Agriculturalist" = "#FDAE61",
      "Pastoralist" = "#40E0D0")) +
    facet_ragged_rows(vars(region), vars(pop), axes = "all_x") +
    theme_classic() +
    labs(title = paste(gsub("_", " ", m), "genes Tajima's enrichment"), x = "Percentile", y = expression(paste("log"[2], " Enrichment Ratio")), color = "Subsistence") +
    theme(strip.text.y.right = element_text(angle = 0),
          strip.background = element_blank(),
          axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 5))
  print(p)
}
dev.off()
```

**Step 3) Retrieved the 12,111 SNPs associated with height from the latest multi-ancestry GWAS meta-analysis of height by Yengo et al. (**[https://www.nature.com/articles/s41586-022-05275-y](https://www.nature.com/articles/s41586-022-05275-y)**) and plotted iHS values of these SNPs to see if they differ from the global distribution**

Calculating the iHS summary statistics of all SNPs `ihs_global_summary.R`

```r
library(parallel)
library(data.table)

cores <- 7
ihs_files <- list.files("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/outputs", full.names = T)

# Calculate per population genome-wide iHS summaries
names <- substr(basename(ihs_files), 1, 3)
ihs_outputs <- mclapply(ihs_files, function(x) fread(x, select = 7, header = T), mc.cores = cores)
ihs_summaries <- mclapply(ihs_outputs, function(x) summary(na.omit(x$`Std iHS`)), mc.cores = cores)
names(ihs_summaries) <- names
abs_summaries <- mclapply(ihs_outputs, function(x) summary(na.omit(abs(x$`Std iHS`))), mc.cores = cores)
names(abs_summaries) <- names

# Summary across all populations
ihs_summaries[["global"]] <- summary(na.omit(unlist(ihs_outputs, use.names = F)))
abs_summaries[["global"]] <- summary(na.omit(abs(unlist(ihs_outputs, use.names = F))))

save(ihs_summaries, abs_summaries, file = "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/fst/ihs/ihs_summaries.RData")

```

[yengo_snp_ihs_region.pdf](yengo_snp_ihs_region.pdf)

```r
popinfo <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
popinfo$region <- factor(popinfo$region, levels = c(levels = c("West Central Africa", "East Central Africa", "East Africa", "Southern Africa", "South Asia", "East Asia", "Southeast Asia", "Oceania")))
popinfo$subsistence <- factor(popinfo$subsistence, levels = c("Hunter-gatherer", "Agriculturalist", "Pastoralist"))
popinfo <- popinfo %>% arrange(region, subsistence)

# yengo_snps <- fread(paste0(wd, "data/gwas-association-downloaded_2024-10-22-accessionId_GCST90245848.tsv"),
#                     #select = c("CHR_ID", "CHR_POS", "MAPPED_GENE", "UPSTREAM_GENE_ID", "DOWNSTREAM_GENE_ID", "SNP_GENE_IDS", "UPSTREAM_GENE_DISTANCE", "DOWNSTREAM_GENE_DISTANCE", "STRONGEST SNP-RISK ALLELE", "CONTEXT", "INTERGENIC", "RISK ALLELE FREQUENCY", "P-VALUE")
#                     )
#
# snp_mart <- useEnsembl(biomart = "snps", dataset = "hsapiens_snp", version = "GRCh37")
# hg19_pos <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"), filters = "snp_filter", values = yengo_snps$SNPS, mart = snp_mart)
# hg19_names <- hg19_pos %>% filter(chr_name %in% as.character(1:22))
# hg19_names$ID <- paste0("chr", hg19_names$chr_name, "pos", hg19_names$chrom_start)
# saveRDS(hg19_names, paste0(wd, "data/yengo_hg19_names.rds"))
hg19_names <- readRDS(paste0(wd, "data/yengo_hg19_names.rds"))

# Get iHS outputs for populations for intersecting snps
ihs_files <- list.files("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/inputs", pattern = "iHS", full.names = T)
names(ihs_files) <- lapply(basename(ihs_files), function(x) strsplit(x, "_")[[1]][1])
ihs_height <- lapply(names(ihs_files), function(x) fread(ihs_files[`x`], select = c("ID", "Std iHS")) %>% filter(ID %in% hg19_names$ID) %>% rename({{x}} := "Std iHS"))
ihs_height <- purrr::reduce(ihs_height, full_join, by = "ID")
saveRDS(ihs_height, paste0(wd, "data/yengo_height_ihs.rds")) # 11,483 unique SNPs
ihs_height <- readRDS(paste0(wd, "data/yengo_height_ihs.rds"))

# Get global distribution of iHS
load(paste0(wd, "data/ihs_summaries.RData"))

# Function to get percentiles per population with respect to each per-population abs iHS distribution
ihs_percentile_proportions <- function(snps_ihs, cutoffs) {
  proportions <- lapply(cutoffs, function(x)
    data.frame(
      cutoff = x,
      pop = sort(popinfo$pop_code),
      proportion = colSums(
        abs(snps_ihs[,2:35]) >
        sapply(abs_summaries[1:34], function(y) quantile(y, x)),
        na.rm = T) / colSums(!is.na(snps_ihs[,2:35]))
    ))
  proportions <- bind_rows(proportions)
  proportions <- proportions %>% left_join(popinfo, by = c("pop" = "pop_code"))
  proportions$cutoff <- as.factor(proportions$cutoff)
  return(proportions)
}

cutoffs <- c(0.5, 0.75, 0.9, 0.95)
height_proportions <- ihs_percentile_proportions(ihs_height, cutoffs)

pdf(paste0(wd, "figures/height_ihs_barplots.pdf"), height = 8, width = 10)
ggplot(height_proportions, aes(y = pop, x = proportion, fill = cutoff)) +
  geom_col(position = "identity") +
  geom_vline(xintercept = 1-cutoffs, linetype = 2) +
  facet_wrap(~region + subsistence, scales = "free_y") +
  labs(y = "", x = "Proportion of SNPs", fill = "Global |iHS| percentile", title = "Height GWAS SNPs") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "right")
dev.off()

# Plot iHS distribution vs background per population
ihs_height <- ihs_height %>% pivot_longer(cols = popinfo$pop_code, names_to = "pop_code", values_to = "ihs") %>% left_join(popinfo, by = "pop_code") %>% mutate(abs_ihs = abs(ihs))
# ihs_height$continent <- ifelse(grepl("Africa", ihs_height$region), "Africa", "Asia")
ihs_height <- na.omit(ihs_height)

# Tails by population or region
pop_ihs_tails <- ihs_height %>% group_by(pop_code) %>% summarize(top10 = quantile(abs_ihs, probs = 0.9),
                                                                 top5 = quantile(abs_ihs, probs = 0.95))
pop_ihs_tails <- column_to_rownames(pop_ihs_tails, var = "pop_code")

ihs_height_top10 <- list()
for (p in popinfo$pop_code) {
  ihs_height_top10[[`p`]] <- ihs_height %>% filter(pop_code == p & abs_ihs > pop_ihs_tails[p, "top10"])
}
ihs_height_top10 <- bind_rows(ihs_height_top10)
ihs_height_top10$pop_code <- factor(ihs_height_top10$pop_code, levels = popinfo$pop_code)
ihs_height_top10$region <- factor(ihs_height_top10$region, levels = levels(popinfo$region))

# No difference per-region between global and height iHS
# Look at specifically the tail e.g. 95 percentile and whether they're higher between hg/ag
pdf(paste0(wd, "figures/yengo_snp_ihs_pop_tail.pdf"), height = 8, width = 12)
ihs_height_top10 %>% filter(subsistence != "Pastoralist") %>%
  ggplot(aes(x = pop_code, y = abs_ihs, fill = subsistence, color = subsistence)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  scale_discrete_manual(aesthetics = c("color", "fill"), values = c("#F8766D", "#00BA38")) +
  #geom_point(position = "jitter") +
  geom_hline(yintercept = 1.57944, linetype = 1) +
  geom_hline(yintercept = 1.93818, linetype = 2) +
  geom_hline(yintercept = 1.72758, linetype = 2, color = "darkgray") +
  geom_hline(yintercept = 2.31746, linetype = 2, color = "darkgray") +
  labs(x = "Region", y = "|iHS|") +
  facet_nested(~region, scales = "free_x") +
  theme_classic() +
  theme(strip.background = element_blank(),
        #strip.text.y.right = element_text(angle = 0, hjust = 0),
        panel.spacing = unit(0, "cm"),
        axis.text.x = element_text(angle = 90))
dev.off()

```

**Step 4) Metabolic disorder GWAS SNPs**

```bash
metabolic_disorder_snps <- fread(paste0(wd, "data/gwas-association-downloaded_2025-01-31-EFO_0000589-withChildTraits.tsv")) # 17,401 SNPs

# # Map GWAS SNPs to hg19
# snp_mart <- useEnsembl(biomart = "snps", dataset = "hsapiens_snp", version = "GRCh37")
# hg19_pos <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"), filters = "snp_filter", values = metabolic_disorder_snps$SNPS, mart = snp_mart)
# hg19_names <- hg19_pos %>% filter(chr_name %in% as.character(1:22))
# hg19_names$ID <- paste0("chr", hg19_names$chr_name, "pos", hg19_names$chrom_start)
# hg19_names <- hg19_names[!duplicated(hg19_names$ID),]
# saveRDS(hg19_names, paste0(wd, "data/metabolic_hg19_names.rds")) # 10,307 SNPs mapped to hg19
hg19_names <- readRDS(paste0(wd, "data/metabolic_hg19_names.rds"))

# Get iHS outputs for populations for intersecting snps
ihs_files <- list.files("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/inputs", pattern = "iHS", full.names = T)
names(ihs_files) <- lapply(basename(ihs_files), function(x) strsplit(x, "_")[[1]][1])
ihs_metabolic <- lapply(names(ihs_files), function(x) fread(ihs_files[`x`], select = c("ID", "Std iHS")) %>% filter(ID %in% hg19_names$ID) %>% rename({{x}} := "Std iHS"))
ihs_metabolic <- purrr::reduce(ihs_metabolic, full_join, by = "ID")
ihs_metabolic <- ihs_metabolic %>% left_join(hg19_names[,c(1,4)], by = "ID", relationship = "one-to-one")
saveRDS(ihs_metabolic, paste0(wd, "data/metabolic_disorder_ihs.rds")) # 8,696 unique SNPs
ihs_metabolic <- readRDS(paste0(wd, "data/metabolic_disorder_ihs.rds"))

# Filter just for t1d and t2d
diabetes_snps <- metabolic_disorder_snps$SNPS[grepl("diabetes|diabetic", metabolic_disorder_snps$`DISEASE/TRAIT`, ignore.case = T)] %>% unique() # 7,098
ihs_diabetes <- ihs_metabolic %>% filter(refsnp_id %in% diabetes_snps) # 4,864

obesity_snps <- metabolic_disorder_snps$SNPS[grepl("obesity", metabolic_disorder_snps$`DISEASE/TRAIT`, ignore.case = T)] %>% unique() # 614
ihs_obesity <- ihs_metabolic %>% filter(refsnp_id %in% obesity_snps) # 374

# Calculate proportion of SNPs that exceed the global thresholds per-population for barplot
cutoffs <- c(0.5, 0.75, 0.9, 0.95)

metabolic_proportions <- ihs_percentile_proportions(ihs_metabolic, cutoffs)
diabetes_proportions <- ihs_percentile_proportions(ihs_diabetes, cutoffs)
obesity_proportions <- ihs_percentile_proportions(ihs_obesity, cutoffs)

pdf(paste0(wd, "figures/metabolic_ihs_barplots.pdf"), height = 8, width = 10)
ggplot(metabolic_proportions, aes(y = pop, x = proportion, fill = cutoff)) +
  geom_col(position = "identity") +
  geom_vline(xintercept = 1-cutoffs, linetype = 2) +
  facet_wrap(~region + subsistence, scales = "free_y") +
  labs(y = "", x = "Proportion of SNPs", fill = "Global |iHS| percentile", title = "Metabolic disorder GWAS SNPs") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "right")

ggplot(diabetes_proportions, aes(y = pop, x = proportion, fill = cutoff)) +
  geom_col(position = "identity") +
  geom_vline(xintercept = 1-cutoffs, linetype = 2) +
  facet_wrap(~region + subsistence, scales = "free_y") +
  labs(y = "", x = "Proportion of SNPs", fill = "Global |iHS| percentile", title = "Diabetes-related GWAS SNPs") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "right")

ggplot(obesity_proportions, aes(y = pop, x = proportion, fill = cutoff)) +
  geom_col(position = "identity") +
  geom_vline(xintercept = 1-cutoffs, linetype = 2) +
  facet_wrap(~region + subsistence, scales = "free_y") +
  labs(y = "", x = "Proportion of SNPs", fill = "Global |iHS| percentile", title = "Obesity-related GWAS SNPs") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "right")
dev.off()

```

**Step 5) Look at enrichment for SNPs from reproductive traits**

```r
reproductive_snps <- lapply(paste0(wd, "data/",
                                   c("EFO_0004703-age_at_menarche.tsv",
                                     "EFO_0004704-age_at_menopause.tsv",
                                     "EFO_0005677-puberty_onset_measurement.tsv",
                                     "EFO_0006923-fertility_measurement.tsv",
                                     "EFO_0011008-sex_hormone_measurement.tsv")),
                            function(x) fread(x) %>% mutate(across(everything(), as.character)))
reproductive_snps <- bind_rows(reproductive_snps) # 10,377 SNPs
#
# # Map GWAS SNPs to hg19
# snp_mart <- useEnsembl(biomart = "snps", dataset = "hsapiens_snp", version = "GRCh37")
# hg19_pos <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"), filters = "snp_filter", values = reproductive_snps$SNPS, mart = snp_mart)
# hg19_names <- hg19_pos %>% filter(chr_name %in% as.character(1:22))
# hg19_names$ID <- paste0("chr", hg19_names$chr_name, "pos", hg19_names$chrom_start)
# hg19_names <- hg19_names[!duplicated(hg19_names$ID),]
# saveRDS(hg19_names, paste0(wd, "data/reproductive_hg19_names.rds")) # 4,481 SNPs mapped to hg19
hg19_names <- readRDS(paste0(wd, "data/reproductive_hg19_names.rds"))

# Get iHS outputs for populations for intersecting snps
ihs_files <- list.files("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/inputs", pattern = "iHS", full.names = T)
names(ihs_files) <- lapply(basename(ihs_files), function(x) strsplit(x, "_")[[1]][1])
ihs_reproductive <- lapply(names(ihs_files), function(x) fread(ihs_files[`x`], select = c("ID", "Std iHS")) %>% filter(ID %in% hg19_names$ID) %>% rename({{x}} := "Std iHS"))
ihs_reproductive <- purrr::reduce(ihs_reproductive, full_join, by = "ID")
ihs_reproductive <- ihs_reproductive %>% left_join(hg19_names[,c(1,4)], by = "ID", relationship = "one-to-one")
saveRDS(ihs_reproductive, paste0(wd, "data/reproductive_ihs.rds")) # 3,378 unique SNPs
ihs_reproductive <- readRDS(paste0(wd, "data/reproductive_ihs.rds"))

# Filter just for menarche/menopausal and sex hormone levels
menstrual_snps <- reproductive_snps$SNPS[grepl("menarche|menopause|menstrual", reproductive_snps$`DISEASE/TRAIT`, ignore.case = T)] %>% unique() # 924
ihs_menstrual <- ihs_reproductive %>% filter(refsnp_id %in% menstrual_snps) # 803

hormone_snps <- reproductive_snps$SNPS[grepl("levels", reproductive_snps$`DISEASE/TRAIT`, ignore.case = T)] %>% unique() # 3,793
ihs_hormone <- ihs_reproductive %>% filter(refsnp_id %in% hormone_snps) # 2,243

# Calculate proportion of SNPs that exceed the global thresholds per-population for barplot
cutoffs <- c(0.5, 0.75, 0.9, 0.95)

reproductive_proportions <- ihs_percentile_proportions(ihs_reproductive, cutoffs)
menstrual_proportions <- ihs_percentile_proportions(ihs_menstrual, cutoffs)
hormone_proportions <- ihs_percentile_proportions(ihs_hormone, cutoffs)

pdf(paste0(wd, "figures/reproductive_ihs_barplots.pdf"), height = 8, width = 10)
ggplot(reproductive_proportions, aes(y = pop, x = proportion, fill = cutoff)) +
  geom_col(position = "identity") +
  geom_vline(xintercept = 1-cutoffs, linetype = 2) +
  facet_wrap(~region + subsistence, scales = "free_y") +
  labs(y = "", x = "Proportion of SNPs", fill = "Global |iHS| percentile", title = "Reproduction and fertility GWAS SNPs") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "right")

ggplot(menstrual_proportions, aes(y = pop, x = proportion, fill = cutoff)) +
  geom_col(position = "identity") +
  geom_vline(xintercept = 1-cutoffs, linetype = 2) +
  facet_wrap(~region + subsistence, scales = "free_y") +
  labs(y = "", x = "Proportion of SNPs", fill = "Global |iHS| percentile", title = "Menstrual-related GWAS SNPs") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "right")

ggplot(hormone_proportions, aes(y = pop, x = proportion, fill = cutoff)) +
  geom_col(position = "identity") +
  geom_vline(xintercept = 1-cutoffs, linetype = 2) +
  facet_wrap(~region + subsistence, scales = "free_y") +
  labs(y = "", x = "Proportion of SNPs", fill = "Global |iHS| percentile", title = "Sex hormone-related GWAS SNPs") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "right")
dev.off()

```

**Step 6) VIP SNP distribution**

```bash
vip_genes <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze1/metadata/from_providers/viral_effects_file_april2023.txt", header = FALSE)

# Map GWAS SNPs to hg19
snp_mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "GRCh37")
hg19_pos <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"), filters = "ensembl_gene_id", values = vip_genes$V1, mart = snp_mart) # 5661/5892 genes
hg19_names <- hg19_pos %>% filter(chromosome_name %in% as.character(1:22))
hg19_names <- makeGRangesFromDataFrame(hg19_names, keep.extra.columns = T, ignore.strand = T, start.field = "start_position", end.field = "end_position")
saveRDS(hg19_names, paste0(wd, "data/vips_hg19_names.rds")) # 5,221 genes mapped to hg19
hg19_names <- readRDS(paste0(wd, "data/vips_hg19_names.rds"))

# Get iHS outputs for populations for intersecting snps
ihs_files <- list.files("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/inputs", pattern = "iHS", full.names = T)
names(ihs_files) <- lapply(basename(ihs_files), function(x) strsplit(x, "_")[[1]][1])
ihs_vip <- lapply(names(ihs_files), function(x) {
  tmp <- fread(ihs_files[`x`], select = c("ID", "Std iHS"))
  chrpos <- str_extract_all(tmp$ID, "\\d+")
  tmp$chr <- sapply(chrpos, function(x) x[1])
  tmp$start <- sapply(chrpos, function(x) x[2])
  tmp$end <- sapply(chrpos, function(x) x[2])
  tmp <- makeGRangesFromDataFrame(tmp, ignore.strand = T, keep.extra.columns = T)
  res <- find_overlaps(hg19_names, tmp)
  tmp <- as.data.frame(res) %>% rename({{x}} := "Std.iHS") %>% dplyr::select(6:7)
  })
ihs_vip <- lapply(ihs_vip, function(x) unique(na.omit(x))) # remove NA entries, SNPs overlapping multiple genes
ihs_vip <- purrr::reduce(ihs_vip, full_join, by = "ID")
saveRDS(ihs_vip, paste0(wd, "data/vip_ihs.rds")) # 1,961,037 unique SNPs
ihs_vip <- readRDS(paste0(wd, "data/vip_ihs.rds"))

# Calculate proportion of SNPs that exceed the global thresholds per-population for barplot
cutoffs <- c(0.5, 0.75, 0.9, 0.95)
vip_proportions <- ihs_percentile_proportions(ihs_vip, cutoffs)

pdf(paste0(wd, "figures/vip_ihs_barplots.pdf"), height = 8, width = 10)
ggplot(vip_proportions, aes(y = pop, x = proportion, fill = cutoff)) +
  geom_col(position = "identity") +
  geom_vline(xintercept = 1-cutoffs, linetype = 2) +
  facet_wrap(~region + subsistence, scales = "free_y") +
  labs(y = "", x = "Proportion of SNPs", fill = "Global |iHS| percentile", title = "VIP-overlapping SNPs") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "right")
dev.off()

```

**Step 7) ALDH SNP distribution**

```bash
hg19_genes <- readRDS("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/gs_data/hg19_genes.rds")
background <- unique(sort(hg19_genes$gene_name))
go_kegg_hallmark_term2gene <- readRDS("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/gs_data/go_kegg_hallmark_term2gene.rds")
go_kegg_hallmark_term2name <- readRDS("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/gs_data/go_kegg_hallmark_term2name.rds")
go_child2 <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/gs_data/GO_child2.csv", header = F)
colnames(go_child2) <- c("term_id", "term_name")

adh_genes <- go_kegg_hallmark_term2gene %>% filter(term == "GO:0004022") # 9 genes
alc_metabolism_genes <- go_kegg_hallmark_term2gene %>% filter(term == "GO:0006066") # 371 genes

# Map GWAS SNPs to hg19
hg19_adh <- hg19_genes %>% filter(gene_name %in% adh_genes$gene)
hg19_alc <- hg19_genes %>% filter(gene_name %in% alc_metabolism_genes$gene) # 357/371

# Get iHS outputs for populations for intersecting snps
ihs_files <- list.files("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/new_iHS/inputs", pattern = "iHS", full.names = T)
names(ihs_files) <- lapply(basename(ihs_files), function(x) strsplit(x, "_")[[1]][1])

ihs_alc <- lapply(names(ihs_files), function(x) {
  tmp <- fread(ihs_files[`x`], select = c("ID", "Std iHS"))
  chrpos <- str_extract_all(tmp$ID, "\\d+")
  tmp$chr <- sapply(chrpos, function(x) x[1])
  tmp$start <- sapply(chrpos, function(x) x[2])
  tmp$end <- sapply(chrpos, function(x) x[2])
  tmp <- makeGRangesFromDataFrame(tmp, ignore.strand = T, keep.extra.columns = T)
  res <- find_overlaps(hg19_alc, tmp)
  tmp <- as.data.frame(res) %>% rename({{x}} := "Std.iHS") #%>% dplyr::select(8:9)
  })

ihs_adh <- lapply(ihs_alc, function(x) filter(x, gene_name %in% hg19_adh@elementMetadata$gene_name))
ihs_adh <- lapply(ihs_adh, function(x) unique(na.omit(x[,8:9]))) # remove NA entries, SNPs overlapping multiple genes
ihs_adh <- purrr::reduce(ihs_adh, full_join, by = "ID")

ihs_alc <- lapply(ihs_alc, function(x) unique(na.omit(x[,8:9]))) # remove NA entries, SNPs overlapping multiple genes
ihs_alc <- purrr::reduce(ihs_alc, full_join, by = "ID")
saveRDS(ihs_alc, paste0(wd, "data/alcohol_metabolism_ihs.rds")) # 134,671 unique SNPs

# Calculate proportion of SNPs that exceed the global thresholds per-population for barplot
cutoffs <- c(0.5, 0.75, 0.9, 0.95)
adh_proportions <- ihs_percentile_proportions(ihs_adh, cutoffs)
alc_proportions <- ihs_percentile_proportions(ihs_alc, cutoffs)

pdf(paste0(wd, "figures/alc_ihs_barplots.pdf"), height = 8, width = 10)
ggplot(adh_proportions, aes(y = pop, x = proportion, fill = cutoff)) +
  geom_col(position = "identity") +
  geom_vline(xintercept = 1-cutoffs, linetype = 2) +
  facet_wrap(~region + subsistence, scales = "free_y") +
  labs(y = "", x = "Proportion of SNPs", fill = "|iHS| percentile", title = "ADH-overlapping SNPs") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "right")

ggplot(alc_proportions, aes(y = pop, x = proportion, fill = cutoff)) +
  geom_col(position = "identity") +
  geom_vline(xintercept = 1-cutoffs, linetype = 2) +
  facet_wrap(~region + subsistence, scales = "free_y") +
  labs(y = "", x = "Proportion of SNPs", fill = "|iHS| percentile", title = "Alcohol metabolism-overlapping SNPs") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "right")
dev.off()

```

[Comparisons with Akbari et al. ](Comparisons%20with%20Akbari%20et%20al%2017c749e8fe65801c98e6ca94cff731f3.md)

[Immune regulatory annotations](Immune%20regulatory%20annotations%2017c749e8fe6580b18984f2ace36706be.md)

[Brain regulatory annotations](Brain%20regulatory%20annotations%2017c749e8fe6580dcb3c1c7216109ce41.md)