# 2) Ancestry and demographic analyses 

This page details the code needed to generate PCA, ADMIXTURE, MSMC2, and ROH results from the Level 3 VCF, visualized in Figures 1 and 2 and Supplemental Figures 2, 3, 4, and 5. 

## Population map 

**Step 1) Create map of sampling locations for all populations to create Figure 1a.**

```r
######################
####    READ ME   ####
######################
# Generating figure for dist. study pops worldwide

library(ggmap)
library(ggplot2)
library(ggforce)
library(dplyr)
library(tidygeocoder)

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/paper/figs")
freeze2_metadata_v3 <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v3.txt")
freeze2_popinfo <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")


#########################
####    CLEAN DATA   ####
#########################
sample_counts <- aggregate(study_ID ~ pop_code, data = freeze2_metadata_v2, FUN = length)
sample_counts <- merge(sample_counts, freeze2_popinfo, by = "pop_code")

# Plot proportion of individuals w/each subsistence by country
summary_data <- sample_counts %>%
  group_by(country, subsistence) %>%
  summarise(count = sum(study_ID), .groups = "drop") %>%
  group_by(country) %>%
  mutate(total_by_country = sum(count), proportion = count / total_by_country)

# Geocode countries
coordinates <- summary_data %>%
  distinct(country) %>%
  geocode(country, method = "osm")

# Merge coordinates with summary data
summary_data <- merge(summary_data, coordinates, by = "country", all.x = TRUE)
summary_data <- summary_data %>%
  mutate(end_angle = cumsum(proportion) * 2 * pi,  # Cumulative end angle
         start_angle = lag(end_angle, default = 0))

####################
####    PLOT    ####
####################
world_map <- borders("world", colour = "gray50", fill = "gray80") # Base world map

ggplot() +
  world_map +
  geom_arc_bar(data = summary_data,
               aes(x0 = long, y0 = lat, r0 = 0, r = 2, start = start_angle, end = end_angle, fill = subsistence),
               alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("Hunter-gatherer" = "blue", "Agriculturalist" = "green")) +
  theme_void() +
  labs(title = "Proportion of Hunter-Gatherers vs Agriculturalists by Country")


save.image("world_pops.RData")
```

## Global and continental PCA

**Step 1) Convert VCF to PLINK binary format.** 

`plink_makebed.sh` → `freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.bed`

```bash
#!/bin/bash

# Load Level 3 VCF into PLINK format
plink --vcf freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.vcf.gz --make-bed --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin

```

**Step 2) Update PLINK famID to population labels.**

`rename_fam.sh` → `freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.fam`

```r
# Rename first column using population labels
ref <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v3.txt")
fam <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.fam")
fam$V1 <- ref$pop_code[match(fam$V1, ref$study_ID)]
fam$V6 <- 1

write.table(fam, file = "freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.fam", col.names =F, row.names=F, quote=F)
```

**Step 3) Filter BED by MAF > 0.01, and run PCA for all samples (global) as well as each continent (Africa, Asia, Oceania).**

`plink_pca.sh` → `freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_${CONT}.{eigenval, eigenvec}`

```bash
#!/bin/bash

# Filter by MAF > 0.01
plink --bfile freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin --maf 0.01 --make-bed --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_global

# Subset to individuals within continent
for CONT in africa asia oceania; do
  mkdir ${CONT}
  plink --bfile freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_global --keep ${CONT}.fam --make-bed --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_${CONT}
done

# Run PCA for each continent
for CONT in global africa asia oceania; do
  plink --bfile freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_${CONT} --pca --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_${CONT}
done
```

**Step 4) Plot global PCA to create Figure 1b.**

`freeze2_global_pca.R`

```r
#######################
####    READ ME    ####
#######################

# Script for plotting global pca for Fig1 of the manuscript

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/pca")

library(tidyverse)
library(ggforce)
library(ggrepel)
library(ggh4x)
library(plotly)
library(gridExtra)
library(patchwork)
library(colorspace)
library(paletteer)
library(data.table)
library(here)
library(concaveman)
library(cowplot)
library(ggnewscale)

out <- here("Fig1_PCA_plots")
dir.create(out)

eigenvec <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf.eigenvec",sep=" ",header=F)
colnames(eigenvec) <- c("pop_code","study_ID", paste0("PC",c(1:(ncol(eigenvec)-2))))

eigenval <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf.eigenval",sep=" ",header=F)
colnames(eigenval) <- paste0("eigenval")

# calculate percent variance explained
pve <- (eigenval/sum(eigenval)) * 100
pve <- round(pve, 2)

meta <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v3.txt")

# define colours in the region
region_colors <- c(
  #Oceania
  sequential_hcl(n = 5, palette = "Blues 2")[c(2)],
  #Asia - 3 subregions
  sequential_hcl(n = 14, palette = "Reds")[c(2,5,9)], #"#DD6157" "#E99A2C" "#E6D050"
  #Africa - 4 subregions
  sequential_hcl(n = 15, palette = "Mint")[c(1,5,10,15)]
)

names(region_colors) <- c("Oceania",
                          "East Asia",
                          "Southeast Asia",
                          "South Asia",
                          "East Africa",
                          "East Central Africa",
                          "West Central Africa",
                          "Southern Africa"
)

pca_dat <- merge(eigenvec, meta, by = c("study_ID", "pop_code"),) %>% mutate(region = factor(region, levels=names(region_colors)))

subsistence_shape <- c(25, 21, 22)
names(subsistence_shape) <- unique(pca_dat$subsistence)

# define label colors
labs <- filter(pca_dat, !duplicated(pca_dat$region)) %>% select(region)
rgbs <- col2rgb(region_colors[labs$region])
labs$txt_color <- apply(rgbs, 2, function(x) {
  as.character(x[1] * 0.299 + x[2] * 0.587 + x[3] * 0.114 <= 150)
})
labs$rtxt <- as.character(labs$region)
labs$rtxt[c(4,7)] <- c("West Central\nAfrica", "East Central\nAfrica")
labs$x <- c(-0.018, -0.03, -0.03, 0.05, 0.018, 0.03, 0.05, -0.035)
labs$y <- c(0.01, -0.01, -0.02, -0.035, 0.0, 0.11, -0.015, 0.01)

subs_region_bar <- pca_dat %>%
  mutate(region= factor(region, levels = c(names(region_colors)))) %>%
  ggplot(aes(y= fct_rev(subsistence), fill = region)) +
  geom_bar(color= "darkgray") +
  scale_fill_manual(values = region_colors) +
  ylab("") + xlab("Number of samples") +
  theme_classic(base_size = 7) +
  theme(legend.position = "none")
ggsave("test.pdf", subs_region_bar, "pdf", width = 4, height = 4, units = "in")
bar_grob <- ggplotGrob(subs_region_bar)

main_pc_plot <- pca_dat %>%
  ggplot(aes(x= PC1, y=PC2, fill = region)) +
  geom_point(aes(shape = subsistence), color= "black",
             size= 1, alpha = 0.7) +
  annotation_custom(bar_grob, xmin = -0.04, xmax = 0.015, ymin = 0.05, ymax = 0.125) +
  scale_shape_manual(values = subsistence_shape) +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  new_scale_color() +
  geom_label(aes(x = x, y = y, label = rtxt, color = txt_color), data = labs,
             size = 5, size.unit = "pt", linewidth = 0) +
  labs(x = NULL, y = NULL) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white")) +
  theme_classic(base_size = 7) +
  theme(legend.position = "none") +
  labs(x = paste0("PC1 (", pve[1,], "%)"), y = paste0("PC2 (", pve[2,], "%)"))
```

**Step 5) Plot continental PCAs to create Figure 1c.**

`freeze2_regional_pcas.R`

```r
#######################
####    READ ME    ####
#######################
# Regional PCAs
# Baseline (QC+ 10% site missingness), regional ops except inds with outlier heterozygosity and/or related, 1% maf

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/pca")

library(tidyverse)
library(cowplot)
library(gridExtra)

# read in common metadata for merging with each continental PCA df
metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v3.txt")

# SHARED VARIABLES
shapes_for_subsistence <- c("Hunter-gatherer" = 25,
                            "Agriculturalist" = 21,
                            "Pastoralist" = 23)

# define colours for population
colours_for_pop <- c("AET" = "firebrick1", "AGT" = "firebrick4", "ATI" = "brown1", "PHP" = "coral", "RAM" = "darkred",
                     "KWI" = "orangered2", "MWA" = "deeppink", "MNO" = "indianred1", "MTW" = "indianred3", "SMB" = "maroon",
                     "ADV" = "chocolate1", "SHL" = "sienna4", "ATY" = "magenta", "PAI" = "orchid",
                     "PNG" = "slateblue2", "VAN" = "lightblue2", "POL" = "mediumblue", "SCI" = "skyblue4", "SLI" = "navy",
                     "BCH" = "mediumturquoise", "CBU" = "darkgreen", "SHK" = "darkseagreen",
                     "KIG" = "lemonchiffon3", "TWA" = "green",
                     "HMB" = "olivedrab3", "NAM" = "mediumaquamarine", "SAN" = "lightseagreen",
                     "BAK" = "turquoise4", "FAN" = "yellowgreen", "GAL" = "palegreen4",
                     "KOT" = "palegreen", "NZE" = "olivedrab1", "TEK" = "mediumspringgreen", "TSO" = "khaki4")

PC_names <- paste0("PC_",c(1:20))

########################
####    OCEANIA     ####
########################
# read in plink eigenvec and eigenval inputs
eigenvec <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_oceania_0.01_maf.eigenvec",sep=" ",header=F)
colnames(eigenvec) <- c("population_code", "sample_name", PC_names)
eigenvec <- eigenvec[order(eigenvec$sample_name),]
eigenval <- scan("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_oceania_0.01_maf.eigenval")

# calculate percent variance explained
pve <- (eigenval/sum(eigenval)) * 100
pve <- round(eigenval, 2)

# merge with imported metadata
oceania <- merge(eigenvec, metadata, by.x = "sample_name", by.y = "study_ID", all.x = TRUE)

# Ensure the dataset includes the pop_sub_label column and regions are ordered
oceania <- oceania %>%
  select(sample, PC_1, PC_2, pop_code, subsistence) %>% 
  mutate(
    subsistence = factor(subsistence), # Ensure regions are ordered
    cont = "Oceania"
  ) %>% arrange(subsistence, population_code)


oceania_plot <- ggplot(filter(all_conts, cont == "Oceania"), aes(PC_1, PC_2, color = pop_code, shape = subsistence, fill = pop_code)) +
  geom_point(size = 1, alpha = 0.7) + # Adjust size if needed
  scale_shape_manual(values = shapes_for_subsistence) +
  scale_color_manual(values = colours_for_pop) +
  scale_fill_manual(values = colours_for_pop) +
  labs(x = paste0("PC1 (", pve[1], "%)"),
       y = paste0("PC2 (", pve[2], "%)"),
       color = NULL, shape = NULL, fill = NULL) +
  theme_classic(base_size = 7) +
  facet_wrap(~ cont) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )


######################
####     ASIA     ####
######################

# read in plink eigenvec and eigenval inputs
eigenvec <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_asia_0.01_maf.eigenvec",sep=" ",header=F)
PC_names <- paste0("PC_",c(1:20))
colnames(eigenvec) <- c("population_code", "sample_name", PC_names)
eigenvec <- eigenvec[order(eigenvec$sample_name),]
eigenval <- scan("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_asia_0.01_maf.eigenval")

# calculate percent variance explained
pve <- (eigenval/sum(eigenval)) * 100
pve <- round(eigenval, 2)

# merge with imported metadata
asia <- merge(eigenvec, metadata, by.x = "sample_name", by.y = "study_ID", all.x = TRUE)

# Ensure the dataset includes the pop_sub_label column and regions are ordered
asia <- asia %>%
  select(sample, PC_1, PC_2, pop_code, subsistence) %>% 
  mutate(
    subsistence = factor(subsistence), # Ensure regions are ordered
    cont = "Asia"
  ) %>% arrange(subsistence, population_code)


asia_plot <- ggplot(filter(all_conts, cont == "Asia"), aes(PC_1, PC_2, color = pop_code, shape = subsistence, fill = pop_code)) +
  geom_point(size = 1, alpha = 0.7) + # Adjust size if needed
  scale_shape_manual(values = shapes_for_subsistence) +
  scale_color_manual(values = colours_for_pop) +
  scale_fill_manual(values = colours_for_pop) +
  labs(x = paste0("PC1 (", pve[1], "%)"),
       y = paste0("PC2 (", pve[2], "%)"),
       color = NULL, shape = NULL, fill = NULL) +
  theme_classic(base_size = 7) +
  facet_wrap(~ cont) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )


########################
####     AFRICA     ####
########################

# read in plink eigenvec and eigenval inputs
eigenvec <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_africa_0.01_maf.eigenvec",sep=" ",header=F)
PC_names <- paste0("PC_",c(1:20))
colnames(eigenvec) <- c("population_code", "sample_name", PC_names)
eigenvec <- eigenvec[order(eigenvec$sample_name),]
eigenval <- scan("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_africa_0.01_maf.eigenval")

# calculate percent variance explained
pve <- (eigenval/sum(eigenval)) * 100
pve <- round(eigenval, 2)

# merge with imported metadata
africa <- merge(eigenvec, metadata, by.x = "sample_name", by.y = "study_ID", all.x = TRUE)

# Ensure the dataset includes the pop_sub_label column and regions are ordered
africa <- africa %>%
  select(sample, PC_1, PC_2, pop_code, subsistence) %>% 
  mutate(
    subsistence = factor(subsistence), # Ensure regions are ordered
    cont = "Africa"
  ) %>% arrange(subsistence, population_code)


africa <- ggplot(filter(all_conts, cont == "Africa"), aes(PC_1, PC_2, color = pop_code, shape = subsistence, fill = pop_code)) +
  geom_point(size = 1, alpha = 0.7) + # Adjust size if needed
  scale_shape_manual(values = shapes_for_subsistence) +
  scale_color_manual(values = colours_for_pop) +
  scale_fill_manual(values = colours_for_pop) +
  labs(x = paste0("PC1 (", pve[1], "%)"),
       y = paste0("PC2 (", pve[2], "%)"),
       color = NULL, shape = NULL, fill = NULL) +
  theme_classic(base_size = 7) +
  facet_wrap(~ cont) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )


######################################
###   COMBINING THE 3 IN A PANEL   ###
######################################
library(cowplot)

# Combine the plots side by side
all_conts <- rbind(africa, asia, oceania)
composite <- plot_grid(asia_plot, oceania_plot, africa_plot, ncol = 3, align = "v")

out <- here("Fig1_PCA_plots")
three_regions_pca_file <- here(out,"three_regions_pca_plot.pdf")
ggsave(three_regions_pca_file, plot = composite, width = 7.8, height = 2.5, units = "in")

save.image("freeze2_regional_pcas_v2.RData")
```

**Step 6) Plot global PCA at the population level to create Supplementary Figure 2.**

`population_global_pca.R`

```r
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(ggtext)
library(plotly)
library(cowplot)

# Set up metadata and colors
meta <- read.delim("freeze2_metadata_v3.txt")
popinfo <- read.delim("freeze2_popinfo.txt")
eigenvec <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf.eigenvec",sep=" ",header=F)
colnames(eigenvec) <- c("pop_code","study_ID", paste0("PC",c(1:(ncol(eigenvec)-2))))
eigenval <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf.eigenval",sep=" ",header=F)
colnames(eigenval) <- paste0("eigenval")

# calculate percent variance explained
pve <- (eigenval/sum(eigenval)) * 100
pve <- round(pve, 2)

popinfo$region <- factor(popinfo$region, 
                         levels = c("Oceania",
                                    "East Asia",
                                    "Southeast Asia",
                                    "South Asia",
                                    "East Africa",
                                    "East Central Africa",
                                    "West Central Africa",
                                    "Southern Africa"))
popinfo <- arrange(popinfo, region)
popinfo$pop_code <- factor(popinfo$pop_code, levels = popinfo$pop_code)

pca_dat <- left_join(eigenvec, meta, by = c("study_ID", "pop_code"))
pca_dat$region <- factor(pca_dat$region, levels = levels(popinfo$region))
pca_dat$pop_code <- factor(pca_dat$pop_code, levels = levels(popinfo$pop_code))
pop_colors <- c("#5A5156","#F6222E","#FE00FA","#16FF32","#3283FE","#FEAF16","#B00068","#1CFFCE","#90AD1C","#2ED9FF","#DEA0FD","#AA0DFE","#F8A19F","#325A9B","#C4451C","#1C8356","#85660D","#B10DA1","#FBE426","#1CBE4F","#FA0087","#F7E1A0","#C075A6","#782AB6","#AAF400","#BDCDFF","#822E1C","#B5EFB5","#7ED7D1","#1C7F93","#D85FF7","#683B79","#66B0FF","#3B00FB") # 36 cols from Polychrome
names(pop_colors) <- sort(as.character(popinfo$pop_code))

subsistence_shape <- c(25, 21, 22)
names(subsistence_shape) <- unique(pca_dat$subsistence)

# admixture colors: colorPallete(13, "#F41C35")
#"#FB00E0" "#90AD1C" "#00FDDF" "#0D60FF"  "#165A00" "#F1E50D" "#26FF22" "#F41C35" "#53C5FF" "#F8D6D7" "#E48232" "#1C7F93" "#BA0DF5"

# Dummy dataframe for labels
labs <- pca_dat[!duplicated(pca_dat$pop_code),]
rgbs <- col2rgb(pop_colors)
labs$txt_color <- apply(rgbs, 2, function(x) {
  x[1] * 0.299 + x[2] * 0.587 + x[3] * 0.114 <= 150
})

# PC1 vs 2 Africa
afr_sub <- c("KIG", "FAN", "GAL", "KOT", "NZE", "TEK", "TSO", "HMB", "BAK", "TWA")
afr_zoom <- pca_dat %>% 
  filter(pop_code %in% afr_sub) %>% 
  ggplot(aes(x = PC1, y = PC2, fill = pop_code)) +
  geom_point(aes(shape = subsistence), size = 1, alpha = 0.7, color = "black") +
  stat_ellipse(aes(color = pop_code), type = "norm", level = 0.85, linetype = 1) +
  scale_shape_manual(values = subsistence_shape) +
  scale_discrete_manual(c("color", "fill"), values = pop_colors) +
  new_scale_color() +
  geom_label(aes(label = population, color = txt_color), 
             data = filter(labs, pop_code %in% afr_sub), 
             size = 6, size.unit = "pt", label.size = 0) +
  labs(x = NULL, y = NULL) +
  scale_color_manual(values = c("black", "white")) +
  theme_classic() +
  theme(plot.background = element_rect(color = "black", fill = 0, linetype = 1),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))
afr_grob <- ggplotGrob(afr_zoom)

ann_coord <- c(-0.038, 0.018, 0.01, 0.13)
rect_coord <- c(0.029, 0.043, -0.04, -0.01)

p1 <- pca_dat %>%  
  ggplot(aes(x = PC1, y = PC2, fill = pop_code)) +
  geom_point(aes(shape = subsistence), size = 2.25, alpha = 0.7, color = "black") +
  annotation_custom(afr_grob, xmin = ann_coord[1], xmax = ann_coord[2], ymin = ann_coord[3], ymax = ann_coord[4]) +
  geom_rect(aes(xmin = rect_coord[1], xmax = rect_coord[2], ymin = rect_coord[3], ymax = rect_coord[4]), 
            data = . %>% slice_head(n = 1), inherit.aes = F,
            color = "black", linewidth = 0.5, alpha = 0) +
  geom_path(
    data = data.frame(
      x = c(ann_coord[1], rect_coord[1], ann_coord[2], rect_coord[2]), 
      y = c(ann_coord[3], rect_coord[3], ann_coord[4], rect_coord[4]),
      grp = c(1, 1, 2, 2)), inherit.aes = F,
      aes(x, y, group = grp), color = "black", linetype = 2) +
  geom_label(aes(label = population, color = txt_color), data = filter(labs, grepl("Africa", region) & !(pop_code %in% afr_sub)), 
             size = 6, size.unit = "pt", label.size = 0) +
  scale_shape_manual(values = subsistence_shape) +
  scale_fill_manual(values = pop_colors) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = "none", fill = "none") +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "top") +
  labs(x = paste0("PC1 (", pve[1,], "%)"), 
       y = paste0("PC2 (", pve[2,], "%)"), shape = "")

# PC1 vs 3 Asia/Oceania
oce_sub <- c("PNG", "SCI", "SLI", "VAN", "KWI")
oce_zoom <- pca_dat %>% 
  filter(pop_code %in% oce_sub) %>%
  ggplot(aes(x = PC1, y = PC3, fill = pop_code)) +
  geom_point(aes(shape = subsistence), size = 1, alpha = 0.7, color = "black") +
  stat_ellipse(aes(color = pop_code), type = "norm", level = 0.85, linetype = 1) +
  scale_shape_manual(values = subsistence_shape) +
  scale_discrete_manual(c("color", "fill"), values = pop_colors) +
  new_scale_color() +
  geom_label(aes(label = population, color = txt_color), 
             data = filter(labs, pop_code %in% oce_sub), 
             size = 6, size.unit = "pt", label.size = 0) +
  scale_color_manual(values = c("black", "white")) +
  scale_x_continuous(n.breaks = 3) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(plot.background = element_rect(color = "black", fill = 0, linetype = 1),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))
oce_grob <- ggplotGrob(oce_zoom)

asi_sub <- c("ATY", "PAI", "AET", "AGT", "ATI", "MNO", "MTW", "MWA", "PHP")
asi_zoom <- pca_dat %>% 
  filter(pop_code %in% asi_sub) %>%
  ggplot(aes(x = PC1, y = PC3, fill = pop_code)) +
  geom_point(aes(shape = subsistence), size = 1, alpha = 0.7, color = "black") +
  stat_ellipse(aes(color = pop_code), type = "norm", level = 0.85, linetype = 1) +
  scale_shape_manual(values = subsistence_shape) +
  scale_discrete_manual(c("color", "fill"), values = pop_colors) +
  new_scale_color() +
  geom_label(aes(label = population, color = txt_color), 
             data = filter(labs, pop_code %in% asi_sub), 
             size = 6, size.unit = "pt", label.size = 0) +
  scale_color_manual(values = c("black", "white")) +
  scale_x_continuous(n.breaks = 3) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(plot.background = element_rect(color = "black", fill = 0, linetype = 1),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))
asi_grob <- ggplotGrob(asi_zoom)

ann_coord1 <- c(-0.025, 0.02, -0.08, 0)
rect_coord1 <- c(-0.038, -0.029, -0.072, -0.007)
ann_coord2 <- c(0, 0.056, 0.02, 0.11)
rect_coord2 <- c(-0.033, -0.026, 0.0305, 0.0665)

p2 <- pca_dat %>% 
  ggplot(aes(x = PC1, y = PC3, fill = pop_code)) +
  geom_point(aes(shape = subsistence), size = 2.25, alpha = 0.7, color = "black", show.legend = F) +
  # OCE
  annotation_custom(oce_grob, xmin = ann_coord1[1], xmax = ann_coord1[2], ymin = ann_coord1[3], ymax = ann_coord1[4]) +
  geom_rect(aes(xmin = rect_coord1[1], xmax = rect_coord1[2], ymin = rect_coord1[3], ymax = rect_coord1[4]), 
            data = . %>% slice_head(n = 1), inherit.aes = F,
            color = "black", linewidth = 0.5, alpha = 0) +
  geom_path(
    data= data.frame(
      x = c(ann_coord1[1], rect_coord1[2], ann_coord1[1], rect_coord1[2]), 
      y = c(ann_coord1[3], rect_coord1[3], ann_coord1[4], rect_coord1[4]),
      grp = c(1, 1, 2, 2)), inherit.aes = F,
      aes(x, y, group = grp), color = "black", linetype = 2) +
  # ASIA
  annotation_custom(asi_grob, xmin = ann_coord2[1], xmax = ann_coord2[2], ymin = ann_coord2[3], ymax = ann_coord2[4]) +
  geom_rect(aes(xmin = rect_coord2[1], xmax = rect_coord2[2], ymin = rect_coord2[3], ymax = rect_coord2[4]),
            data = . %>% slice_head(n = 1), inherit.aes = F,
            color = "black", linewidth = 0.5, alpha = 0) +
  geom_path(
    data= data.frame(
      x = c(ann_coord2[1], rect_coord2[2], ann_coord2[1], rect_coord2[2]), 
      y = c(ann_coord2[3], rect_coord2[3], ann_coord2[4], rect_coord2[4]),
      grp = c(3, 3, 4, 4)), inherit.aes = F,
      aes(x, y, group = grp), color = "black", linetype = 2) +
  
  geom_label(aes(label = population, color = txt_color), 
             data = filter(labs, !grepl("Africa", region) & !(pop_code %in% c(oce_sub, asi_sub))), 
             size = 6, size.unit = "pt", label.size = 0) +
  scale_shape_manual(values = subsistence_shape) +
  scale_fill_manual("fill", values = pop_colors) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = "none", fill = "none") +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "right") +
  labs(x = paste0("PC1 (", pve[1,], "%)"), 
       y = paste0("PC3 (", pve[3,], "%)"), shape = "")

plot <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1.075, 1), labels = c("a", "b"))

pdf("global_pca_v2.pdf", width = 7.5, height = 8)
print(plot)
dev.off()

```

## Global and continental admixture proportions

**Step 1) Starting from PLINK outputs generated above for PCA, LD prune SNPs and subset samples by continent.**

`plink_make_bed.sh` → `freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned_${CONT}`

```bash
#!/bin/bash

module load plink

PRE="freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf"

# Load in previous bims and LD-prune w/ 50 SNPS, slide 10 SNPs, r2 = 0.3 per-window, then make BED of SNPs
plink --bfile ${PRE} --indep-pairwise 50 10 0.3 --out ${PRE}

plink --bfile ${PRE} --extract ${PRE}.prune.in --make-bed --out ${PRE}_pruned_global

# Subset to individuals within continent
for CONT in africa asia oceania; do
  mkdir ${CONT}
  plink --bfile ${PRE}_pruned --keep ${CONT}.fam --make-bed --out ${PRE}_pruned_${CONT}
done 
```

**Step 2) Run ADMIXTURE for K = 1-13 with 10 replicates. Specify the sample subset with the `CONT` variable.**

`admixture.sh ${K} ${r} ${CONT}` → `log.freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned_${CONT}.K${K}.r${r}.out`

```bash
#!/bin/bash

K=$1 # Get the value of K from the first command line argument
r=$2 # Get the value of r from the second command line argument
CONT=$3

INPUTDIR=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/updated/baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned
PRE=freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned_${CONT}

# Run ADMIXTURE at specified K and r 
mkdir -p "${CONT}/R${r}"
/project/lbarreiro/SHARED/PROGRAMS/admixture_linux-1.3.0/admixture -s ${RANDOM} --cv ${INPUTDIR}/${PRE}.bed $K -j16 | tee log.${PRE}.K${K}.r${r}.out
mv ${PRE}.${K}.Q ${PRE}.K${K}r${r}.Q
```

**Step 3) Extract cross-validation error rates from ADMIXTURE output files and prepare ind2pop (i.e. population label) for pong.**

`create_ind2pop.sh` → `CV_error_${CONT}.txt` and ` ${CONT}.ind2pop`

```bash
#!/bin/bash

for CONT in global, africa, asia, oceania; do
  grep -h CV ${CONT}/**/log*.out | awk '{gsub(/[^0-9]+/, "", $3); print $3, $4}' > CV_error_${CONT}.txt`
  cut -d' ' -f1 /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/updated/baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned/freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned_${CONT}.fam > ${CONT}.ind2pop
done 
```

**Step 4) Assemble the Qfilemap file that points pong towards which “.Q” files to display.**

`create_Qmap.sh` → `${CONT}.multiplerun.Qfilemap`

```bash
#!/bin/bash

PRE=freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned

# Function to extract necessary values from each run to Qfilemap
createQmap(){
local r=$1
local K=$2
local CONT=$3
awk -v K=${K} -v r=${r} -v file=${CONT}/R${r}/${PRE}.K${K}r${r} 'BEGIN{ \
printf("K%dr%d\t%d\t%s.Q\n",K,r,K,file)
}' >> ${CONT}.multiplerun.Qfilemap
}

export -f createQmap

# Loop across all K, r, runs
for K in {2..13}; do for r in {1..10}; for CONT in global, africa, asia, oceania; do createQmap ${r} ${K} ${CONT}; done; done; done
```

**Step 5) Install pong through Conda and run via command line to get a visualization in-browser. Additional options `-n` and `-l` specify order of populations and colors for Supplementary Figure 3.**

```bash
#!/bin/bash

pong -m ${CONT}.multiplerun.Qfilemap -i ${CONT}.ind2pop -n ${CONT}.poporder.txt -l ${CONT}.colors.txt -v -o pong_${CONT}
```

**Step 6) Plot ADMIXTURE outputs for each continent manually to create Figure 1d.**

`freeze2_admixture.R`

```r
#####################################
#### Set dir and import packages ####
#####################################

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/admixture/manuscript_plots") # Copy all files to be plotted here
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)

#############################################
#### Import cross-validation error rates ####
#############################################

CV_error_asia <- read.table("CV_error_asia.txt", quote="\"", comment.char="")
names <- c("K", "cross_validation_error")
colnames(CV_error_asia) <- names

ggplot(CV_error_asia, aes(x = K, y = cross_validation_error)) +
  geom_point() + ylab("Cross validation error") +
  ggtitle("Asia run") + scale_x_continuous(breaks=seq(2, 12, 1))

CV_error_oceania <- read.table("CV_error_oceania.txt", quote="\"", comment.char="")
names <- c("K", "cross_validation_error")
colnames(CV_error_oceania) <- names

ggplot(CV_error_oceania, aes(x = K, y = cross_validation_error)) +
  geom_point() + ylab("Cross validation error") +
  ggtitle("Oceania run") + scale_x_continuous(breaks=seq(2, 12, 1))

CV_error_africa <- read.table("CV_error_africa.txt", quote="\"", comment.char="")
names <- c("K", "cross_validation_error")
colnames(CV_error_africa) <- names

ggplot(CV_error_africa, aes(x = K, y = cross_validation_error)) +
  geom_point() + ylab("Cross validation error") +
  ggtitle("Africa run") + scale_x_continuous(breaks=seq(2, 12, 1))

#####################################
#### Plot ADMIXTURE at lowest CV ####
#####################################

# ASIA (K=3)
q <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned_asia.K3.r2.Q")
k <- 3
m <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned_asia.fam")
colnames(m) <- c("famid", "sampleid", "pat", "mat", "sex", "phenotype")

# combine famid and sampleid with the admixture q matrix
l_df <- cbind(q, m[,c(1,2)])
colnames(l_df) <- c(seq(1:k), 'pop', 'id')
l_df <- l_df %>% drop_na(pop)
asia <- l_df %>% gather(K, value, -id, -pop) # gather/transform df
asia <- asia %>% filter(pop != "ORR")

# order populations using metadata region and country
metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v3.txt")
metadata <- merge(x = metadata, y = l_df, by.x = "study_ID", by.y = "id") # only retaining samples in l_df

region_order <- c("East Asia", "Southeast Asia", "South Asia")

# convert the region column to a factor with the custom order — you wont see the results just yet
metadata$region <- factor(metadata$region, levels = region_order, ordered = TRUE)

# sort the entire metadata frame now based on the custom order of the region column, then extract unique pop_code
metadata <- metadata[order(metadata$region), ]
pop_order <- unique(metadata$pop_code)

# define however many colours there are for K value
pop_color <- c("sienna4", "firebrick1", "orchid")

asia_p <- ggplot(data=asia, aes(x=reorder(id, value, function(x){max(x)}), y=value, fill=reorder(K, sort(as.integer(K))))) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = pop_color) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_x_discrete(expand=c(-1, 0)) +
  facet_grid(. ~ factor(pop, levels=pop_order), scales = "free", space="free", switch="both") +
  theme_classic() +
  theme(panel.spacing = unit(0.2, "lines"),
        strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size = 10, colour = "black", angle = 90, hjust = 0.8),
        strip.placement = "outside",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.position="bottom") +
  ylab("") +
  xlab("") +
  labs(fill="K") +
  guides(fill="none")


# OCEANIA (K=2)
q <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_oceania_0.01_maf_pruned.K2r2.Q")
k <- 2
m <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_oceania_0.01_maf_pruned.fam")
colnames(m) <- c("famid", "sampleid", "pat", "mat", "sex", "phenotype")

# combine famid and sampleid with the admixture q matrix
l_df <- cbind(q, m[,c(1,2)])
colnames(l_df) <- c(seq(1:k), 'pop', 'id')
l_df <- l_df %>% drop_na(pop)
oceania <- l_df %>% gather(K, value, -id, -pop) # gather/transform df

# order populations using metadata region and country
metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v3.txt")
metadata <- merge(x = metadata, y = l_df, by.x = "study_ID", by.y = "id") # only retaining samples in l_df

region_order <- c("Oceania")

# convert the region column to a factor with the custom order — you wont see the results just yet
metadata$region <- factor(metadata$region, levels = region_order, ordered = TRUE)

# sort the entire metadata frame now based on the custom order of the region column, then extract unique pop_code
metadata <- metadata[order(metadata$region), ]
pop_order <- unique(metadata$pop_code)

# define however many colours there are for K value
pop_color <- c("lightblue2", "slateblue2", "black")

oceania_p <- ggplot(data=oceania, aes(x=reorder(id, value, function(x){max(x)}), y=value, fill=reorder(K, sort(as.integer(K))))) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = pop_color) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_x_discrete(expand=c(-1, 0)) +
  facet_grid(. ~ factor(pop, levels=pop_order), scales = "free", space="free", switch="both") +
  theme_classic() +
  theme(panel.spacing = unit(0.2, "lines"),
        strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size = 10, colour = "black", angle = 90, hjust = 0.8),
        strip.placement = "outside",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.position="bottom") +
  ylab("") +
  xlab("") +
  labs(fill="K") +
  guides(fill="none")


# AFRICA (K=4)
q <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_africa_0.01_maf_pruned.K4r2.Q")
k <- 4
m <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_africa_0.01_maf_pruned.fam")
colnames(m) <- c("famid", "sampleid", "pat", "mat", "sex", "phenotype")

# combine famid and sampleid with the admixture q matrix
l_df <- cbind(q, m[,c(1,2)])
colnames(l_df) <- c(seq(1:k), 'pop', 'id')
l_df <- l_df %>% drop_na(pop)
africa <- l_df %>% gather(K, value, -id, -pop) # gather/transform df

# order populations using metadata region and country
metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v3.txt")
metadata <- merge(x = metadata, y = l_df, by.x = "study_ID", by.y = "id") # only retaining samples in l_df

region_order <- c("East Africa", "Southern Africa", "East Central Africa", "West Central Africa")

# convert the region column to a factor with the custom order — you wont see the results just yet
metadata$region <- factor(metadata$region, levels = region_order, ordered = TRUE)

# sort the entire metadata frame now based on the custom order of the region column, then extract unique pop_code
metadata <- metadata[order(metadata$region), ]
pop_order <- unique(metadata$pop_code)

# define however many colours there are for K value
pop_color <- c("darkgreen", "lightseagreen", "green", "khaki4")

africa_init <- ggplot(data=africa, aes(x=reorder(id, value, function(x){max(x)}), y=value, fill=reorder(K, sort(as.integer(K))))) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = pop_color) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_x_discrete(expand=c(0, 0)) +
  facet_grid(. ~ factor(pop, levels=pop_order), scales = "free", space="free", switch="both") +
  theme_classic() +
  theme(panel.spacing = unit(0.2, "lines"),
        strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size = 10, colour = "black", angle = 90, hjust = 0.8),
        strip.placement = "outside",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.position="bottom",
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5)) +
  ylab("") +
  xlab("") +
  labs(fill="K") +
  guides(fill="none")


save.image("regional_admixture_plots.RData")
```

## Historical effective population sizes

**Step 1) We selected 4 individuals with good coverage from each population to estimate population sizes; the following code outlines the pipeline for 2 populations listed below as an example. First, get average sequencing depth and call genotypes at all reference panel sites.**
- KWI001 (37.8x), KWI002 (38.11x), KWI004 (36.46x), KWI005 (36.36x)
- SMB005 (33.84x), SMB010 (37.92x), SMB019 (36.72x), SMB026 (35.98x)

`msmc2_generate_masks_vcf.sh ${sample}` → `${sample}.mask.chr${CHR}.bed.gz` and `${sample}.chr${CHR}.vcf.gz`

```bash
#!/bin/bash

sample=$1 # Get the value of sample from the first command line argument
BAM_dir=/project/lbarreiro/GLOBUS/msmc2_bams
BAM=${BAM_dir}/${sample}.md.bam
reference_fasta=/project/lbarreiro/USERS/bridget/references/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta

# estimating the average sequencing depth using all sites on chromosome 20
DEPTH=$(samtools depth -r 20 ${BAM} | awk '{sum += $3} END {print sum / NR}')

for CHR in {1..22}; do samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $CHR -f ${reference_fasta} ${BAM} | bcftools call -c -V indels | bamCaller.py ${DEPTH} ${sample}.mask.chr${CHR}.bed.gz | gzip -c > ${sample}.chr${CHR}.vcf.gz; done
```

**Step 2) Phase with SHAPEIT2 and reheader output VCFs (can't use SHAPEIT5 due to need for AC field).**

`phase_indv_samples_shapeit2_per_chr.sh ${sample} ${chr}` → `${sample}.chr${chr}.shapeit2.1KGref.vcf.gz`

```bash
#!/bin/bash

sample=${sample}  # Use the sample passed from the master script
chr=${chr}        # Use the chromosome number passed from the master script

output_dir="/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/msmc2/inputs/${sample}"
temp_dir=${output_dir}/chr${chr}  # Temporary directory for the job

# Create the output and temp directories if they do not exist
mkdir -p ${output_dir}
mkdir -p ${temp_dir}

log_align=${temp_dir}/${sample}.chr${chr}.alignments
snp_exclude_list=${log_align}.snp.strand.exclude
log_main=${temp_dir}/${sample}.chr${chr}.main

vcf=${sample}.chr${chr}.vcf.gz
phased_haps=${temp_dir}/${sample}.chr${chr}.shapeit2.1KGref.haps
phased_sample=${temp_dir}/${sample}.chr${chr}.shapeit2.1KGref.sample
phased_vcf=${temp_dir}/${sample}.chr${chr}.shapeit2.1KGref.vcf
log_convert=${temp_dir}/${sample}.chr${chr}.convert

# initial shapeit2 checks
/project/lbarreiro/SHARED/PROGRAMS/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit -check \
--input-vcf ${vcf} \
-M /project/lbarreiro/USERS/bridget/references/maps/impute2-1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt \
--input-ref /project/lbarreiro/USERS/bridget/references/maps/impute2-1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz /project/lbarreiro/USERS/bridget/references/maps/impute2-1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz /project/lbarreiro/USERS/bridget/references/maps/impute2-1000GP_Phase3/1000GP_Phase3.sample \
--output-log ${log_align}

# shapeit2 phasing run
/project/lbarreiro/SHARED/PROGRAMS/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit \
--input-vcf ${vcf} \
-M /project/lbarreiro/USERS/bridget/references/maps/impute2-1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt \
--input-ref /project/lbarreiro/USERS/bridget/references/maps/impute2-1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz /project/lbarreiro/USERS/bridget/references/maps/impute2-1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz /project/lbarreiro/USERS/bridget/references/maps/impute2-1000GP_Phase3/1000GP_Phase3.sample \
--exclude-snp ${snp_exclude_list} \
--thread 10 \
-O ${temp_dir}/${sample}.chr${chr}.shapeit2.1KGref \
--output-log ${log_main}

# shapeit2 convert haps to vcf
/project/lbarreiro/SHARED/PROGRAMS/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit -convert \
--input-haps ${phased_haps} ${phased_sample} \
--output-vcf ${phased_vcf} \
--output-log ${log_convert}

# zip vcf
bgzip ${phased_vcf}

# Reheader the VCF file
bcftools reheader -s ${output_dir}/${sample}.txt -o ${temp_dir}/${sample}.chr${chr}.shapeit2.1KGref.newID.vcf.gz ${temp_dir}/${sample}.chr${chr}.shapeit2.1KGref.vcf.gz

# Move all relevant files into the sample subdirectory
mv ${temp_dir}/* ${output_dir}/

# Remove the temporary directory
rm -rf ${temp_dir}
```

**Step 3) Generate `multihetsep` input files for MSMC2 using a pair of populations (HG-AG representative pair per-region).**

`gen_multihetsep_per_chr.sh ${chr}` → `KWI4_SMB4.chr${chr}.multihetsep.txt`

```bash
#!/bin/bash
#SBATCH --partition=caslake
#SBATCH --account=pi-\lbarreiro
#SBATCH --time=3:00:00
#SBATCH --mem=5G
#SBATCH --nodes=1

INDIR=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/msmc2/inputs

SAMPLE1=KWI001
SAMPLE2=KWI002
SAMPLE3=KWI004
SAMPLE4=KWI005
SAMPLE5=SMB005
SAMPLE6=SMB010
SAMPLE7=SMB019
SAMPLE8=SMB026

OUTDIR=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/msmc2/outputs/multihetsep/KWI4_SMB4
MASKDIR=/project/lbarreiro/USERS/bridget/programs/msmc-tools/ref_masks
OUTPRE=KWI4_SMB4

# this script combines input from multiple individuals into one input file (multihetsep) per chromosome
# you are using the phased VCFs and masks per sample, plus the chromosome mappability mask, which gives all regions on the chromosome, on which short sequencing reads can be uniquely mapped

/project/lbarreiro/USERS/bridget/programs/msmc-tools/generate_multihetsep.py --chr ${chr} \
--mask=$INDIR/$SAMPLE1/$SAMPLE1.mask.chr${chr}.bed.gz --mask=$INDIR/$SAMPLE2/$SAMPLE2.mask.chr${chr}.bed.gz --mask=$INDIR/$SAMPLE3/$SAMPLE3.mask.chr${chr}.bed.gz --mask=$INDIR/$SAMPLE4/$SAMPLE4.mask.chr${chr}.bed.gz --mask=$INDIR/$SAMPLE5/$SAMPLE5.mask.chr${chr}.bed.gz --mask=$INDIR/$SAMPLE6/$SAMPLE6.mask.chr${chr}.bed.gz --mask=$INDIR/$SAMPLE7/$SAMPLE7.mask.chr${chr}.bed.gz --mask=$INDIR/$SAMPLE8/$SAMPLE8.mask.chr${chr}.bed.gz --mask=$MASKDIR/hs37d5_chr${chr}.mask.bed \
$INDIR/$SAMPLE1/$SAMPLE1.chr${chr}.shapeit2.1KGref.newID.vcf.gz $INDIR/$SAMPLE2/$SAMPLE2.chr${chr}.shapeit2.1KGref.newID.vcf.gz $INDIR/$SAMPLE3/$SAMPLE3.chr${chr}.shapeit2.1KGref.newID.vcf.gz $INDIR/$SAMPLE4/$SAMPLE4.chr${chr}.shapeit2.1KGref.newID.vcf.gz $INDIR/$SAMPLE5/$SAMPLE5.chr${chr}.shapeit2.1KGref.newID.vcf.gz $INDIR/$SAMPLE6/$SAMPLE6.chr${chr}.shapeit2.1KGref.newID.vcf.gz $INDIR/$SAMPLE7/$SAMPLE7.chr${chr}.shapeit2.1KGref.newID.vcf.gz $INDIR/$SAMPLE8/$SAMPLE8.chr${chr}.shapeit2.1KGref.newID.vcf.gz \
> $OUTDIR/$OUTPRE.chr${chr}.multihetsep.txt
```

**Step 4) Run MSMC2 with all chromosomes.**

`msmc_KWI4_SMB4.sh` → `KWI.4ind.out ` and `SMB.4ind.out`

```bash
#!/bin/bash

input_prefix=KWI4_SMB4
multihetsep_dir=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/msmc2/outputs/multihetsep/${input_prefix}
output_dir=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/msmc2/outputs/coal

cd ${output_dir}
source activate msmc2

# for 4 Korowai individuals, all possible pairs of haplotypes
msmc2_Linux -t 12 -s -I 0,1,2,3,4,5,6,7 -o KWI.4ind.out ${multihetsep_dir}/${input_prefix}.chr*.multihetsep.txt

# for 4 Sumba, all possible pairs of haplotypes
msmc2_Linux -t 12 -s -I 8,9,10,11,12,13,14,15 -o SMB.4ind.out ${multihetsep_dir}/${input_prefix}.chr*.multihetsep.txt
```

**Step 5) Create bootstrapped datasets. Columns 1 (`rep`) and 2 (`bootstrap_dir`) of `bootstrap_paths.txt` will used to submit the script for the following step.**

`bootstrap_KWI4_SMB4.sh` → `bootstrap_paths.txt`

```bash
#!/bin/bash

# Run bootstrapping
multihetsep_bootstrap.py --nr_chromosomes 22 bootstrap_dir KWI4_SMB4.chr*.multihetsep.txt

# Write path of bootstrap directories
`ls -d $PWD/bootstrap_dir* | sort -V > bootstrap_paths.txt`

`awk -F'_' '{print $NF " " $0}' bootstrap_paths.txt > temp.txt && mv temp.txt bootstrap_paths.txt`
```

**Step 6) Bootstrap for all KWI and SMB haplotypes separately.**

`msmc2_KWI_4ind_bootstrap.sh ${rep} ${bootstrap_dir}` → `KWI4_SMB4.chr${chr}.multihetsep.txt`

```bash
#!/bin/bash

cd ${bootstrap_dir}
source activate msmc2

# KWI bootstrap
msmc2_Linux -t 6 -s -I 0,1,2,3,4,5,6,7 -o KWI.4ind.bootstrap${rep}.allchr.out bootstrap_multihetsep.chr*.txt

# SMB bootstrap
msmc2_Linux -t 6 -s -I 8,9,10,11,12,13,14,15 -o SMB.4ind.bootstrap${rep}.allchr.out bootstrap_multihetsep.chr*.txt
```

**Step 7) Plot MSMC2 results to create Figure 2a and Supplementary Figure 4.**

`ne_for_manuscript.R`

```r
library(ggplot2)
library(dplyr)
library(scales)
library(stringr)
library(ggh4x)
library(cowplot)

metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
mu <- 1.25e-8 # mutation rate
gen <- 28 # generation time

##############################
####     Data cleaning    ####
##############################

# List all *final.txt files
files <- list.files("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/msmc2/outputs/coal", pattern = "*.final.txt", full.names = TRUE)

# Initialize empty dataframe
alldfs <- data.frame()

# Loop over each file
for (f in files) {
  df <- read.csv(f, sep = "")

  # Scale time and lambda to real units
  df$left_time_boundary <- df$left_time_boundary / mu * gen
  df$right_time_boundary <- df$right_time_boundary / mu * gen
  df$lambda <- (1 / df$lambda) / (2 * mu)

  # Add dataset name from filename
  dataset_name <- str_extract(basename(f), "^[^\\.]+")  # extracts e.g., "BAK"
  df$pop <- dataset_name

  alldfs <- rbind(alldfs, df)
}

# Remove infinite values
alldfs <- alldfs %>% filter_all(all_vars(!is.infinite(.)))

# Remove zero boundaries
alldfs <- alldfs %>% filter(left_time_boundary > 0)

# Merge with metadata
alldfs <- alldfs %>%
  left_join(metadata, by = c("pop" = "pop_code"))

############################################
####     Integrate metadata and plot    ####
############################################

# do not plot recent time points
alldfs <- alldfs %>%
  filter(!is.infinite(lambda) & left_time_boundary > 5000 & right_time_boundary > 5000)

# reorder regions for faceting
region_order <- c("East Asia", "Southeast Asia", "South Asia", "Oceania", "East Africa",
                  "Southern Africa", "East Central Africa", "West Central Africa")
alldfs <- alldfs %>%
  mutate(region = factor(region, levels = region_order))

subsistence_colors <- c(
  "Hunter-gatherer" = "#C21F84",
  "Agriculturalist" = "#FDAE61",
  "Pastoralist" = "#40E0D0")


p <- ggplot(alldfs, aes(x = left_time_boundary, y = lambda, color = subsistence, group = pop)) +
  geom_step(size = 0.5, alpha = 0.7, key_glyph = draw_key_rect) +
  scale_x_log10(labels = scales::scientific, expand = c(0, 0), limits = c(5e3, 1e7)) +
  scale_y_continuous(
    limits = c(0, 35000),
    breaks = seq(0, 35000, by = 5000),  # adjust tick spacing as desired
    labels = scales::comma_format()
  ) +
  scale_color_manual(values = subsistence_colors) +
  facet_wrap(~ region, scales = "free_y", ncol = 4) +
  labs(
    x = "Years before present",
    y = expression(paste("Effective population size ", (N[e]))),
    color = "Subsistence"
  ) +
  theme_classic(base_size = 7) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
  )


#############################################
####     Plot only recent ~100k years    ####
#############################################
recentdfs <- alldfs %>%
  filter(
    !is.infinite(lambda),
    left_time_boundary >= 0,
    left_time_boundary <= 100000,
    right_time_boundary >= 0,
    right_time_boundary <= 100000
  )

p <- ggplot(alldfs, aes(x = left_time_boundary, y = lambda, color = subsistence, group = pop)) +
  geom_step(size = 0.5, alpha = 0.7, key_glyph = draw_key_rect) +
  scale_x_log10(labels = scales::scientific, expand = c(0, 0), limits = c(5e3, 1e7)) +
  scale_y_continuous(
    limits = c(0, 35000),
    breaks = seq(0, 35000, by = 5000),  # adjust tick spacing as desired
    labels = scales::comma_format()
  ) +
  scale_color_manual(values = subsistence_colors) +
  facet_wrap(~ region, ncol = 4, axes = "all_x") +
  labs(
    x = "Years before present",
    y = expression(paste("Effective population size ", (N[e]))),
    color = "Subsistence"
  ) +
  theme_classic(base_size = 7) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
  )

agri_lines <- data.frame(
  region = c("East Asia","Southeast Asia","South Asia","Oceania","East Africa","Southern Africa","East Central Africa","West Central Africa"),
  xintercept = c(9000,4500,9000,10000,7000,3000,5000,5000)
)

q <- ggplot(recentdfs, aes(x = left_time_boundary, y = lambda, color = subsistence, group = pop)) +
  geom_vline(
    data = agri_lines,
    aes(xintercept = xintercept),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "grey40",
           alpha = 0.6
  ) +
  geom_step(size = 0.4, alpha = 0.8) +
  scale_x_log10(
    limits = c(1000, 100000),
    labels = scales::scientific,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(min(recentdfs$lambda, na.rm = TRUE),
               max(recentdfs$lambda, na.rm = TRUE)),
    labels = scales::comma_format()
  ) +
  scale_color_manual(values = subsistence_colors) +
  facet_wrap(~region, scales = "free_y", ncol = 4) +
  labs(
    x = "Years before present",
    y = "Effective population size (Ne)"
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
  )

ggsave("Supp_msmc2_zoomed.pdf", plot = q, width = 8, height = 3.5, units = "in")
save.image("ne_for_manuscript.RData")          
```

## Historical genetic isolation

**Step 1) Use PLINK and the Level 2 VCF to calculate runs of homozygosity for each individual.**

`plink_roh.sh` → `freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_0.05_miss_increasemaxmissing.hom.indiv`

```bash
#!/bin/bash

plink --make-bed --vcf freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.vcf.gz --maf 0.01 --geno 0.05 --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_0.05_miss

plink --bfile freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_0.05_miss  \
  --homozyg --homozyg-snp 50 --homozyg-kb 1000 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 \
  --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_0.05_miss_increasemaxmissing
```

**Step 2) Use ROH outputs to create Figure 2b and Supplementary Figure 5.**

`roh_figure.R` 

```r
###################
####  READ ME  ####
###################
# visualizing ROH for manuscript fig 2C. I'm calculating ROH using the Level 2 VCF
# + 1% MAF filter + 5% missingness filter + relax the number of missing sites permitted per ROH window

library(dplyr)
library(ggplot2)
library(readr)
library(cowplot)
library(ggbeeswarm)

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/roh")


#############################
####    DATA CLEANING    ####
#############################

roh <- read.csv("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/roh/freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf_0.05_miss_increasemaxmissing.hom.indiv", sep="")
metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v3.txt")

# Merge with metadata for population and region info
roh_merged <- roh %>%
  left_join(metadata, by = c("IID" = "study_ID")) %>%
  mutate(ROH_MB = KB / 1000)  # convert to MB

# order by region
region_order <- c("East Asia", "Southeast Asia", "South Asia", "Oceania", "East Africa",
                  "Southern Africa", "East Central Africa", "West Central Africa")


metadata$region <- factor(metadata$region, levels = region_order, ordered = TRUE)
metadata <- metadata[order(metadata$region), ]
pop_order <- unique(metadata$pop_code)
roh_merged <- roh_merged %>%
  mutate(pop_code = factor(pop_code, levels = pop_order))


#################################
####   MANUSCRIPT 2B PLOT    ####
#################################

pops_to_plot <- c("ATY", "PAI", "VAN", "SLI", "AET", "PHP", "ADV", "SHL", "CBU", "BCH", "TWA", "KIG", "SAN", "NAM", "BAK", "FAN")

roh_merged_subset <- roh_merged %>%
  filter(pop_code %in% pops_to_plot) %>%
  mutate(pop_code = factor(pop_code, levels = pops_to_plot))

region_order <- c(
  "East Asia", "Southeast Asia", "South Asia", "Oceania",
  "East Africa", "Southern Africa", "East Central Africa", "West Central Africa"
)

roh_merged_subset <- roh_merged_subset %>%
  mutate(region = factor(region, levels = region_order))

# Plot Fig 2b
plot2b_v2 <- ggplot(roh_merged_subset, aes(x = pop_code, y = ROH_MB, fill = subsistence)) +
  geom_boxplot(
    outlier.shape = NA,
    lwd = 0.2,
    width = 0.5,
    alpha = 0.7
  ) +
  geom_beeswarm(
    aes(color = subsistence),
    size = 0.5,
    priority = "density"
  ) +
  geom_signif(data = annot_df, aes(xmin = xmin, xmax = xmax, y_position = y_position, annotations = annotations), manual = T, inherit.aes = F) +
  facet_wrap(~ region, scales = "free_x", ncol = 4) +
  scale_y_continuous(labels = scales::label_number(accuracy = 10), limits = c(0, 500)) +
  scale_fill_manual(values = subsistence_colors) +
  scale_color_manual(values = subsistence_colors) +
  labs(
    x = "Population",
    y = "Total ROH Length (Mbp)"
  ) +
  theme_classic(base_size = 7) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
  )

ggsave("Fig_2B_ROH_total_length_only.pdf",plot = plot2b_v2,width = 8, height = 4)


# Compute p-value to measure significance in subsistence-based ROH difference
# 1. Recode to binary: HG vs nonHG
roh_merged_subset <- roh_merged_subset %>%
  mutate(subsistence_group = ifelse(subsistence == "Hunter-gatherer", "HG", "nonHG"))

# 2. Run Wilcoxon tests by region
region_pvals <- roh_merged_subset %>%
  group_by(region) %>%
  summarise(
    n_HG = sum(subsistence_group == "HG"),
    n_nonHG = sum(subsistence_group == "nonHG"),
    p_value = if (n_HG > 1 & n_nonHG > 1)
      wilcox.test(ROH_MB ~ subsistence_group)$p.value
    else
      NA_real_,
    .groups = "drop"
  )

# Do this for roh_merged (all pops not just the subset)
roh_all <- roh_merged %>%
  mutate(subsistence_group = ifelse(subsistence == "Hunter-gatherer", "HG", "nonHG"))

# 1) Helper: Wilcoxon per region with safe ties handling + summaries
wilcox_by_region <- function(df) {
  n_HG    <- sum(df$subsistence_group == "HG")
  n_nonHG <- sum(df$subsistence_group == "nonHG")

  med_HG    <- median(df$ROH_MB[df$subsistence_group == "HG"], na.rm = TRUE)
  med_nonHG <- median(df$ROH_MB[df$subsistence_group == "nonHG"], na.rm = TRUE)

  if (n_HG > 1 && n_nonHG > 1) {
    wt <- wilcox.test(ROH_MB ~ subsistence_group, data = df,
                      exact = FALSE, correct = TRUE)  # ties-safe
    data.frame(
      n_HG = n_HG,
      n_nonHG = n_nonHG,
      median_HG = med_HG,
      median_nonHG = med_nonHG,
      median_diff = med_HG - med_nonHG,
      W = unname(wt$statistic),
      p_value = wt$p.value
    )
  } else {
    data.frame(
      n_HG = n_HG, n_nonHG = n_nonHG,
      median_HG = med_HG, median_nonHG = med_nonHG,
      median_diff = med_HG - med_nonHG,
      W = NA_real_, p_value = NA_real_
    )
  }
}

# 2) Run per region
region_tests_all <- roh_all %>%
  group_by(region) %>%
  group_modify(~ wilcox_by_region(.x)) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    stars = case_when(
      !is.na(p_value) & p_value < 0.0005 ~ "***",
      !is.na(p_value) & p_value < 0.005  ~ "**",
      !is.na(p_value) & p_value < 0.05   ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(region)

region_tests_all
readr::write_csv(region_tests_all, "supp_table_ROH_Wilcoxon_by_region.csv")
write.table(
  region_tests_all,
  file = "supp_table_ROH_Wilcoxon_by_region.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Plot all pops for supplementary figure 5
plot_supp_roh <- ggplot(roh_merged, aes(x = pop_code, y = ROH_MB, fill = subsistence)) +
  geom_boxplot(
    outlier.shape = NA,
    lwd = 0.2,
    width = 0.5,
    alpha = 0.3
  ) +
  ggbeeswarm::geom_beeswarm(
    aes(color = subsistence),
    size = 0.1,
    priority = "density"
  ) +
  facet_wrap(~ region, scales = "free_x", ncol = 4) +
  scale_y_continuous(labels = scales::label_number(accuracy = 10)) +
  scale_fill_manual(values = subsistence_colors) +
  scale_color_manual(values = subsistence_colors) +
  labs(
    y = "Total ROH Length (MB)"
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.4),
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
  )

ggsave("Supp_Fig_ROH_total_length_only_all_pops.pdf", plot = plot_supp_roh, width = 11, height = 7)


#############################
####    FIG 2 LEGENDS    ####
#############################
library(patchwork)

subsistence_colors <- c(
  "Hunter-gatherer" = "#C11C84",
  "Agriculturalist" = "#FDAE61",
  "Pastoralist" = "#40E0D0"
)

# Region dummy dataframe
region_order <- c(
  "East Asia", "Southeast Asia", "South Asia", "Oceania",
  "East Africa", "Southern Africa", "East Central Africa", "West Central Africa"
)

region_legend_df <- data.frame(
  region = factor(names(region_colors), levels = names(region_colors)),
  x = 1:length(region_colors),
  y = 1
)
region_legend_df <- data.frame(
  region = factor(region_order, levels = region_order),
  x = 1:length(region_order),
  y = 1
)

# Subsistence dummy dataframe
subsistence_legend_df <- data.frame(
  subsistence = factor(names(subsistence_colors), levels = names(subsistence_colors)),
  x = 1:length(subsistence_colors),
  y = 1
)

p_region_legend <- ggplot(region_legend_df, aes(x = factor(x), y = 1, fill = region)) +
  geom_tile(color = "black", width = 0.9, height = 0.9) +
  scale_fill_manual(values = region_colors, breaks = region_order) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.7, "cm")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave("Fig_2_region_legend.pdf",plot = p_region_legend, width = 6, height = 5)

p_subsistence_legend <- ggplot(subsistence_legend_df, aes(x = factor(x), y = 1, fill = subsistence)) +
  geom_tile(color = "black", width = 0.9, height = 0.9) +
  scale_fill_manual(values = subsistence_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.7, "cm")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("Fig_2_subsistence_legend.pdf",plot = p_subsistence_legend, width = 6, height = 5)


save.image("roh_figure.RData")
```
