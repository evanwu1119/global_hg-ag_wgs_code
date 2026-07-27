# 2) Ancestry and demographic analysis 

This page details the code needed to generate PCA, ADMIXTURE, and MSMC2 results from the Level 3 VCF, corresponding to Figure 1, Supplemental Figures 2, . 


## Global PCAs

**Step 1) Generate PLINK fileset using VCF.** `baseline_filters_site_0.1_cohort_missing_nohet_nokin`

```bash
plink --vcf freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.vcf.gz --make-bed --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin

```

**Step 2) Update PLINK famID to population labels**

```r
setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/updated/baseline_filters_site_0.1_cohort_missing_nohet_nokin")

ref <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v3.txt")
fam <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.fam")
fam$V1 <- ref$pop_code[match(fam$V1, ref$study_ID)]
fam$V6 <- 1

write.table(fam, file = "freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.fam", col.names =F, row.names=F, quote=F)
```

**Step 3) Add MAF filter (1%) and run PCA**

```bash
plink --bfile freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin --maf 0.01 --make-bed --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf

plink --bfile freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf --pca --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf
```

**Step 4) Plot PCA with `freeze2_pcas_updated.R` in `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/pca/`**

## ADMIXTURE

**Step 1) Starting from PLINK outputs generated above for PCA, LD pruning SNPs and make BED file for all samples. Additionally, subset samples by continent.**

`plink_make_bed.sh` -> `freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned_${CONT}`

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

`admixture.sh ${K} ${r} ${CONT}` -> `log.freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned_${CONT}.K${K}.r${r}.out`

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

`create_ind2pop.sh` -> `CV_error_${CONT}.txt` and ` ${CONT}.ind2pop`

```bash
#!/bin/bash

for CONT in global, africa, asia, oceania; do
  grep -h CV ${CONT}/**/log*.out | awk '{gsub(/[^0-9]+/, "", $3); print $3, $4}' > CV_error_${CONT}.txt`
  cut -d' ' -f1 /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/updated/baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned/freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned_${CONT}.fam > ${CONT}.ind2pop
done 
```

**Step 4) Assemble the Qfilemap file that points pong towards which “.Q” files to display.**

`create_Qmap.sh` -> `${CONT}.multiplerun.Qfilemap`

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

```R
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
gath_L_df <- l_df %>% gather(K, value, -id, -pop) # gather/transform df
gath_L_df <- gath_L_df %>% filter(pop != "ORR")

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

asia_init <- ggplot(data=gath_L_df, aes(x=reorder(id, value, function(x){max(x)}), y=value,
                                          fill=reorder(K, sort(as.integer(K)))))
asia_p <- asia_init +
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

asia_p

# OCEANIA (K=2)

q <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_oceania_0.01_maf_pruned.K2r2.Q")
k <- 2
m <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_oceania_0.01_maf_pruned.fam")
colnames(m) <- c("famid", "sampleid", "pat", "mat", "sex", "phenotype")

# combine famid and sampleid with the admixture q matrix
l_df <- cbind(q, m[,c(1,2)])
colnames(l_df) <- c(seq(1:k), 'pop', 'id')
l_df <- l_df %>% drop_na(pop)
gath_L_df <- l_df %>% gather(K, value, -id, -pop) # gather/transform df

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

oceania_init <- ggplot(data=gath_L_df, aes(x=reorder(id, value, function(x){max(x)}), y=value,
                                        fill=reorder(K, sort(as.integer(K)))))
oceania_p <- oceania_init +
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

oceania_p

# AFRICA (K=4)
q <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_africa_0.01_maf_pruned.K4r2.Q")
k <- 4
m <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_africa_0.01_maf_pruned.fam")
colnames(m) <- c("famid", "sampleid", "pat", "mat", "sex", "phenotype")

# combine famid and sampleid with the admixture q matrix
l_df <- cbind(q, m[,c(1,2)])
colnames(l_df) <- c(seq(1:k), 'pop', 'id')
l_df <- l_df %>% drop_na(pop)
gath_L_df <- l_df %>% gather(K, value, -id, -pop) # gather/transform df

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

africa_init <- ggplot(data=gath_L_df, aes(x=reorder(id, value, function(x){max(x)}), y=value,
                                        fill=reorder(K, sort(as.integer(K)))))
africa_p <- africa_init +
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

africa_p

save.image("regional_admixture_plots.RData")
```

## MSMC2

