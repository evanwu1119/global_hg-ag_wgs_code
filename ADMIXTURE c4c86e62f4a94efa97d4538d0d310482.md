# ADMIXTURE

**Note: Do this with baseline QC filters, 10% site missingness, 1% maf, LD-pruned, and no heterozygosity outliers and related individuals**

---

**Start with PLINK fileset:** `freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf`

**Step 1. LD pruning (windows with 50 snps, sliding 10 snps, r2=0.3)**

```bash
plink --bfile freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf --indep-pairwise 50 10 0.3 --out freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf

Total genotyping rate is 0.982132.
14190641 variants and 918 people pass filters and QC.
Note: No phenotypes present.
Pruned 833593 variants from chromosome 1, leaving 271052.
Pruned 891882 variants from chromosome 2, leaving 280413.
Pruned 774804 variants from chromosome 3, leaving 245400.
Pruned 784341 variants from chromosome 4, leaving 236979.
Pruned 697506 variants from chromosome 5, leaving 221427.
Pruned 706050 variants from chromosome 6, leaving 215663.
Pruned 628043 variants from chromosome 7, leaving 198657.
Pruned 599424 variants from chromosome 8, leaving 187553.
Pruned 463204 variants from chromosome 9, leaving 152169.
Pruned 527555 variants from chromosome 10, leaving 168004.
Pruned 519350 variants from chromosome 11, leaving 163079.
Pruned 504551 variants from chromosome 12, leaving 162628.
Pruned 383739 variants from chromosome 13, leaving 121050.
Pruned 352947 variants from chromosome 14, leaving 115525.
Pruned 311566 variants from chromosome 15, leaving 107279.
Pruned 336122 variants from chromosome 16, leaving 115946.
Pruned 360286 variants from chromosome 17, leaving 109046.
Pruned 304065 variants from chromosome 18, leaving 101310.
Pruned 238454 variants from chromosome 19, leaving 82527.
Pruned 237470 variants from chromosome 20, leaving 80139.
Pruned 152437 variants from chromosome 21, leaving 51828.
Pruned 143377 variants from chromosome 22, leaving 52201.
Pruning complete.  10750766 of 14190641 variants removed.
Marker lists written to
freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf.prune.in
and
freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf.prune.out

plink --bfile freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf --extract freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf.prune.in --make-bed --out freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned
```

PLINK fileset Directory:

`/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned`

- freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned.bed

**Step 2. Run ADMIXTURE**

```bash
#!/bin/bash
#SBATCH --partition=caslake
#SBATCH --account=pi-\lbarreiro
#SBATCH --time=36:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1

INPUTDIR=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned
prefix=freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned

K=$1 # Get the value of K from the first command line argument
r=$2  # Get the value of r from the second command line argument

/project/lbarreiro/SHARED/PROGRAMS/admixture_linux-1.3.0/admixture -s ${RANDOM} --cv ${INPUTDIR}/${prefix}.bed $K | tee log.${prefix}.K${K}.r${r}.out
mv ${prefix}.${K}.Q ${prefix}.K${K}r${r}.Q
```

```bash
sbatch admixture.sh 2 1
sbatch admixture.sh 2 2
sbatch admixture.sh 2 3
sbatch admixture.sh 2 4
sbatch admixture.sh 2 5
```

Extract cross-validation error rates from log.*.out files

`grep -h CV log*.out | awk '{gsub(/[^0-9]+/, "", $3); print $3, $4}' > CV_error.txt`

## Visualizing results

---

**Step 1. Assemble pong input files**

To run pong requires setting up a few files: (1) an ind2pop file that maps individuals to populations; (2) a Qfilemap file that points pong towards which “.Q” files to display

(1) ind2pop file

```bash
cut -d' ' -f1 /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned/freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned.fam > freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned.ind2pop
```

(2) Qfilemap

```bash
prefix=freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned

createQmap(){
local r=$1
local K=$2
awk -v K=$K -v r=$r -v file=${prefix}.K${K}r${r} 'BEGIN{ \
printf("K%dr%d\t%d\t%s.Q\n",K,r,K,file)
}' >> ${prefix}.multiplerun.Qfilemap
}

export -f createQmap

for K in {2..12}; do for r in {1..5}; do createQmap $r $K; done; done
```

**Step 2. Plotting on local computer within pong_python2 conda environment**

```bash
conda activate pong_python2

# all pops
pong -m freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned.multiplerun.Qfilemap -i freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned.ind2pop
```

## Figures for all-pops runs

---

[freeze2_allpops_ADMIXTURE.docx](freeze2_allpops_ADMIXTURE.docx)

*Cross-validation error*

![Untitled](Untitled%201.png)

*K=9 (lowest CVE)*

![Untitled](Untitled%202.png)

```bash
#####################################
#### Set dir and import packages ####
#####################################

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/admixture")
library(ggplot2)
library(tidyr)
library(dplyr)

#############################################
#### Import cross-validation error rates ####
#############################################

# ALL POPS >> K=9 is the lowest
allpops_CV <- read.table("CV_error.txt", quote="\"", comment.char="")
names <- c("K", "cross_validation_error")
colnames(allpops_CV) <- names

ggplot(allpops_CV, aes(x = K, y = cross_validation_error)) + geom_point() + ylab("Cross validation error") + ggtitle("All pops run") + scale_x_continuous(breaks=seq(2, 12, 1))

############################################
#### Plot global ADMIXTURE at lowest CV ####
############################################

# ALL POPS
q <- read.table("freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned.K9r1.Q") 
k <- 9
m <- read.table("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned/freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf_pruned.fam")
colnames(m) <- c("famid", "sampleid", "pat", "mat", "sex", "phenotype")

# combine famid and sampleid with the admixture q matrix
l_df <- cbind(q, m[,c(1,2)])
colnames(l_df) <- c(seq(1:k), 'pop', 'id')
l_df <- l_df %>% drop_na(pop)
gath_L_df <- l_df %>% gather(K, value, -id, -pop) # gather/transform df
gath_L_df <- gath_L_df %>% filter(pop != "ORR")

# order populations using metadata region and country
metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v2.txt")
metadata <- merge(x = metadata, y = l_df, by.x = "study_ID", by.y = "id") # only retaining samples in l_df
metadata <- metadata %>%
  mutate(region = ifelse(pop_code %in% c("ADV", "SHL"), "South Asia", region))

region_order <- c("South Asia", "East Asia", "Southeast Asia", "Oceania", "East Africa", 
                  "Southern Africa", "East Central Africa", "West Central Africa")

# convert the region column to a factor with the custom order — you wont see the results just yet
metadata$region <- factor(metadata$region, levels = region_order, ordered = TRUE)

# sort the entire metadata frame now based on the custom order of the region column, then extract unique pop_code
metadata <- metadata[order(metadata$region), ]
pop_order <- unique(metadata$pop_code)

# define however many colours there are for K value
pop_color <- c("firebrick", "lightblue", "wheat3","coral2", "yellow", "orange", "dark green", "red", "violet", "black", "light green")

global_init <- ggplot(data=gath_L_df, aes(x=reorder(id, value, function(x){max(x)}), y=value,
                                          fill=reorder(K, sort(as.integer(K)))))
global_p <- global_init +
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

global_p

save.image("freeze2_admixture.RData")

```

## Asian + Oceanian-only ADMIXTURE

---

![pong_asia.png](pong_asia.png)

Average CV error per-K (raw: `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/admixture/CV_error_asia.txt`)

| 1 | 0.1677 |
| --- | --- |
| 2 | 0.1555 |
| 3 | 0.1540 |
| 4 | 0.1529 |
| 5 | 0.1522 |
| 6 | 0.1528 |
| 7 | 0.1531 |
| 8 | 0.1536 |
| 9 | 0.1544 |
| 10 | 0.1551 |
| 11 | 0.1563 |
| 12 | 0.1571 |

K = 5

![freeze2_admixture_asia.png](freeze2_admixture_asia.png)

## African-only ADMIXTURE

---

![pong_africa.png](pong_africa.png)

Average CV error per-K (raw: `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/admixture/CV_error_africa.txt`)

| 1 | 0.2745 |
| --- | --- |
| 2 | 0.2593 |
| 3 | 0.2522 |
| 4 | 0.2501 |
| 5 | 0.2515 |
| 6 | 0.2529 |
| 7 | 0.2544 |
| 8 | 0.2556 |
| 9 | 0.2594 |
| 10 | 0.2628 |
| 11 | 0.2664 |
| 12 | 0.2700 |

K=4

![freeze2_admixture_africa.png](freeze2_admixture_africa.png)