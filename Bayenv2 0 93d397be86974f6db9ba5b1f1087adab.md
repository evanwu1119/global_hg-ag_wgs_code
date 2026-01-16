# Bayenv2.0

**References:**

[https://bitbucket.org/tguenther/bayenv2_public/src/master/](https://bitbucket.org/tguenther/bayenv2_public/src/master/)
[https://doi.org/10.1073/pnas.0914625107](https://doi.org/10.1073/pnas.0914625107)
[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3761302/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3761302/)

**Summary of Method:**

> “For each SNP and each environmental variable, we contrasted allele frequencies between the two sets of populations using a Bayesian linear model method that controls for the covariance of allele frequencies between populations due to population history and accounts for differences in sample sizes among populations. **The statistic resulting from this method is a Bayes factor (BF), which is a measure of the support for a model in which a SNP allele frequency distribution is dependent on an environmental variable in addition to population structure, relative to a model in which the allele frequency distribution is dependent on population structure alone.** For subsequent analyses, we use a transformed rank statistic based on the location of each SNP in the overall distribution of BFs. Because we rank each SNP relative to SNPs within the same allele frequency range and from the same ascertainment panel, this transformed rank statistic allows us to make comparisons across SNP sets. To conduct analyses for the two types of variables (ecoregion and subsistence) as a whole, we also calculated for each SNP a minimum rank statistic across all of the variables within each category, which results in a summary statistic for ecoregion and subsistence, respectively.”
> 

## Info from Bridget

---

**Step 1. Get allele count information per-population** 

- Get allele-frequency / counts information per-population from the PLINK dataset in `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/baseline_filters_site_0.1_cohort_missing_nohet_nokin`

```bash
# reformat PLINK files' A1 A2 alleles according to the REF + ALT allele
"By default, the minor allele is coded A1 and the major allele is coded A2 (this is used in many output files, e.g. from --freq or --assoc). "

plink --vcf freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.vcf.gz --keep-allele-order --make-bed --out baseline_filters_site_0.1_cohort_missing_nohet_nokin_keepREFALT

# give the PLINK .bim files an ID based on chrXposXXXX
awk 'BEGIN{OFS="\t"} {$2="chr"$1"pos"$4; print}' baseline_filters_site_0.1_cohort_missing_nohet_nokin_keepREFALT.bim > temp_file.bim && mv temp_file.bim baseline_filters_site_0.1_cohort_missing_nohet_nokin_keepREFALT.bim

# to get alternate-allele frequency , you must use the option --keep-allele-order, it will generate the counts and frequency of the A1 allele (otherwise it does that for the minor allele and places the minor allele as A1 regardless of whether it is REF/ALT)
plink --bfile baseline_filters_site_0.1_cohort_missing_nohet_nokin_keepREFALT --freq --family --keep-allele-order --out baseline_filters_site_0.1_cohort_missing_nohet_nokin_perpop_keepREFALT
```

**Note:** the reason why I am choosing to be consistent with the REF/ALT assignments and in getting the counts for the ALT alleles specifically is because I also intend to do a LOF variant count analysis, and the LOF variants is always the ALT allele rather than the REF.

The output looks something like this where A1 allele (T) is the ALT allele, and A2 allele (C) is the REF allele. 

![Screenshot 2024-06-17 at 17.46.41.png](Screenshot_2024-06-17_at_17.46.41.png)

**Step 2. Parse the ALT allele count data that’s stratified by population**

Working directory:

`/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/baseline_filters_site_0.1_cohort_missing_nohet_nokin_keepREFALT`

- This is the file `baseline_filters_site_0.1_cohort_missing_nohet_nokin_perpop_keepREFALT.frq.strat`

Looks like this

```bash
 CHR               SNP     CLST   A1   A2      MAF    MAC  NCHROBS
   1      chr1pos16280      ADV    C    T        0      0       24 
   1      chr1pos16280      AET    C    T        0      0       56 
   1      chr1pos16280      AGT    C    T        0      0       34 
   1      chr1pos16280      ATI    C    T        0      0       36 
   1      chr1pos16280      ATY    C    T        0      0       18 
   1      chr1pos16280      BAK    C    T  0.09091      4       44 
   1      chr1pos16280      BCH    C    T        0      0       70 
   1      chr1pos16280      CBU    C    T        0      0       44 
   1      chr1pos16280      FAN    C    T     0.08      4       50
```

Use `clean_and_split_frq_stats.sh`

```bash
#!/bin/bash
#SBATCH --partition=caslake
#SBATCH --account=pi-\lbarreiro
#SBATCH --time=36:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1

input_file="baseline_filters_site_0.1_cohort_missing_nohet_nokin_perpop_keepREFALT.frq.strat"

# Step 1: Modify the 2nd column (SNP -> POS) and strip the characters before "pos"
awk 'BEGIN {OFS="\t"} NR==1 {print $1, "POS", $3, $4, $5, $6, $7, $8} NR>1 {gsub(/^.*pos/, "", $2); print}' "$input_file" > modified_file.txt

# Step 2: Separate the file based on unique values in the 3rd column (CLST)
awk 'NR>1 {print > $3"_alt_allele_freq.txt"}' modified_file.txt
```

Now you have a modified_file with cleaned up coordinates, and a population-specific alt_allele_freq file.

```bash
# Strip the modified_file down to just the following columns so it takes less time to read into R (note that this doesn't have the number of chromosomes so you can't extrapolate from the MAF to get LoF counts)

CHR	POS	CLST	MAF
1	16280	ADV	0
1	16280	AET	0
1	16280	AGT	0
1	16280	ATI	0
1	16280	ATY	0
1	16280	BAK	0.09091
1	16280	BCH	0
1	16280	CBU	0
1	16280	FAN	0.08

cut -f1,2,3,6 modified_file.txt > allpops_alt_allele_freq.txt
```

## Summary of files

---

- `modified_file.txt` is the comprehensive file for the allele frequencies of each position within each population
- `allpops_alt_allele_freq.txt` is the stripped version that just has chr, pos, pop, and maf info; if you don’t need the allele info and # of chromosomes, this is much quicker to load
- `<XXX>_alt_allele_freq.txt` contains the information from `modified_file.txt` but it’s split into per-population

## Running Bayenv

---

Working directory: `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv` 

Bayenv2.0 manual:

[bayenv2_manual.pdf](bayenv2_manual.pdf)

These steps contain the code to run the Bayenv software genome-wide across all 34 populations with 5 replicate runs per SNP. For the purposes of the environmental correlation agriculturalist and pastoralist populations were coded as 0 while hunter-gatherers were coded as 1. Based on QC from our first test runs we found that Bayenv significant SNPs tend to be common especially compared to the global SNP AF distribution. In light of this, we decided to filter for SNPs with a MAF > 0.05 and that were present in at least 2 populations, leaving us with 7.1M SNPs to test.  

**Step 1) Filter `modified_file.txt` to remove monomorphic and multiallelic sites**

`python_filter_mac.py`  → `modified_file_filtered.txt`

```python
#!/project/lbarreiro/USERS/evanwu/software/conda_envs/vcf/bin/python
import sys
import os

IN = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/baseline_filters_site_0.1_cohort_missing_nohet_nokin_keepREFALT/modified_file.txt'
REF = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/isec/sites.txt'
OUT = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/modified_file_filtered.txt'

with open(IN, 'r') as in_file, open(REF, 'r') as exclude_sites, open(OUT, 'a') as out_file:
    ex_next = exclude_sites.readline().split()

    for line in in_file:
        tmp = line.split()

        if tmp[0] == ex_next[0] and tmp[1] == ex_next[1]:
            if tmp[2] == "VAN":
                ex_next = exclude_sites.readline().split()
                if len(ex_next) == 0:
                    ex_next = ['0', '0']
            continue

        out_file.write(line)
```

**Step 2) Convert filtered file of allele counts into Bayenv’s SNPSFILE format (tab-separated allele counts for reference and alternate allele, one line each)**

`python_mac2bayenv.sh`  → `bayenv_snpsfile.txt`  and `bayenv_loci.txt` 

```python
#!/project/lbarreiro/USERS/evanwu/software/conda_envs/vcf/bin/python
import sys
import os
from itertools import islice
import pandas as pd

IN = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/modified_file_filtered.txt'
OUT = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_snpsfile.txt'
LOCI = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_loci.txt'

with open(IN, 'r') as in_file, open(LOCI, 'a') as sites, open(OUT, 'a') as out_file:
    # read in loci information for 34 pops
    mac = in_file.readline() # skip header col

    while True:
        mac = pd.Series(islice(in_file, 34)).str.split()
        if mac.empty:
            break

        sites.write(mac[0][0] + '_' + mac[0][1] + '\n')

        aac = mac.str.get(6)
        nac = mac.str.get(7)
        rac = nac.astype(int) - aac.astype(int)
        out_file.write(rac.astype(str).str.cat(sep = '\t') + '\n')
        out_file.write(aac.str.cat(sep = '\t') + '\n')
```

**Step 3) Sample ~100K random SNPs from the SNPSFILE to estimate the covariance matrix**

`sample_random_snps.sh`  → `bayenv_random_snps.txt`

```bash
#!/bin/bash

awk '
  BEGIN {srand(); n = 0}
  NR%2 {if (rand() < 0.002) {print $0; n++; next;}}
  n==1 {print $0; n--}
  ' bayenv_snpsfile.txt > bayenv_random_snps.txt
```

**Step 4) Estimate the covariance matrix with 100K SNPs and 100K iterations and confirm similarity to pairwise FST estimates** 

`bayenv_matrix.sh`  → `bayenv_matrixfile.txt` 

```bash
#!/bin/bash
#SBATCH --partition=caslake
#SBATCH --account=pi-lbarreiro
#SBATCH --time=24:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=bayenv_matrix
#SBATCH --chdir=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv

BAYENV=/project/lbarreiro/USERS/evanwu/software/bayenv2_public/bayenv2

# Estimate covariance matrix using ~100k random snps
$BAYENV -i bayenv_random_snps.txt -p 34 -k 100000 -r 2024 > matrix.out
```

`bayenv_analysis.Rmd {check covariance matrix}`

[bayenv_corr_mat.pdf](bayenv_corr_mat.pdf)

```r
library(ComplexHeatmap)

# Correlation matrix should be similar to pairwise FST estimates
covmat <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_matrixfile.txt", select = 1:34) %>% as.matrix()
covmat <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/OLD_bayenv_matrixfile.txt", select = 1:34) %>% as.matrix()
popinfo <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
popinfo$region <- factor(popinfo$region, levels = c("West Central Africa", "East Central Africa", "East Africa", "Southern Africa", "South Asia", "East Asia", "Southeast Asia", "Oceania"))
popinfo$subsistence <- factor(popinfo$subsistence, levels = c("Hunter-gatherer", "Agriculturalist", "Pastoralist"))

colnames(covmat) <- popinfo$pop_code
rownames(covmat) <- popinfo$pop_code
cormat <- 1 - cov2cor(covmat)
clust <- hclust(dist(cormat))

subsistence <- data.frame(subsistence = popinfo$subsistence, row.names=popinfo$pop_code)
subsistence$subsistence <- factor(subsistence$subsistence)
region <- data.frame(region = popinfo$region, row.names=popinfo$pop_code)
region$region <- factor(popinfo$region, levels = c("East Africa", "East Central Africa", "West Central Africa", "Southern Africa", "East Asia", "Southeast Asia", "South Asia", "Oceania"))

ann_col <- list(subsistence = c("Agriculturalist" = "lightgreen",
                                "Hunter-gatherer" = "lightcoral",
                                "Pastoralist" = "lightblue"),
                region = c("East Africa" = "firebrick4",
                           "East Central Africa" = "red3",
                           "West Central Africa" = "darkorange",
                           "Southern Africa" = "gold",
                           "East Asia" = "orchid",
                           "South Asia" = "purple2",
                           "Southeast Asia" = "dodgerblue3",
                           "Oceania" = "turquoise2"))

pdf(paste0(wd, "figures/bayenv_corr_mat.pdf"), height = 10, width = 12)
pheatmap(cormat, cluster_rows=clust, cluster_cols=clust,
         cellwidth=15, cellheight=15,
         annotation_row = subsistence, annotation_names_row = F,
         annotation_col = region, annotation_names_col = F,
         annotation_colors = ann_col,
         heatmap_legend_param=list(title="1 - AF correlation"),
         display_numbers=T, number_format="%.2f", fontsize_number=6)
dev.off()
```

**Step 5) Calculate the AF stats for all SNPs (`bayenv_snp_aaf.txt`). Then filter loci that are singletons, doubletons, population-specific, or global MAF < 0.05 based on the AAF file and select sites (`bayenv_maf5_2pops_snp_aaf.txt`  sites to keep: 57,349,118 → 7,124,658 SNPs)**

`python_filter_snps.py` → `bayenv_maf5_2pops_snpsfile.txt`

```bash
#!/project/lbarreiro/USERS/evanwu/software/conda_envs/vcf/bin/python
import sys
import os
import numpy as np

# Filter the snpsfile
IN = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_snpsfile.txt'
REF = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_maf5_2pops_snp_aaf.txt'
OUT = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_maf5_2pops_snpsfile.txt'

# Calculate stats for all loci
#IN = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_snpsfile.txt'
#REF = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_snps.txt'
#AAF = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_snp_aaf.txt'

# Given line number n to select from BF info file, get lines 2n-2 and 2n-1 from snps file as python uses n-1 for nth line
# Calculate AAF for each pop, HG/AG, and overall,
# Column order is MAC, MAF, # pops showing MA, HG, AG, pop-specific AAFs
hg_idx = np.array([1,2,3,4,6,8,14,17,25,33]) - 1
ag_idx = np.array([5,7,9,10,12,13,15,16,19,20,21,22,23,24,25,27,28,29,30,31,32,34]) - 1

# open(OUT, 'a') as out_file
with open(IN, 'r') as in_file, open(REF, 'r') as sites, open(OUT, 'a') as out_file: #, open(AAF, 'a') as aaf_file:
    next_site = sites.readline().split()[0]
    next_site = int(next_site)

    for l,line in enumerate(in_file):
        # Just print the SNP counts
        if l == (next_site * 2 - 2):
            out_file.write(line)

            #ref_counts = np.array(list(map(int, line.split())))

        elif l == (next_site * 2 - 1):
            out_file.write(line)

            # Calculate AAF
            #alt_counts = np.array(list(map(int, line.split())))
            #pop_counts = ref_counts + alt_counts
            #
            #aaf = np.divide(alt_counts, pop_counts, out = np.zeros(alt_counts.shape, dtype = float), where = pop_counts!=0)
            #hg_aaf = alt_counts[hg_idx].sum() / pop_counts[hg_idx].sum()
            #ag_aaf = alt_counts[ag_idx].sum() / pop_counts[ag_idx].sum()
            #maf = min(sum(alt_counts) / sum(pop_counts), sum(ref_counts) / sum(pop_counts))
            #
            ## Append total AC and number of pops present
            #mpop = min(sum(ref_counts != 0), sum(alt_counts != 0))
            #mac = min(sum(ref_counts), sum(alt_counts))

            #snp_aaf = np.insert(aaf, 0, [mac, maf, mpop, hg_aaf, ag_aaf])
            #np.savetxt(aaf_file, snp_aaf[np.newaxis], fmt = '%.4f', delimiter = '\t')

            next_site = sites.readline()

            if len(next_site) == 0:
                break

            next_site = next_site.split()[0]
            next_site = int(next_site)
```

**Step 6) Calculate Bayes Factors with 100K MCMC iterations to obtain demographic-neutral correlations with subsistence, repeat 5 times per SNP and output calculated correlations; 10,000 SNPs per split (`bayenv_environfile.txt` , 0 = AG and 1 = HG)**

[`submit.sh`](http://submit.sh) submission script

```bash
#!/bin/bash

# split -a 3 -l 14000 bayenv_noSDS_snpsfile.txt splits/split.
# split -a 3 -l 14000 r1_outlier_snpsfile.txt splits/split.
# split -a 3 -l 20000 bayenv_maf5_2pops_snpsfile.txt splits/split.
readarray -d '' splits < <(find splits -name split.* -print0)
#readarray -t splits < splits.txt

for s in "${splits[@]}"; do
        suf="${s##*.}"
        sbatch --export=suffix=$suf bayenv_calc_bf.sh
done
```

`bayenv_calc_bf.sh`  → `bfs/r{1..5}/${suffix}.bf` 

```bash
#!/bin/bash
#SBATCH --partition=caslake
#SBATCH --account=pi-lbarreiro
#SBATCH --time=36:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=bayenv_calc_bf
#SBATCH --chdir=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bfs/r5

SNPFILE="../splits/split.${suffix}"
#SNPFILE="../outlier_snpsfile.txt"
ENVFILE="../bayenv_environfile.txt"
MATFILE="../bayenv_matrixfile.txt"
POPNUM=34
ITNUM=100000
ENVNUM=1

BAYENV=/project/lbarreiro/USERS/evanwu/software/bayenv2_public/bayenv2
c=0

# Read in each line of the file, for every loci (=2 lines) calculate bayes factor and replace the temporary file
while IFS='' read -r line; do
        if [ $c -eq 0 ]; then
                printf '%s\t' $line > "tmp.${suffix}"
                printf '\n' >> "tmp.${suffix}"
                ((c++))
        else
                printf '%s\t' $line >> "tmp.${suffix}"
                printf '\n' >> "tmp.${suffix}"
                ((c--))
                $BAYENV -i "tmp.${suffix}" -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t -c -o ${suffix}
        fi
done < $SNPFILE
```

**Step 7) Script to finish up runs where a particular split was terminated prematurely due to walltime (`bads_r{1..5}.txt` list of unfinished loci)**

`get_bads.sh → bads_r{1..5}.txt`

```bash
for i in {1..5}; do
        echo "Bads $i"
        wc -l r${i}/*.bf | head -n -1 | grep -v -E "10000|4658" > bads_r${i}.txt
done
```

`finish_bfs.sh`  → `bfs/tmp{1..5}/${suffix}.bf`

```bash
#!/bin/bash
#SBATCH --partition=caslake
#SBATCH --account=pi-lbarreiro
#SBATCH --time=6:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=10
#SBATCH --job-name=finish_bfs
#SBATCH --chdir=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/tmp5

BAD="../bads_r5.txt"

ENVFILE="../bayenv_environfile.txt"
MATFILE="../bayenv_matrixfile.txt"
POPNUM=34
ITNUM=100000
ENVNUM=1

BAYENV=/project/lbarreiro/USERS/evanwu/software/bayenv2_public/bayenv2

# Read in each line of the file, for every loci (=2 lines) calculate bayes factor and replace the temporary file
while IFS=' ' read -r n file; do
        suffix=${file##*/}
        suffix=${suffix%%.*} # split identifier
        if [ $suffix = "bbk" ]; then
                finish="$(( (4658 - $n) * 2 ))" # number of loci for last split
        else
                finish="$(( (10000 - $n) * 2 ))" # number of loci to finish
        fi

        # Get the last few snps needed
        SNPFILE="split.${suffix}"
        tail -n $finish "../splits/split.${suffix}" > $SNPFILE
        c=0

        # Run loop to calculate BF
        while IFS='' read line; do
                if [ $c -eq 0 ]; then
                        printf '%s\t' $line > "tmp.${suffix}"
                        printf '\n' >> "tmp.${suffix}"
                        ((c++))
                else
                        printf '%s\t' $line >> "tmp.${suffix}"
                        printf '\n' >> "tmp.${suffix}"
                        ((c--))
                        $BAYENV -i "tmp.${suffix}" -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t -c -o ${suffix}
                fi
        done < $SNPFILE

        # Append to finished BFs
        #cat ${suffix}.bf >> ../bfs/${suffix}.bf
done < $BAD

```

**Step 8) Move complete splits into same folders, concatenate, then paste the replicate runs to get our final data with 5 replicates of BFs and Spearman correlations per SNP**

`cat_bfs.sh` → `bayenv_maf5_2pops_{bfs, spearman_p, pearson_r}.txt`

```bash
for i in {1..5}; do 
	tmp${i}/*.bf | grep -v finish | xargs -I {} mv {} r${i} # move to correct folders
	cat r${i}/*.bf > r${i}_bfs.txt
done

paste r{1..5}_bfs.txt > tmp.txt
awk '{print $2,$6,$10,$14,$18}' tmp.txt > bayenv_maf5_2pops_bfs.txt
awk '{print $3,$7,$11,$15,$19}' tmp.txt > bayenv_maf5_2pops_spearman_p.txt
awk '{print $4,$8,$12,$16,$20}' tmp.txt > bayenv_maf5_2pops_pearson_r.txt
rm tmp.txt
```

# Downstream analysis

---

From the first run, we found that replicates could vary widely with respect to the magnitude of p-values. With this completed run, we also had the computed environmental correlations across replicates. With BFs and effect sizes/correlations across 5 replicates for the 7M set of tested SNPs, we did some exploratory analysis to decide reasonable cutoffs for significance:

- We tried aggregating BFs across replicates using:
    - Average of log10 BFs
    - Median of log10 BFs
    - At least 4 replicates where BF > cutoff
    - At least 5 replicates where BF > cutoff
- We tried filtering based on absolute value of Spearman correlations at levels of [0, 0.05, 0.1, 0.2]

Takeaways from exploratory analysis:

- logBFs vary more than correlations, with correlations across runs at ~0.55 for BFs and ~0.75 for
- For logBFs the median is more stringent suggesting that large outliers may be unbalancing the average, with no correlation cutoff we get 75K significant SNPs from average and 50K from median
- Strangely there is very little correlation ~0.07 between logBFs and correlations, though you would expect that larger effect sizes have higher p-value
- For the enrichments, keeping a correlation threshold of at least 0.05 preserves interesting signals and likely reduces outlier effects (see above), with increasing strictness our enrichments become a bit too sparse

Ultimately we decided to use the following: 

- Median log10 BF across replicates > 1.5, this significance cutoff was used by these studies:
    - https://www.nature.com/articles/s41467-021-27510-2#Sec7
    - https://onlinelibrary.wiley.com/doi/epdf/10.1111/mec.13722
    - Ultimately they cite a textbook with guidelines for interpretation of Bayes factors: [https://doi.org/10.1093/oso/9780198503682.001.0001](https://doi.org/10.1093/oso/9780198503682.001.0001)
- Median Spearman rho > 0.05
- 32,089 SNPs / 7.2M significant after these thresholds

**Step 0) Exploratory analysis of Bayenv stats across the replicate runs**

`bayenv_analysis.Rmd {explore bf stats}`

```r
#library(biomaRt)
library(clusterProfiler)
library(data.table)
library(GenomicRanges)
library(ggplot2)
#library(ggpattern)
library(parallel)
library(plyranges)
library(rtracklayer)
library(tidyverse)

wd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/"

rep_bfs <- fread(paste0(wd, "bayenv_maf5_2pops_bfs.txt"), header = F)
rep_bfs <- log10(rep_bfs)
rep_bfs[2757637] <- 295.1854 # infinite, set to max of that replicate col

rep_rho <- fread(paste0(wd, "bayenv_maf5_2pops_spearman_p.txt"), header = F)
#rep_cor <- fread(paste0(wd, "bayenv_maf5_2pops_pearson_r.txt"), header = F) # Pearson correlation not recommended to use
cutoff <- 1.5

# BF summary statistics
# bf_cors <- cor(rep_bfs)
# print(bf_cors)
#           V1        V2        V3        V4        V5
# V1 1.0000000 0.5342143 0.5708442 0.5535897 0.5708177
# V2 0.5342143 1.0000000 0.5409288 0.5521865 0.5499941
# V3 0.5708442 0.5409288 1.0000000 0.5656822 0.5718234
# V4 0.5535897 0.5521865 0.5656822 1.0000000 0.5786261
# V5 0.5708177 0.5499941 0.5718234 0.5786261 1.0000000

bf_stats <- data.frame(
  avg = rowMeans(rep_bfs),
  med = apply(rep_bfs, 1, median),
  rep = rowSums(rep_bfs > cutoff),
  var = apply(rep_bfs, 1, var)
)

# > summary(bf_stats$avg)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  -1.1576  -0.9092  -0.8212  -0.6867  -0.6730 295.1854
# > summary(bf_stats$med)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  -1.1921  -0.9474  -0.8738  -0.7632  -0.7511 295.1854

bf_stats$sig_avg <- bf_stats$avg > cutoff # 76,262,
bf_stats$sig_med <- bf_stats$med > cutoff # 52,674, 2285 unique vs avg
bf_stats$sig_rep <- bf_stats$rep >= 4 # 4 = 28218, 19 unique vs avg
bf_stats$sig_rep5 <- bf_stats$rep >= 5 # 5 = 11,388

# Check the cumulative distribution of aggregated BFs
bf_avg_ecdf <- ecdf(bf_stats$avg)
bf_med_ecdf <- ecdf(bf_stats$med)

# Correlation between BFs and abs rho/effect size is small ~0.07, lowers when looking at more stringent cutoffs al ~ 0.06, med ~ 0.045, rep4 ~ 0.035

# Variance across replicates is better controlled by median than average
# summary(bf_stats$var_bfs[bf_stats$sig_avg])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#     0.000     7.324    14.559    32.016    31.026 15176.274
# summary(bf_stats$var_bfs[bf_stats$sig_med]) # median slightly lowers variance
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#     0.000     5.981    13.276    31.342    30.345 15176.274
# summary(bf_stats$var_bfs[bf_stats$sig_rep]) # even higher variance for most significant
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#     0.000     7.873    17.911    40.594    39.917 15176.274

# Check distribution of the correlations
# rho_cors <- cor(rep_rho)
# print(rho_cors)
#           V1        V2        V3        V4        V5
# V1 1.0000000 0.7500832 0.7555189 0.7374246 0.7512302
# V2 0.7500832 1.0000000 0.7510880 0.7639200 0.7606848
# V3 0.7555189 0.7510880 1.0000000 0.7545652 0.7645635
# V4 0.7374246 0.7639200 0.7545652 1.0000000 0.7620621
# V5 0.7512302 0.7606848 0.7645635 0.7620621 1.0000000

rho_stats <- data.frame(
  rho_avg = rowMeans(rep_rho),
  rho_med = apply(rep_rho, 1, median),
  rho_var = apply(rep_rho, 1, var)
)

# Average and median don't look too different for the correlations
# > summary(rho_stats$avg)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -0.33964 -0.02596  0.01713  0.01669  0.06016  0.30015
# > summary(rho_stats$med)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -0.34023 -0.02893  0.01685  0.01665  0.06302  0.29747

abs_rho_avg_ecdf <- ecdf(abs(rho_stats$avg))
abs_rho_med_ecdf <- ecdf(abs(rho_stats$med))

ggplot(data.frame(x = c(-0.5, 2)), aes(x = x)) +
  geom_hline(yintercept = c(0.9, 0.95, 0.99), linetype = 2) +
  geom_vline(xintercept = c(1, 1.5), linetype = 2) +
  stat_function(fun = bf_avg_ecdf, geom = "line", color = "red") +
  stat_function(fun = bf_med_ecdf, geom = "line", color = "blue") +
  labs(x = "Log10 Bayes factor", y = "Cumulative proportion of SNPs") +
  annotate("text", x = 1.15, y = c(0.97, 0.96), label = c("Average", "Median"), color = c("red", "blue")) +
  theme_classic()

ggplot(data.frame(x = c(0.1, 0.3)), aes(x = x)) +
  geom_hline(yintercept = c(0.9, 0.95, 0.99), linetype = 2) +
  geom_vline(xintercept = seq(0.1, 0.3, 0.05), linetype = 2) +
  stat_function(fun = abs_rho_avg_ecdf, geom = "line", color = "red") +
  stat_function(fun = abs_rho_med_ecdf, geom = "line", color = "blue") +
  labs(x = "Abs Spearman rho", y = "Cumulative proportion of SNPs") +
  annotate("text", x = 0.15, y = c(0.94, 0.93), label = c("Average", "Median"), color = c("red", "blue")) +
  theme_classic()

# Number of replicates where signs consistent (positive)
reps_consistent <- data.frame(
  bfs = rowSums(rep_bfs > 0),
  cors = rowSums(rep_rho > 0)
)

num_consistent <- table(reps_consistent$bfs)
#       0       1       2       3       4       5
# 5900359  883282  187258   76288   46009   31462

num_consistent <- table(reps_consistent$cors)
#       0       1       2       3       4       5
# 1597115  733545  536088  567578  901743 2788589

bf_stats <- cbind(bf_stats[, c(1:3,5:8)], rho_stats[1:2])
fwrite(bf_stats, paste0(wd, "bf_stats.txt.gz")) # median 0.05 = 32089, 0.1 = 14981, look at how many independent within these

```

**Step 1) Using Bridget’s Annovar runs, get the SNP function and consequence annotations for the SNPs we tested**

`python_outlier_annos.py` → `bayenv_maf5_2pops_annos.txt`

```python
#!/project/lbarreiro/USERS/evanwu/software/conda_envs/vcf/bin/python
import sys
import os
import numpy as np

ANNO_DIR = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing/annotation'
SNPS = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_maf5_2pops_snp_aaf.txt'
OUT = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_maf5_2pops_annos.txt'

#SNPS = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/r1_bfs.txt'
#OUT = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/r1_outlier_annos.txt'

# For each chromosome, open Bridget's anno file and run through list of SNPs to extract
with open(SNPS, 'r') as sites, open (OUT, 'a') as out_file:
    next_site = sites.readline().split()[1].split("_")

    # While there are still entries to extract
    while len(next_site) != 0:
        # Open the current chromosome anno file
        anno_file = ANNO_DIR + "/chr" + next_site[0] + "_anno.hg19_multianno_nogenotypes_filtered.txt"
        with open(anno_file) as annos:
            # Keep parsing the chromosome anno file until we have found all SNPs in that chromosome
            for line in annos:
                pos = line.split()[1]
                if pos == next_site[1]:
                    out_file.write(line)

                    old_chr = next_site[0]

                    # Catch if end of snp file
                    try:
                        next_site = sites.readline().split()[1].split("_")
                    except:
                        print("No more lines")
                        sys.exit(0)

                    if next_site[0] != old_chr:
                        break

```

**Step 2) Do the same for Bridget’s annotations of immune cCREs**  

python_snp_ccres.py → `bayenv_maf5_2pops_ccres.txt`

```bash
#!/project/lbarreiro/USERS/evanwu/software/conda_envs/vcf/bin/python
import sys
import os
import numpy as np

CCRES = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/i-cCREs_GRCh37.tsv'
SNPS = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_maf5_2pops_snp_aaf.txt'
OUT = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/bayenv_maf5_2pops_ccres.txt'

#SNPS = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/r1_outlier_bfs_aaf.txt'
#OUT = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/r1_outlier_ccres.txt'

#SNPS = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/test_bfs.txt'
#OUT = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/test_ccres.txt'

# Stream through Bridget's CCREs file and get annotation if SNP overlaps a regulatory region, else print "."
with open(SNPS, 'r') as sites, open(CCRES, 'r') as ccres, open (OUT, 'a') as out_file:
    next_site = sites.readline()
    next_site = sites.readline().split()
    chrom, pos = next_site[1].split("_")
    next_ccre = ccres.readline()
    next_ccre = ccres.readline().split()

    chrom = "chr" + chrom
    pos = int(pos)

    # While there are still entries to extract
    while len(next_site) > 0:
        # Catch if end of ccre file
        if next_ccre[0] == "0":
            out_file.write(".\n")
            try:
                next_site = sites.readline().split()
                chrom, pos = next_site[1].split("_")
                chrom = "chr" + chrom
                pos = int(pos)
            except:
                print("No more lines")
                sys.exit(0)
        if chrom == next_ccre[0]:
            # If position is less than CCRE start output "." and go to next SNP
            if pos > int(next_ccre[2]):
                try:
                    next_ccre = ccres.readline().split()
                except:
                    print("No more CCREs")
                    sys.exit(0)
            else:
                # If within range output annotation
                if pos >= int(next_ccre[1]):
                    out_file.write(next_ccre[5] + "\n")
                else:
                    out_file.write(".\n")

                # Move to next SNP after writing output, catch EOF
                next_site = sites.readline().split()
                chrom, pos = next_site[1].split("_")
                chrom = "chr" + chrom
                pos = int(pos)
        else:
            if int(chrom[3:]) < int(next_ccre[0][3:]):
                out_file.write(".\n")
                next_site = sites.readline().split()
                chrom, pos = next_site[1].split("_")
                chrom = "chr" + chrom
                pos = int(pos)
            else:
                next_ccre = ccres.readline().split()
                if next_ccre[0] == "chrX":
                    next_ccre = ["0", "0", "0"]
```

RMarkdown for the rest of the analysis: `bayenv_analysis.Rmd`

**Step 3) Plot the top and bottom SNPs based on BF and correlations for different aggregation methods**

Figures for top and bottom SNPs based on median: 

[top15_med_bf_bot_af.pdf](top15_med_bf_bot_af.pdf)

[top15_med_bf_top_af.pdf](top15_med_bf_top_af.pdf)

[top15_med_cor_bot_af.pdf](top15_med_cor_bot_af.pdf)

[top15_med_cor_top_af.pdf](top15_med_cor_top_af.pdf)

```r
popinfo <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
popinfo$region <- factor(popinfo$region, levels = c("West Central Africa", "East Central Africa", "East Africa", "Southern Africa", "South Asia", "East Asia", "Southeast Asia", "Oceania"))
popinfo$subsistence <- factor(popinfo$subsistence, levels = c("Hunter-gatherer", "Agriculturalist", "Pastoralist"))

cutoff <- 1.5 # log BF threshold

# Create dataframe of all tested SNPs and annotations
all_bfs <- fread(paste0(wd, "bayenv_maf5_2pops_snp_aaf.txt"), select = 3:41, header = F)
colnames(all_bfs) <- c("mac", "maf", "npop_ma", "hg_aaf", "ag_aaf", popinfo$pop_code)

all_annos <- fread(paste0(wd, "bayenv_maf5_2pops_annos.txt"), select = c(1,2,4,5,6,7,9,22,24,25))
colnames(all_annos) <- c("chr", "pos", "ref", "alt", "varregion", "gene", "vartype", "clindiag", "clinrevstat", "clinsig")

ccres <- fread(paste0(wd, "bayenv_maf5_2pops_ccres.txt"), header = F)
ccres <- rbind(ccres, data.frame(V1 = rep(".", times = 13))) # script didn't work so add "." for the last 13 SNPs in chr22 not in regulatory regions
colnames(ccres) <- "ccres"

all_bfs <- cbind(all_bfs, bf_stats, all_annos, ccres)
all_bfs <- all_bfs %>% select(chr, pos, ref, alt, avg, med, rep, sig_avg, sig_med, sig_rep4, sig_rep5, rho_avg, rho_med, varregion, gene, vartype, clindiag, clinrevstat, clinsig, ccres, maf, hg_aaf, ag_aaf, popinfo$pop_code)
saveRDS(all_bfs, paste0(wd, "all_bfs_annotated.rds"))

outliers_list <- lapply(colnames(all_bfs)[8:11], function(x) all_bfs[all_bfs[[x]],])
names(outliers_list) <- colnames(all_bfs)[8:11]
saveRDS(outliers_list, paste0(wd, "outliers_list.rds"))

# Function to plot the aafs given a DF of snps
snp_aaf_plot <- function(df, nsnps, name) {
  pdf(paste0("tests/top", nsnps, "_", name, "_af.pdf"), height = nsnps*1.5, width = 10)
  plot <- df %>% filter(subsistence != "Pastoralist") %>%
    ggplot(aes(x = subsistence, y = pop_aaf, fill = subsistence, color = subsistence)) +
    geom_hline(aes(yintercept = hg_aaf), color = "#F8766D", linetype = 2) +
    geom_hline(aes(yintercept = ag_aaf), color = "#00BA38", linetype = 2) +
    geom_boxplot(alpha = 0.5) +
    scale_discrete_manual(aesthetics = c("color", "fill"), values = c("#F8766D", "#00BA38")) +
    geom_point(position = "jitter", size = 1, alpha = 0.5) +
    labs(x = "Region", y = "Alternate allele frequency", fill = "Subsistence", color = "Subsistence") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    facet_grid(chr_pos~region, switch = "x") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "top",
          strip.placement = "bottom",
          strip.text.x = element_text(angle = 0, size = 7),
          strip.text.y = element_text(angle = 0, hjust = 0),
          strip.background = element_blank(),
          panel.spacing.x = unit(0, "cm"),
          axis.line.x = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 5))
  print(plot)
  dev.off()
}

# Wrapper to plot different outlier dfs sorted by BF and correlation/effect size
plot_top_snps <- function(dfs, nsnps) {
  for (n in names(dfs)) {
    sig <- dfs[[`n`]]
    name <- substr(n, 5, 8)
    filter <- ifelse(name == "avg", name, "med")

    top_list <- list(
      bf_top = sig %>% arrange(desc(!!sym(filter))) %>% slice_head(n = nsnps),
      bf_bot = sig %>% arrange(desc(!!sym(filter))) %>% slice_tail(n = nsnps),
      cor_top = sig %>% arrange(desc(abs(rho_med))) %>% slice_head(n = nsnps),
      cor_bot = sig %>% arrange(desc(abs(rho_med))) %>% slice_tail(n = nsnps)
    )

    for (t in 1:length(top_list)) {
      df <- top_list[[t]]
      df$chr_pos <- paste0(df$chr, ":", df$pos, " ", df$ref, ":", df$alt, ", BF = ", round(df[[`filter`]], 3), ", rho = ", round(df$rho_med, 3),
                                ";\n", gsub("_", " ", df$varregion), " ", gsub(";", ",", df$gene))
      df <- pivot_longer(df, cols = popinfo$pop_code, names_to = "pop_code", values_to = "pop_aaf")
      df <- left_join(df, popinfo[,c(1,3,5)], by = "pop_code")
      df$chr_pos <- factor(df$chr_pos, levels = unique(df$chr_pos))
      top_list[[t]] <- df
    }

    lapply(names(top_list), function(x) snp_aaf_plot(top_list[[x]], nsnps, paste0(name, "_", x)))
  }
}

plot_top_snps(outliers_list, 15)

```

**Step 4) Calculate Chi-squared enrichment of annotation categories for Bayenv significant SNPs**

Encode cCRE definitions: [http://screen-beta.wenglab.org/about](http://screen-beta.wenglab.org/about)

Comparison across different aggregation and cutoffs:

[annotation_chisq_enrich_comparison.pdf](annotation_chisq_enrich_comparison.pdf)

```r
# Wrapper to test all chisq enrichments for aggregation method at different effect size cutoffs
snp_chisq_enrich <- function(method, corr_cutoffs) {
  funcs <- c("varregion", "vartype", "ccres")
  enrich <- list()

  for (f in funcs) {
    for (c in corr_cutoffs) {
      test <- all_bfs[[method]] & abs(all_bfs$rho_med) > c
      cont_tbl <- table(all_bfs[[f]], test)
      cont_tbl <- pivot_wider(as.data.frame(cont_tbl), names_from = 2, values_from = 3)
      colnames(cont_tbl) <- c("category", "background", "significant")
      cont_tbl <- column_to_rownames(cont_tbl, var = "category")

      # Simplify counts depending on annotation
      if (f == "varregion") {
        simp_tbl <- rbind(
          colSums(cont_tbl[grep("intronic", rownames(cont_tbl)),]),
          colSums(cont_tbl[grep("exonic", rownames(cont_tbl)),]),
          colSums(cont_tbl[grep("UTR", rownames(cont_tbl)),]),
          colSums(cont_tbl[c(10, 11, grep("intergenic|stream", rownames(cont_tbl))),]) # Assuming "splicing" and "ncRNA_splicing" are intergenic
        )
        rownames(simp_tbl) <- c("intronic", "exonic", "UTR", "intergenic")
      } else if (f == "vartype") {
        simp_tbl <- rbind(
          colSums(cont_tbl[grep("synonymous", rownames(cont_tbl)),]),
          colSums(cont_tbl[grep("nonsynonymous|gain|loss", rownames(cont_tbl)),]), # coding consequence variants
          colSums(cont_tbl[grep("[.]", rownames(cont_tbl)),])
          )
        rownames(simp_tbl) <- c("synonymous SNV", "nonsynonymous SNV", "noncoding")
      } else {
        simp_tbl <- cont_tbl
        simp_tbl[5,] <- simp_tbl[5,] + simp_tbl[9,]
        simp_tbl <- simp_tbl[-9,] # Combine the TF overlaps
      }

      name <- paste0(f, "_", c)
      tryCatch(
        expr = {
          enrich[[`name`]] <- chisq.test(simp_tbl[, c("background", "significant")])
        }
      )
    }
  }
  return(enrich)
}

corr_cutoffs <- c(0, 0.05, 0.1, 0.2)
avg_enrich <- snp_chisq_enrich("sig_avg", corr_cutoffs)
med_enrich <- snp_chisq_enrich("sig_med", corr_cutoffs)
rep4_enrich <- snp_chisq_enrich("sig_rep4", corr_cutoffs)
rep5_enrich <- snp_chisq_enrich("sig_rep5", corr_cutoffs)

all_enrich <- list()
i <- 1
for (method in list(avg_enrich, med_enrich, rep4_enrich, rep5_enrich)) {
  for (func in c("varregion", "vartype", "ccres")) {
    for (cutoff in corr_cutoffs) {
      name <- paste0(func, "_", cutoff)
      all_enrich[[i]] <- data.frame(
        #method = deparse(substitute(method)),
        func = func,
        cutoff = cutoff,
        pval = method[[`name`]]$p.value,
        stdres = method[[`name`]]$stdres[,2])
      i <- i + 1
    }
  }
}
all_enrich <- bind_rows(all_enrich)
all_enrich$method <- rep(c("avg", "med", "rep4", "rep5"), each = 60)
all_enrich$sig <- all_enrich$pval < 0.05
all_enrich <- rownames_to_column(all_enrich, var = "category")
all_enrich$category <- sapply(all_enrich$category, function(x) strsplit(x, "[.]")[[1]][1])
all_enrich$category[all_enrich$category == ""] <- "no cCRE"
all_enrich$cutoff <- factor(all_enrich$cutoff)
all_enrich$func <- case_when(all_enrich$func == "varregion" ~ "Region",
                             all_enrich$func == "vartype" ~ "Consequence",
                             all_enrich$func == "ccres" ~ "cCRE")

library(ggragged)
pdf(paste0("tests/annotation_chisq_enrich_comparison.pdf"), height = 10, width = 10)
all_enrich %>% ggplot(aes(cutoff, stdres, color = method, fill = method, alpha = sig)) +
  geom_point(stat = "identity") +
  geom_line(aes(group = method)) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Correlation cutoff", y = "Standardized residual", fill = "", color = "") +
  facet_ragged_cols(vars(category), vars(func), scales = "free_y") +
  guides(alpha = "none") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0),
        legend.position = "top")
dev.off()

save(avg_enrich, med_enrich, rep4_enrich, rep5_enrich, all_enrich, file = paste0(wd, "tests/annotation_chisq_enrich.RData"))
```

**Step 5) Gene ontology enrichment using overrepresentation analysis of genes overlapped by the significant SNPs** 

Comparison across different aggregation and cutoffs:

[go_enrich_comparison.pdf](go_enrich_comparison.pdf)

```r
hg19_genes <- readRDS("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/gs_data/hg19_genes.rds")
go_kegg_hallmark_term2gene <- readRDS("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/gs_data/go_kegg_hallmark_term2gene.rds")
go_kegg_hallmark_term2name <- readRDS("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/gs_data/go_kegg_hallmark_term2name.rds")
go_child2 <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/gs_data/GO_child2.csv", header = F)
colnames(go_child2) <- c("term_id", "term_name")

all_loci <- all_bfs %>% select(chr, pos) %>% rename(start = pos) %>% mutate(end = start)
all_loci <- makeGRangesFromDataFrame(all_loci, ignore.strand = T)
all_genes <- hg19_genes %>% join_overlap_left(all_loci)
all_genes <- unique(all_genes$gene_name)
saveRDS(all_genes, paste0(wd, "figures/bayenv_background.rds"))

go_enrich <- function(methods, corr_cutoffs, background, term2gene, term2name) {
  enrich <- list()

  for (m in methods) {
    df <- outliers_list[[`m`]]
    for (c in corr_cutoffs) {
      df$test <- abs(df$rho_med) > c
      outlier_loci <- df %>% filter(test) %>% select(chr, pos, test) %>% rename(start = pos) %>% mutate(end = start)
      outlier_loci <- makeGRangesFromDataFrame(outlier_loci, ignore.strand = T, keep.extra.columns = T)
      outlier_genes <- hg19_genes %>% join_overlap_left(outlier_loci)
      outlier_genes <- outlier_genes[!is.na(outlier_genes@elementMetadata$test),]
      sig_genes <- unique(outlier_genes[outlier_genes@elementMetadata$test,]$gene_name)

      name <- paste0(m, "_", c)
      tryCatch(
        expr = {
          sig_gost <- enricher(sig_genes, universe = background, TERM2GENE = term2gene, TERM2NAME = term2name, pvalueCutoff = 2, minGSSize = 1, maxGSSize = 3e4)
          enrich[[`name`]] <- sig_gost
        }
      )
    }
  }
  return(enrich)
}

go_outlier_enrich <- go_enrich(names(outliers_list), corr_cutoffs, all_genes, go_kegg_hallmark_term2gene, go_kegg_hallmark_term2name)

go_outlier_sig <- lapply(go_outlier_enrich, function(x) x@result %>% filter(p.adjust < 0.05 & ID %in% go_child2$term_id))
go_plot <- lapply(names(go_outlier_sig), function(x) data.frame(
  method = rep(unlist(strsplit(x, "_"))[2], nrow(go_outlier_sig[[`x`]])),
  cutoff = rep(unlist(strsplit(x, "_"))[3], nrow(go_outlier_sig[[`x`]])),
  ID = go_outlier_sig[[`x`]]$ID,
  Description = go_outlier_sig[[`x`]]$Description,
  p.adjust = go_outlier_sig[[`x`]]$p.adjust
))
go_plot <- bind_rows(go_plot)
go_plot$logp <- -log10(as.numeric(go_plot$p.adjust))

# Heatmap of results
pdf(paste0("tests/go_enrich_comparison.pdf"), height = 15, width = 10)
go_plot %>% ggplot(aes(cutoff, Description, fill = logp)) +
  geom_tile() +
  scale_fill_viridis_b() +
  facet_wrap(~method, nrow = 1) +
  labs(x = "Correlation cutoff", fill = "-Log10 p value") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 3),
        strip.background = element_blank())
dev.off()

save(go_outlier_sig, go_outlier_enrich, file = paste0(wd, "tests/outlier_GO.RData"))

```

**Step 6) Logistic regression analysis done only for the final significant set of SNPs**

```r
all_bfs <- fread(paste0(wd, "all_bfs_annotated.txt.gz"))
uniq_annots <- c("intergenic", "UTR", "exonic", "intronic",
                 "noncoding", "nonsynonymous", "synonymous",
                 "nonregulatory", "pELS", "dELS", "PLS", "H3K4me3", "CTCF", "TF")

# Prepare data for logistic regression
enrich_data <- data.frame(significant = all_bfs$sig+0,
                        intergenic = ifelse(grepl("intergenic|stream", all_bfs$varregion), 1, 0),
                        UTR = ifelse(grepl("UTR", all_bfs$varregion), 1, 0),
                        exonic = ifelse(grepl("exonic", all_bfs$varregion), 1, 0),
                        intronic = ifelse(grepl("intronic|splicing", all_bfs$varregion) & !grepl("exonic", all_bfs$varregion), 1, 0),
                        noncoding = ifelse(grepl("[.]|unknown", all_bfs$vartype), 1, 0),
                        nonsynonymous = ifelse(grepl("nonsynonymous|start|stop", all_bfs$vartype), 1, 0),
                        synonymous = ifelse(all_bfs$vartype == "synonymous SNV", 1, 0),
                        nonregulatory = ifelse(all_bfs$ccres == ".", 1, 0),
                        pELS = ifelse(all_bfs$ccres == "pELS", 1, 0),
                        dELS = ifelse(all_bfs$ccres == "dELS", 1, 0),
                        PLS = ifelse(all_bfs$ccres == "PLS", 1, 0),
                        H3K4me3 = ifelse(all_bfs$ccres == "CA-H3K4me3", 1, 0),
                        CTCF = ifelse(all_bfs$ccres == "CA-CTCF", 1, 0),
                        TF = ifelse(grepl("TF", all_bfs$ccres), 1, 0))
colnames(enrich_data) <- c("significant", uniq_annots)

# Function to run logistic regression and store results, including confidence intervals
run_logistic_regression <- function(data, annotation_name, element_cols) {
  results <- list()
  for (c in element_cols) {
    # Fit the logistic regression model using the binary high iHS column
    model <- glm(formula = as.formula(paste("significant", "~", c)), data = data, family = binomial)

    # Extract odds ratio and p-value for the element
    odds_ratio <- exp(coef(model)[c])
    p_value <- summary(model)$coefficients[c, "Pr(>|z|)"]

    # Calculate 95% CI for the odds ratio
    conf_int <- exp(confint(model)[c,])
    lower_ci <- conf_int[1]
    upper_ci <- conf_int[2]

    # Append results to the provided results dataframe
    results[[`c`]] <- c(c, odds_ratio, lower_ci, upper_ci, p_value)
  }
  results <- bind_rows(results) %>% t()
  colnames(results) <- c(annotation_name, "OR", "LCI", "UCI", "p")
  return(results)
}

logistic_enrich <- run_logistic_regression(enrich_data, "annotation", uniq_annots) %>% as.data.frame()
logistic_enrich$annotation_category <- rep(c("Variant region", "Variant type", "Immune cCRE"), times = c(4,3,7))
logistic_enrich[,2:5] <- sapply(logistic_enrich[,2:5], as.numeric)

saveRDS(logistic_enrich, file = "figures/logistic_enrich.rds")

```

## **Results summary**

---

Figure components:

- Manhattan plot summarizes 30K SNPs and shows the genes of  25/108  nonsynonymous/missense mutations, tagging immune, metabolic, sperm, and hair related functions
    - We tried using intergenic and genic SNPs but the genes are more difficult to discern, still some immune/metabolic but more developmental/neural stuff
    - Enrichment of the nonsynonymous genes show immune and actin/myosin
- AF plots for a few cool SNPs, should show some global HG-AG diff but also some region-specific signals
- Logistic annotation enrichment shows immune regulatory and genic annotations are enriched among Bayenv significant SNPs
- GO enrichment had >1K total enrichments, we filtered by the second child set to get ~200 and among those display around 50 that are grouped in more interpretable biological processes

Setup

```r
#figures
library(ggfastman)
library(ggrepel)
library(ggnewscale)
library(khroma)
library(cowplot)
```

Manhattan plot

```r
all_bfs <- fread(paste0(wd, "all_bfs_annotated.txt.gz"))

# Manhattan plot of tested BFs
manhattan <- all_bfs %>% dplyr::select(chr, pos, med, sig, ref, alt, gene, vartype, varregion) %>% rename(pvalue = med) %>% mutate(color = ifelse(sig, "#DF536B", NA))
manhattan$chr <- paste0("chr", manhattan$chr)
manhattan$pvalue[manhattan$pvalue > 1.5 & !manhattan$sig] = 1.5

# Annotate genes by GO membership
go_kegg <- left_join(go_kegg_hallmark_term2gene, go_kegg_hallmark_term2name, by = "term")
immune_genes <- go_kegg[grepl("immune", go_kegg$name),]$gene %>% unique()
metabolic_genes <- go_kegg[
  grepl("carbohydrate|lipid|protein|amino acid", go_kegg$name) &
  grepl("metabolic", go_kegg$name),]$gene %>% unique()
sperm_genes <- go_kegg[grepl("sperm", go_kegg$name),]$gene %>% unique()
muscle_genes <- go_kegg[grepl("muscle", go_kegg$name),]$gene %>% unique()
hair_genes <- go_kegg[grepl("hair", go_kegg$name),]$gene %>% unique()

manhattan$annot <- case_when(
  manhattan$gene %in% immune_genes ~ "Immune",
  manhattan$gene %in% sperm_genes ~ "Sperm",
  manhattan$gene %in% metabolic_genes ~ "Metabolic",
  manhattan$gene %in% muscle_genes ~ "Muscle",
  manhattan$gene %in% hair_genes ~ "Hair",
  T ~ "Other"
)
manhattan$annot <- factor(manhattan$annot, levels = c("Other", "Immune", "Metabolic", "Muscle", "Hair", "Sperm"))
saveRDS(manhattan, paste0(wd, "data/manhattan.rds"))
manhattan <- readRDS(paste0(wd, "data/manhattan.rds"))

nonsyn_manhattan_plot <- fast_manhattan(manhattan, build = "hg19", color1 = "dodgerblue", color2 = "lightskyblue", log10p = F, highlight = "color", alpha = 0.7, pointsize = 2, pixels = c(1024,1024), dodge_x = F) +
  geom_hline(yintercept = 1.5, linetype = 2, color = "black") +
  geom_text_repel(data = . %>%
                    filter(sig & grepl("nonsynonymous|gain|loss", vartype)) %>%
                    arrange(desc(pvalue)) %>%
                    distinct(gene, .keep_all = T) %>%
                    slice(c(1:24,46)),
                  aes(label = gene), color = "black", size = 5/.pt,
                  segment.size = 0.25, segment.linetype = 1, min.segment.length = 0.3,
                  ylim = c(2,NA),
                  force = 2) +
  new_scale("shape") +
  geom_point(data = . %>%
               filter(sig & grepl("nonsynonymous|gain|loss", vartype)) %>%
               arrange(desc(pvalue)) %>%
               distinct(gene, .keep_all = T) %>%
               slice(c(1:24,46)),
             aes(shape = annot), size = 1, color = "black") +
  #scale_x_continuous(expand = c(0.1,0), breaks = c(1:22), labels = as.character(c(1:22)), name = "Chromosome") +
  scale_y_continuous(limits = c(-1.5, 100), expand = c(0.01, 0), name = expression(paste("log"[10], " Bayes factor"))) +
  guides(shape = guide_legend(title = NULL, nrow = 1, position = "top"), color = "none") +
  theme_classic(base_size = 7)

pdf(paste0(wd, "figures/bayenv_top25_nonsyn_manhattan.pdf"), height = 4, width = 6)
#print(top_manhattan_plot)
print(nonsyn_manhattan_plot)
dev.off()

```

AF plots

```r
# AF of chosen SNPs
all_top <- bayenv_outlier %>%
  filter(pos %in% c(32602137, 4692233, 33057351) |
           (vartype == "nonsynonymous SNV" & grepl("HLA-B|SPATA|MAST|COL20A1", gene)))
all_top$gene[c(2,7)] <- c("CDYL", "HLA-DQA1") # change to closest gene
all_top <- all_top %>%  mutate(chr_pos = paste0(chr, ":", pos, " ", varregion, #" ",
                                                ifelse(ccres == ".", " ", paste0(" ", ccres, " ")),
                                                gsub(";", ",", gene)))
all_top <- all_top[c(1,4,5,7),] %>% arrange(chr, pos)

# AF of top nonsynonymous SNPs
all_top <- bayenv_outlier %>%
  filter(sig & !grepl("intergenic|stream|ncRNA", varregion)) %>%
           #vartype %in% c("nonsynonymous SNV", "stopgain", "stoploss", "startgain")) %>%
  arrange(desc(med)) %>%
  distinct(gene, .keep_all = T) %>%
  #slice_head(n = 20) %>%
  mutate(chr_pos = paste0(chr, ":", pos, " ", gsub(";", ",", gene)))
all_top <- all_top[c(1:24,46),] # top 25 + FUT6
top25$annot <- case_when(
  top25$gene %in% immune_genes ~ "Immune",
  top25$gene %in% sperm_genes ~ "Sperm",
  top25$gene %in% metabolic_genes ~ "Metabolic",
  top25$gene %in% muscle_genes ~ "Muscle",
  top25$gene %in% hair_genes ~ "Hair",
  top25$gene %in% olfactory_genes ~ "Olfactory",
  T ~ "Other"
)
top25$annot <- factor(top25$annot, levels = c("Other", "Immune", "Metabolic", "Muscle", "Hair", "Sperm"))
all_top <- pivot_longer(all_top, cols = popinfo$pop_code, names_to = "pop_code", values_to = "pop_aaf")
all_top <- left_join(all_top, popinfo[,c(1,3,5)], by = "pop_code")
all_top$chr_pos <- factor(all_top$chr_pos, levels = unique(all_top$chr_pos))
all_top <- all_top %>% group_by(chr_pos, region) %>% mutate(group_mean = mean(pop_aaf), adj_aaf = scale(pop_aaf, scale = F))
all_top$region_abbr <- case_when(
  all_top$region == "West Central Africa" ~ "WC Africa",
  all_top$region == "East Central Africa" ~ "EC Africa",
  all_top$region == "East Africa" ~ "E Africa",
  all_top$region == "Southern Africa" ~ "S Africa",
  all_top$region == "South Asia" ~ "S Asia",
  all_top$region == "Southeast Asia" ~ "SE Asia")
all_top$region_abbr <- factor(all_top$region_abbr, levels = c("WC Africa", "EC Africa", "E Africa", "S Africa", "S Asia", "SE Asia"))

pdf(paste0(wd, "figures/top25_genic_nonsyn_af.pdf"), height = 10, width = 8)
snp_af_plot <- all_top %>%
  filter(region != "Oceania" & region != "East Asia" & group_mean != 0 & group_mean != 1) %>%
  ggplot(aes(x = region_abbr, y = adj_aaf, fill = subsistence, color = subsistence)) +
    geom_hline(yintercept = 0, color = "black", linetype = 2) +
    geom_boxplot(alpha = 0.5, outliers = F, linewidth = 0.5, fatten = 1) +
    scale_discrete_manual(aesthetics = c("color", "fill"), values = c(
      "Hunter-gatherer" = "#FFB6C1",
      "Agriculturalist" = "#FDAE61",
      "Pastoralist" = "#40E0D0")) +
    labs(x = "", y = "Regional AAF deviation", fill = "", color = "") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    facet_wrap(~chr_pos, scales = "free_y", ncol = 1) +
    theme_classic(base_size = 7) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "top",
          strip.placement = "bottom",
          strip.text.x = element_text(size = 5),
          strip.background = element_blank(),
          #panel.spacing.x = unit(0, "cm"),
          axis.line.x = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1))
print(snp_af_plot)
dev.off()

```

Logistic regression plot

```r
muted <- color("muted")

# Plot exp vs. observed?
load(paste0(wd, "tests/annotation_chisq_enrich.RData"))
med_enrich_0.05 <- med_enrich[c(2,6,10)]
chisq_results <- lapply(med_enrich_0.05, function(x) data.frame(
  expected = x$expected[,2],
  observed = x$observed[,2],
  stdres = x$stdres[,2]))
chisq_results <- bind_rows(chisq_results)
chisq_results <- rownames_to_column(chisq_results, var = "annotation")
chisq_results$annotation[8] <- "no cCRE"
chisq_results$annotation <- factor(chisq_results$annotation, levels =
                                     c("intergenic", "UTR", "exonic", "intronic",
                                       "noncoding", "nonsynonymous SNV", "synonymous SNV",
                                       "no cCRE", "CA", "CA-CTCF", "CA-TF", "CA-H3K4me3", "dELS", "pELS", "PLS"))
#chisq_results$annotation <- str_replace(chisq_results$annotation, "^\\w{1}", toupper)
chisq_results$group <- rep(c("SNP region", "SNP consequence", "Immune cCRE overlap"), times = c(4,3,8))
chisq_results$group <- factor(chisq_results$group, levels = rev(c("SNP region", "SNP consequence", "Immune cCRE overlap")))

pdf(paste0(wd, "figures/bayenv_chisq_enrich.pdf"), height = 6, width = 4)
chisq_results %>%
  ggplot(aes(y = annotation, x = stdres, fill = group, color = group)) +
    geom_col(alpha = 0.5, show.legend = F) +
    scale_color_manual(values = muted(8)[1:3]) +
    scale_fill_manual(values = muted(8)[1:3]) +
    labs(x = "Standard residual", y = NULL, color = "", fill = "") +
    facet_wrap(~group, scales = "free_y", ncol = 1) +
    theme_classic() +
    theme(strip.background = element_blank())
dev.off()

# Plot logistic enrichments
logistic_enrich <- readRDS(paste0(wd, "data/logistic_enrich.rds"))
logistic_enrich$annotation[8] <- "no cCRE"
logistic_enrich$annotation <- factor(logistic_enrich$annotation, levels =
                                     c("intergenic", "UTR", "exonic", "intronic",
                                       "noncoding", "nonsynonymous", "synonymous",
                                       "no cCRE", "CTCF", "TF", "H3K4me3", "dELS", "pELS", "PLS"))
logistic_enrich$sig <- factor(ifelse(logistic_enrich$p < 0.05, 1, 0))

pdf(paste0(wd, "figures/bayenv_logistic_enrich.pdf"), height = 5, width = 5)
logistic_enrich_plot <- logistic_enrich %>%
  filter(!(annotation %in% c("CTCF", "TF"))) %>%
  ggplot(aes(x = OR, y = annotation, color = annotation_category, alpha = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2) +
  # Place p-values above the odds ratio points
  #geom_text(aes(label = ifelse(p < 0.05, paste0("p = ", scientific(p, digits = 2), "*"), "")),
            #nudge_y = 0.5, hjust = 0.5, size = 5, size.unit = "pt", color = "black") +  # Position above with nudge_y
  scale_color_manual(values = muted(8)[2:4]) +
  scale_alpha_manual(values = c(0.3, 1)) +
  facet_wrap(~annotation_category, scales = "free_y", ncol = 1) +
  guides(color = "none", alpha = "none") +
  labs(y = "", x = "Odds Ratio", color = "") +
  theme_classic(base_size = 7) +
  theme(legend.position = c(0.8,0.9),
        strip.background = element_blank())
print(logi
```

GO plot

```r
muted <- color("muted")

# Separate enrichments into groups of immune, metabolic, muscle, development/growth, and nervous system as major groups
load(paste0(wd, "tests/outlier_GO.RData"))
go_outlier_sig <- go_outlier_enrich[["sig_med_0.05"]]@result %>% filter(p.adjust < 0.05)
saveRDS(go_outlier_sig, paste0(wd, "data/bayenv_go_enrich.rds"))
go_outlier_sig <- readRDS(paste0(wd, "data/bayenv_go_enrich.rds"))

test <- go_outlier_sig %>% filter(
  #ID %in% go_child2$term_id |
    !grepl("GO:", ID) | grepl("alcohol", Description))
go_outlier_sig$enrich_group <- case_when(
  grepl("immune|leukocyte", go_outlier_sig$Description, ignore.case = T) ~ "Immune",
  grepl("metabolic|carbohydrate|lipid|catabolic", go_outlier_sig$Description, ignore.case = T) ~ "Metabolic",
  grepl("muscle|actin|myosin|Z disc|band|sarcomere|costamere", go_outlier_sig$Description, ignore.case = T) & !(grepl("acting", go_outlier_sig$Description, ignore.case = T)) ~ "Muscle",
  grepl("axon|synap|dendrit|perikaryon|neuron|action potential", go_outlier_sig$Description, ignore.case = T) ~ "Nervous",
  grepl("development|growth|morph", go_outlier_sig$Description, ignore.case = T) ~ "Development",
  T ~ "Other"
)
go_outlier_sig$logp <- -log10(go_outlier_sig$p.adjust)
go_outlier_sig$gene_ratio <- sapply(go_outlier_sig$GeneRatio, function(x) eval(parse(text = x)))
go_outlier_sig$desc <- str_to_sentence(go_outlier_sig$Description)
go_outlier_sig$enrich_group <- factor(go_outlier_sig$enrich_group, levels = rev(c("Immune", "Metabolic", "Muscle", "Nervous", "Development", "Other")))
go_outlier_sig <- go_outlier_sig %>% group_by(enrich_group) %>% arrange(gene_ratio, .by_group = T)
go_outlier_sig$desc <- factor(go_outlier_sig$desc, levels = go_outlier_sig$desc)
go_outlier_sig$enrich_group <- factor(go_outlier_sig$enrich_group, levels = c("Immune", "Metabolic", "Muscle", "Nervous", "Development", "Other"))
saveRDS(go_outlier_sig, paste0(wd, "data/go_outlier_sig.rds"))
readRDS(paste0(wd, "data/go_outlier_sig.rds"))

pdf(paste0(wd, "figures/bayenv_go_enrich.pdf"), height = 6, width = 5)
go_enrich_plot <- go_outlier_sig %>%
  filter(enrich_group != "Other" & enrich_group != "Nervous") %>%
  ggplot(aes(y = desc, x = gene_ratio, color = enrich_group)) +
    geom_segment(aes(x = 0, y = desc, xend = gene_ratio, yend = desc), show.legend = F) +
    geom_point(aes(color = enrich_group, size = logp)) +
    scale_color_manual(values = muted(8)[5:8]) +
    scale_size(breaks = c(10, 25, 40),
               labels = as.character(c(10,25,40)),
               range = c(0.5, 3)) +
    labs(x = "Gene ratio", y = NULL, size = expression(paste("-log"[10], " (p)"))) +
    guides(color = guide_legend(title = NULL, order = 2), size = guide_legend(order = 1)) +
    #facet_wrap(~enrich_group, scales = "free_y", ncol = 1) +
    theme_classic(base_size = 7) +
    theme(legend.position = c(0.7,0.6),
          legend.box = "horizontal",
          legend.spacing = unit(0, "mm"),
          legend.key.spacing = unit(0, "mm"),
          axis.text.y = element_text(size = 5))
print(go_enrich_plot)
dev.off()

```

Arrange panel

```r
save(nonsyn_manhattan_plot, snp_af_plot, logistic_enrich_plot, go_enrich_plot, file = paste0(wd, "figures/bayenv_figures_v2.RData"))
load(paste0(wd, "figures/bayenv_figures.RData"))

bayenv_figures <- plot_grid(
  plot_grid(
    plot_grid(nonsyn_manhattan_plot),
    plot_grid(snp_af_plot, logistic_enrich_plot, labels = c("b", "c"), rel_widths = c(1, 0.8), nrow = 1, label_size = 8),
  rel_heights = c(0.8,1), nrow = 2), go_enrich_plot, rel_widths = c(2, 1), nrow = 1, labels = c("a", "d"), label_size = 8)

```