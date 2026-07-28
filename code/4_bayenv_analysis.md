# Bayenv analysis

These steps contain the code to run Bayenv genome-wide across all 34 populations with 5 replicate runs per SNP. For the purposes of the environmental correlation, agriculturalist and pastoralist populations were coded as 0, while hunter-gatherers were coded as 1. We also reran this pipeline removing pastoralists and removing samples with uncertain allele frequencies (< 15 individuals), this only differs by the removal of 2 or 8 populations i.e. changing the input files/environmental variables, so we don't duplicate the code for those runs here (label = sample filtering: all pops = `bayenv`, no pastoralists = `no_pst`, confident AFs = `nsamples15` + snps filtering: `_maf5_2pops`). Based on QC from our first test runs we found that Bayenv significant SNPs tend to be common especially compared to the global SNP AF distribution. In light of this, we decided to filter for SNPs with a MAF > 0.05 and that were present in at least 2 populations, leaving us with 7.1M SNPs to test. This analysis comprises the results presented in Figure 5, Supplementary Figure 8, and Supplementary Tables 6 and 7.

## Prepare Bayenv input data

**Step 1) Get allele count information per-population starting from the Level 2 VCF.** 

`plink_allelecounts.sh` → `baseline_filters_site_0.1_cohort_missing_nohet_nokin_perpop_keepREFALT.strat`

```bash
#!/bin/bash

# reformat PLINK files' A1 A2 alleles according to the REF + ALT allele
"By default, the minor allele is coded A1 and the major allele is coded A2 (this is used in many output files, e.g. from --freq or --assoc). "

plink --vcf freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.vcf.gz --keep-allele-order --make-bed --out baseline_filters_site_0.1_cohort_missing_nohet_nokin_keepREFALT

# give the PLINK .bim files an ID based on chrXposXXXX
awk 'BEGIN{OFS="\t"} {$2="chr"$1"pos"$4; print}' baseline_filters_site_0.1_cohort_missing_nohet_nokin_keepREFALT.bim > temp_file.bim && mv temp_file.bim baseline_filters_site_0.1_cohort_missing_nohet_nokin_keepREFALT.bim

# to get alternate-allele frequency , you must use the option --keep-allele-order, it will generate the counts and frequency of the A1 allele (otherwise it does that for the minor allele and places the minor allele as A1 regardless of whether it is REF/ALT)
plink --bfile baseline_filters_site_0.1_cohort_missing_nohet_nokin_keepREFALT --freq --family --keep-allele-order --out baseline_filters_site_0.1_cohort_missing_nohet_nokin_perpop_keepREFALT
```

**Step 2) Parse the ALT allele count data.**

`clean_and_split_frq_stats.sh` → `modified_file.txt`

```bash
#!/bin/bash

input_file="baseline_filters_site_0.1_cohort_missing_nohet_nokin_perpop_keepREFALT.frq.strat"

# Step 1: Modify the 2nd column (SNP -> POS) and strip the characters before "pos"
awk 'BEGIN {OFS="\t"} NR==1 {print $1, "POS", $3, $4, $5, $6, $7, $8} NR>1 {gsub(/^.*pos/, "", $2); print}' "$input_file" > modified_file.txt
```


**Step 3) Convert filtered file of allele counts into Bayenv’s SNPSFILE format (tab-separated allele counts for reference and alternate allele, one line each).**

`python_mac2bayenv.sh`  → `bayenv_snpsfile.txt` and `bayenv_snps.txt` 

```python
#!/project/lbarreiro/USERS/evanwu/software/conda_envs/vcf/bin/python
import sys
import os
from itertools import islice
import pandas as pd

NAME = 'bayenv'
DIR = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/data/'
IN = DIR + 'modified_file.txt'
OUT = DIR + NAME '_snpsfile.txt'
LOCI = DIR + NAME + '_snps.txt'

ln = 0 # track SNP number

with open(IN, 'r') as in_file, open(LOCI, 'a') as sites, open(OUT, 'a') as out_file:
    # Read in loci information for 34 pops
    mac = in_file.readline() # skip header col

    while True:
        mac = pd.Series(islice(in_file, 34)).str.split()
        if mac.empty:
            break

				# Increment SNP number and write out SNP chrpos
				ln+=1
        sites.write(ln + '\t' + mac[0][0] + '_' + mac[0][1] + '\n')

				# Write out the SNPSFILE format 
        aac = mac.str.get(6)
        nac = mac.str.get(7)
        rac = nac.astype(int) - aac.astype(int)
        out_file.write(rac.astype(str).str.cat(sep = '\t') + '\n')
        out_file.write(aac.str.cat(sep = '\t') + '\n')
```

**Step 4) Calculate the AF stats for all SNPs.**

`python_calc_aaf.py` → `bayenv_snp_aaf.txt`

```python
#!/project/lbarreiro/USERS/evanwu/software/conda_envs/vcf/bin/python
import sys
import os
import numpy as np

# Calculate stats for all loci
NAME = 'bayenv'
DIR = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/data/'
IN = DIR + NAME + '_snpsfile.txt'
REF = DIR + NAME + '_snps.txt'
AAF = DIR + NAME + '_snp_aaf.txt'

# Given line number n to select from BF info file, get lines 2n-2 and 2n-1 from snps file as python uses n-1 for nth line
# Calculate AAF for each pop, HG/AG, and overall,
# Column order is MAC, MAF, # pops showing MA, HG, AG, pop-specific AAFs
hg_idx = np.array([1,2,3,4,6,8,14,17,25,33]) - 1
ag_idx = np.array([5,7,9,10,12,13,15,16,19,20,21,22,23,24,25,27,28,29,30,31,32,34]) - 1

with open(IN, 'r') as in_file, open(REF, 'r') as sites, open(AAF, 'a') as aaf_file:
    next_site = sites.readline().split()[0]
    next_site = int(next_site)

    for l,line in enumerate(in_file):
        # Just print the SNP counts
        if l == (next_site * 2 - 2):

            ref_counts = np.array(list(map(int, line.split())))

        elif l == (next_site * 2 - 1):

            # Calculate AAF
            alt_counts = np.array(list(map(int, line.split())))
            pop_counts = ref_counts + alt_counts

            aaf = np.divide(alt_counts, pop_counts, out = np.zeros(alt_counts.shape, dtype = float), where = pop_counts!=0)
            hg_aaf = alt_counts[hg_idx].sum() / pop_counts[hg_idx].sum()
            ag_aaf = alt_counts[ag_idx].sum() / pop_counts[ag_idx].sum()
            maf = min(sum(alt_counts) / sum(pop_counts), sum(ref_counts) / sum(pop_counts))

            # Append total AC and number of pops present
            mpop = min(sum(ref_counts != 0), sum(alt_counts != 0))
            mac = min(sum(ref_counts), sum(alt_counts))

            snp_aaf = np.insert(aaf, 0, [mac, maf, mpop, hg_aaf, ag_aaf])
            np.savetxt(aaf_file, snp_aaf[np.newaxis], fmt = '%.4f', delimiter = '\t')

            next_site = sites.readline()

            if len(next_site) == 0:
                break

            next_site = next_site.split()[0]
            next_site = int(next_site)
```

**Step 5) Filter SNPs that are in less than 2 populations or global MAF < 0.05 in `bayenv_snps_aaf.txt` by command line, then use these sites to filter the `snpsfile`.**

`python_filter_snps.py` → `bayenv_maf5_2pops_snpsfile.txt`

```python
#!/project/lbarreiro/USERS/evanwu/software/conda_envs/vcf/bin/python
import sys
import os
import numpy as np

# Filter the snpsfile
NAME = 'bayenv'
DIR = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/data/'
IN = DIR + NAME + '_snpsfile.txt'

# Command to get MAF > 0.05, at least 2 pops: awk '($4 > 0.05 && $5 >= 2){print $0}' bayenv_snp_aaf.txt > bayenv_maf5_2pops_snp_aaf.txt
REF = DIR + NAME + '_maf5_2pops_snp_aaf.txt'
OUT = DIR + NAME + '_maf5_2pops_snpsfile.txt'

# Given line number n to select from BF info file, get lines 2n-2 and 2n-1 from snps file as python uses n-1 for nth line
# Calculate AAF for each pop, HG/AG, and overall,
# Column order is MAC, MAF, # pops showing MA, HG, AG, pop-specific AAFs

hg_idx = np.array([1,2,3,4,6,8,14,17,25,33]) - 1
ag_idx = np.array([5,7,9,10,11,12,13,15,16,18,19,20,21,22,23,24,26,27,28,29,30,31,32,34]) - 1

with open(IN, 'r') as in_file, open(REF, 'r') as sites, open(OUT, 'a') as out_file:
    next_site = sites.readline().split()[0]
    next_site = int(next_site)

    for l,line in enumerate(in_file):
        # Just print the SNP counts
        if l == (next_site * 2 - 2):
            out_file.write(line)

        elif l == (next_site * 2 - 1):
            out_file.write(line)

            next_site = sites.readline()

            if len(next_site) == 0:
                break

            next_site = next_site.split()[0]
            next_site = int(next_site)
```

**Step 6) Sample ~100K random SNPs from the SNPSFILE to estimate the covariance matrix. Do this three separate times to confirm that the covariance matrix has converged.**

`sample_random_snps.sh` → `bayenv_maf5_2popsp_random_snps_{1, 2, 3}.txt`

```bash
#!/bin/bash

name="bayenv_maf5_2pops"
mkdir -p matrix

# Sample 100K random SNPs three times
for i in {1..3}; do
	awk '
  	BEGIN {srand(); n = 0}
  	NR%2 {if (rand() < 0.002) {print $0; n++; next;}}
  	n==1 {print $0; n--}
	' data/${name}_snpsfile.txt > matrix/${name}_random_snps${i}.txt
done
```

**Step 7) Estimate the covariance matrix with 100K SNPs and 100K iterations.** 

`bayenv_matrix.sh` → `bayenv_maf5_2pops_matrixfile.txt` 

```bash
#!/bin/bash

BAYENV=/project/lbarreiro/USERS/evanwu/software/bayenv2_public/bayenv2
POP=34
name="bayenv_maf5_2pops"

# Estimate covariance matrix using ~100k random snps
$BAYENV -i ${name}_random_snps${rep}.txt -p ${POP} -k 100000 -r ${RANDOM} > ${name}_matrix${rep}.out

# Get last converged matrix
tail -n $[ ${POP} + 1 ] ${name}_matrix${rep}.out > ${name}_matrix${rep}.txt
```

**Step 8) Confirm the covariance matrices converged and reflect expected population structure.**

`bayenv_analysis.Rmd {check covariance matrix}`

```R
library(data.table)
library(ggplot2)
library(scales)
library(tidyverse)
library(ComplexHeatmap)

popinfo <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
popinfo$region <- factor(popinfo$region, levels = c("West Central Africa", "East Central Africa", "East Africa", "Southern Africa", "South Asia", "East Asia", "Southeast Asia", "Oceania"))
popinfo$subsistence <- factor(popinfo$subsistence, levels = c("Hunter-gatherer", "Agriculturalist", "Pastoralist"))
# popinfo <- filter(popinfo, samples %in% run_samples) # Define run samples for nopst, nsamples15 to subset popinfo

# Correlation matrix should be similar to pairwise FST estimates
name = "bayenv_maf5_2pops" #  nopst_maf5_2pops, nsamples15_maf5_2pops
npops = 34

maf51 <- fread(paste0("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/matrix/" name, "_matrix1.txt"), select = 1:npops) %>% as.matrix()
rownames(maf51) <- popinfo$pop_code
colnames(maf51) <- popinfo$pop_code

maf52 <- fread(paste0("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/matrix/" name, "_matrix2.txt"), select = 1:26) %>% as.matrix()
rownames(maf52) <- popinfo$pop_code
colnames(maf52) <- popinfo$pop_code

maf53 <- fread(paste0("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/matrix/" name, "_matrix3.txt"), select = 1:26) %>% as.matrix()
rownames(maf53) <- popinfo$pop_code
colnames(maf53) <- popinfo$pop_code

color_fun <- circlize::colorRamp2(breaks = c(0.13, 0.28), colors = c("white","red"))
Heatmap(maf51, col = color_fun, name = "Covariance")
Heatmap(maf52, col = color_fun, name = "Covariance")
Heatmap(maf53, col = color_fun, name = "Covariance")
```

**Step 9) Calculate Bayes Factors with 100K MCMC iterations to obtain demographic-neutral correlations with subsistence. Repeat 5 times per SNP and output calculated correlations and logBFs. Run in batches of 10,000 SNPs per split (`bayenv_environfile.txt` prepared manually, 0 = AG/PST and 1 = HG). `${suffix}` refers to the file extension generated by `split`.**

`bayenv_calc_bf.sh ${suffix} ${r}`  → `bfs/r${r}/${suffix}.bf` 

```bash
#!/bin/bash

# Generate batches of 10K SNPs from the filtered SNPSFILE by command line
# split -a 3 -l 20000 bayenv_maf5_2pops_snpsfile.txt splits/split.

SNPFILE="../../splits/split.${suffix}"
ENVFILE="../../data/bayenv_environfile.txt"
MATFILE="../../matrix/bayenv_maf5_2pops_matrix1.txt"
POPNUM=34
ITNUM=100000
ENVNUM=1

c=0

mkdir -p r${r}
cd r${r}

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
                bayenv2 -i "tmp.${suffix}" -e ${ENVFILE} -m ${MATFILE} -k ${ITNUM} -r ${RANDOM} -p ${POPNUM} -n ${ENVNUM} -t -c -o ${suffix}
        fi
done < ${SNPFILE}
```

**Step 10) Finish up runs where a particular split was terminated prematurely due to walltime. `bads_r{1..5}.txt` are lists of unfinished loci for each replicate run).**

`finish_bfs.sh ${r}`  → `tmp${r}/${suffix}.bf`

```bash
#!/bin/bash

# Generate list of unfinished splits for each replicate run
# for i in {1..5}; do
#         echo "Bads $i"
#         wc -l r${i}/*.bf | head -n -1 | grep -v -E "10000|4658" > bads_r${i}.txt
# done

BAD="../bads_r${r}.txt"

ENVFILE="../../data/bayenv_environfile.txt"
MATFILE="../../matrix/bayenv_matrixfile.txt"
POPNUM=34
ITNUM=100000
ENVNUM=1

LAST="bbk"
LASTN=4658

mkdir -p tmp${r}
cd tmp${r}

# Read in each line of the file, for every loci (=2 lines) calculate bayes factor and replace the temporary file
while IFS=' ' read -r n file; do
        suffix=${file##*/}
        suffix=${suffix%%.*} # split identifier
        if [ $suffix = ${LAST} ]; then
                finish="$(( (${LASTN} - $n) * 2 ))" # number of loci for last split
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

**Step 11) Move complete splits into same folders, concatenate, then paste the replicate runs to get our final data with 5 replicates of BFs and Spearman correlations per SNP**

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

**Step 1) Using Bridget’s Annovar runs, get the SNP function and consequence annotations for the SNPs we tested**

`python_outlier_annos.py` → `bayenv_maf5_2pops_annos.txt`

```python
#!/project/lbarreiro/USERS/evanwu/software/conda_envs/vcf/bin/python
import sys
import os
import numpy as np

ANNO_DIR = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing/annotation'
SNPS = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/data/bayenv_maf5_2pops_snp_aaf.txt'
OUT = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/data/bayenv_maf5_2pops_annos.txt'

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
