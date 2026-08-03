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

POP=34
name="bayenv_maf5_2pops"

# Estimate covariance matrix using ~100k random snps
bayenv2 -i ${name}_random_snps${rep}.txt -p ${POP} -k 100000 -r ${RANDOM} > ${name}_matrix${rep}.out

# Get last converged matrix
tail -n $[ ${POP} + 1 ] ${name}_matrix${rep}.out > ${name}_matrix${rep}.txt
```

**Step 8) Confirm the covariance matrices converged and reflect expected population structure. Save plot for Supplementary Figure 8a later.**

`bayenv_analysis.Rmd {check covariance matrix}` → `figs8b_mat.rds`

```r
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

mat_plot <- pheatmap(maf5, color = rev(hcl.colors(10, "Reds 3")),
         cellwidth = 7, cellheight = 7, fontsize = 7, show_row_dend = F,# show_column_dend = F,
         heatmap_legend_param = list(title = "Covariance",
                                     title_gp = gpar(fontsize = 7),
                                     labels_gp = gpar(fontsize = 7)))

saveRDS("mat_plot", "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/results_bayenv/data/figs8a_mat.rds")
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
                        bayenv2 -i "tmp.${suffix}" -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t -c -o ${suffix}
                fi
        done < $SNPFILE

        # Append to finished BFs
        #cat ${suffix}.bf >> ../bfs/${suffix}.bf
done < ${BAD}
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

`python_snp_annos.py` → `bayenv_maf5_2pops_annos.txt`

```python
#!/project/lbarreiro/USERS/evanwu/software/conda_envs/vcf/bin/python
import sys
import os
import numpy as np

NAME = 'bayenv_maf5_2pops'
DIR = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/data/'
ANNO_DIR = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing/annotation'
SNPS = DIR + NAME + '_snp_aaf.txt'
OUT = DIR + NAME + '_annos.txt'

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

**Step 2) Do the same for immune cCRE annotations**  

`python_snp_ccres.py` → `bayenv_maf5_2pops_ccres.txt`

```python
#!/project/lbarreiro/USERS/evanwu/software/conda_envs/vcf/bin/python
import sys
import os
import numpy as np

NAME = 'bayenv_maf5_2pops'
DIR = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/data/'
CCRES = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/i-cCREs_GRCh37.tsv'
SNPS = DIR + NAME +  '_snp_aaf.txt'
OUT = DIR + NAME + '_ccres.txt'

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

**Step 3) Read all Bayenv outputs, SNP information, and SNP annotations into R as a single dataframe. Decide Bayenv SNPs as having median logBF > 1.5 and median Spearman correlation > 0.05 across all 5 replicates. Combine the significant SNPs across the all population, no PST, and n >= 15 runs to create Supplementary Table 6.**

`bayenv_analysis.Rmd {process outliers}` → `results_${name}/{all, outlier}_bfs_annotated.txt` and `bayenv_sup_table6.xlsx`

```r
library(data.table)
library(scales)
library(tidyverse)
library(openxlsx)

# Path and run params
name <- "bayenv"
wd <- paste0("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/results", name, "/")
dd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/data/"

popinfo <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
popinfo$region <- factor(popinfo$region, levels = c("West Central Africa", "East Central Africa", "East Africa", "Southern Africa", "South Asia", "East Asia", "Southeast Asia", "Oceania"))
popinfo$subsistence <- factor(popinfo$subsistence, levels = c("Hunter-gatherer", "Agriculturalist", "Pastoralist"))

cutoff <- 1.5 # log BF threshold
cor_cut <- 0.05 # rho threshold

# Read BFs and correlations from all replicate runs
rep_bfs <- fread(paste0(wd, "bfs/", name, "_maf5_2pops_bfs.txt"), header = F)
rep_bfs <- log10(rep_bfs)

rep_rho <- fread(paste0(wd, "bfs/", name, "_maf5_2pops_spearman_p.txt"), header = F)

bf_stats <- data.frame(
  med = apply(rep_bfs, 1, median),
  rho_med = apply(rep_rho, 1, median)
)

bf_stats$sig <- bf_stats$med > cutoff & abs(bf_stats$rho_med) > cor_cut
bf_stats <- bf_stats %>% dplyr::select(med, rho_med, sig)

# Create dataframe of all tested SNPs and annotations
all_bfs <- fread(paste0(dd, name, "_maf5_2pops_snp_aaf.txt"), select = 3:41, header = F)
colnames(all_bfs) <- c("mac", "maf", "npop_ma", "hg_aaf", "ag_aaf", popinfo$pop_code)

all_annos <- fread(paste0(dd, name, "_maf5_2pops_annos.txt"), select = c(1,2,4,5,6,7,9,22,24,25))
colnames(all_annos) <- c("chr", "pos", "ref", "alt", "varregion", "gene", "vartype", "clindiag", "clinrevstat", "clinsig")

ccres <- fread(paste0(dd, name, "_maf5_2pops_ccres.txt"), header = F)
ccres <- rbind(ccres, data.frame(V1 = rep(".", times = 13))) # script didn't work so add "." for the last 13 SNPs in chr22 not in regulatory regions
colnames(ccres) <- "ccres"

all_bfs <- cbind(all_bfs, bf_stats, all_annos, ccres)
all_bfs <- all_bfs %>% select(chr, pos, ref, alt, avg, med, rep, sig_avg, sig_med, sig_rep4, sig_rep5, rho_avg, rho_med, varregion, gene, vartype, clindiag, clinrevstat, clinsig, ccres, maf, hg_aaf, ag_aaf, popinfo$pop_code)
fwrite(all_bfs, paste0(wd, "data/all_bfs_annotated.txt.gz"))

# Subset to only significant SNPs
outlier_bfs <- all_bfs[all_bfs$sig,]
fwrite(outlier_bfs, paste0(wd, "data/outlier_bfs_annotated.txt.gz"))

# Combine the significant files across runs to get Supplementary Table 6
pre="/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/"
snps1 <- fread(paste0(pre, "results_bayenv/data/outlier_bfs_annotated.csv"))
snps2 <- fread(paste0(pre, "results_nsamples15/data/outlier_bfs_annotated.csv"))
snps3 <- fread(paste0(pre, "results_nopst/data/outlier_bfs_annotated.csv"))

write.xlsx(list("All populations" = snps1, "N >= 15 populations" = snps2, "No PST populations" = snps3), file = "bayenv_sup_table6.xlsx")
```

**Step 4) Perform over-representation analysis of Bayenv significant SNPs in GO, KEGG, and MSigDB gene sets. We tried overlapping SNPs and genes using both gene bodies and 100kb around the gene center, using gene bodies for the manuscript figures. For more information about the gene set resources see the code page `0_annotation_curation`. Combine the significant enrichments across the all population, no PST, and n >= 15 runs to create Supplementary Table 7.** 

`bayenv_analysis.Rmd {functional enrichment}` → `results_${name}/outlier_GO{_100kb}.csv`, `results_${name}/outlier_GO{_100kb}.RData`, and `bayenv_sup_table7.xlsx`

```r
library(data.table)
library(GenomicRanges)
library(plyranges)
library(rtracklayer)
library(scales)
library(tidyverse)
library(clusterProfiler)
library(openxlsx)

# Path and run params
name <- "bayenv"
wd <- paste0("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/results", name, "/")
dd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/data/"
gd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/gs_data/"

# Load gene resources
hg19_genes <- readRDS(paste0(gd, "hg19_genes.rds"))
go_kegg_hallmark_term2gene <- readRDS(paste0(gd, "go_kegg_hallmark_term2gene.rds"))
go_kegg_hallmark_term2name <- readRDS(paste0(gd, "go_kegg_hallmark_term2name.rds"))
go_child2 <- fread(paste0(gd, "GO_child2.csv"), header = F)
colnames(go_child2) <- c("term_id", "term_name")

# Create genomic background as all genes with at least one BF
all_bfs <- fread(paste0(wd, "all_bfs_annotated.txt.gz"))
all_loci <- all_bfs %>% select(chr, pos) %>% rename(start = pos) %>% mutate(end = start)
all_loci <- makeGRangesFromDataFrame(all_loci, ignore.strand = T)
all_genes <- hg19_genes %>% join_overlap_left(all_loci)
all_genes <- unique(all_genes$gene_name)
saveRDS(all_genes, paste0(wd, "data/bayenv_background.rds"))
all_genes <- readRDS(paste0(wd, "data/bayenv_background.rds"))

# Create coordinates of 100kb regions surrounding genes
 <- as.data.frame(hg19_genes)
g100$middle <- g100$start + floor((g100$end - g100$start) / 2)
g100$start <- g100$middle - 50000 + 1
g100$end <- g100$middle + 50000
g100$width <- 100000
g100 <- as_granges(g100, keep_mcols = T)

# Overlap outlier snps with gene bodies
outlier_bfs <- fread(paste0(wd, "data/outlier_bfs_annotated.txt.gz"))
outlier_loci <- outlier_bfs %>% dplyr::select(chr, pos, med, rho_med, sig) %>% rename(start = pos, logbf_med = med) %>% dplyr::mutate(end = start)
outlier_loci <- makeGRangesFromDataFrame(outlier_loci, ignore.strand = T, keep.extra.columns = T)

outlier_genes <- hg19_genes %>% join_overlap_left(outlier_loci)
outlier_genes <- outlier_genes[!is.na(outlier_genes@elementMetadata$sig),]
sig_genes <- unique(outlier_genes[outlier_genes@elementMetadata$sig,]$gene_name) 

sig_gost <- enricher(sig_genes, universe = all_genes, TERM2GENE = go_kegg_hallmark_term2gene, TERM2NAME = go_kegg_hallmark_term2name, pvalueCutoff = 2, minGSSize = 1, maxGSSize = 3e4)

sig_sets <- sig_gost@result %>% filter(p.adjust < 0.05) 

save(sig_genes, sig_gost, sig_sets, file = paste0(wd, "data/outlier_GO.RData"))
fwrite(sig_sets, file = paste0(wd, "data/outlier_GO.csv"))
fwrite(sig_gost@result, file = paste0(wd, "data/outlier_GO_all.csv"))

# Repeat enrichment for 100kb regions
outlier_genes <- g100 %>% join_overlap_left(outlier_loci)
outlier_genes <- outlier_genes[!is.na(outlier_genes@elementMetadata$sig),]
fwrite(as.data.frame(outlier_genes), paste0(wd, "data/bayenv_sig_outlier_genes.txt"), sep = "\t")
sig_genes <- unique(outlier_genes[outlier_genes@elementMetadata$sig,]$gene_name)

sig_gost <- enricher(sig_genes, universe = all_genes, TERM2GENE = go_kegg_hallmark_term2gene, TERM2NAME = go_kegg_hallmark_term2name, pvalueCutoff = 2, minGSSize = 1, maxGSSize = 3e4)

sig_sets <- sig_gost@result %>% filter(p.adjust < 0.05) 

save(sig_genes, sig_gost, sig_sets, file = paste0(wd, "data/outlier_GO_100kb.RData"))
fwrite(sig_sets, file = paste0(wd, "data/outlier_GO_100kb.csv"))
fwrite(sig_gost@result, file = paste0(wd, "data/outlier_GO_100kb_all.csv"))

# Write out all enrichment results across runs as Supplementary Table 7
write.xlsx(list("All populations" = snps1, "N >= 15 populations" = snps2, "No PST populations" = snps3), file = "bayenv_sup_table5.xlsx")

go1 <- fread(paste0(pre, "results_bayenv/data/outlier_GO.csv"))
go2 <- fread(paste0(pre, "results_bayenv/data/outlier_GO_100kb.csv"))
go3 <- fread(paste0(pre, "results_nsamples15/data/outlier_GO.csv"))
go4 <- fread(paste0(pre, "results_nsamples15/data/outlier_GO_100kb.csv"))
go5 <- fread(paste0(pre, "results_nopst/data/outlier_GO.csv"))
go6 <- fread(paste0(pre, "results_nopst/data/outlier_GO_100kb.csv"))

write.xlsx(list("All populations" = go1, "N >= 15 populations" = go2, "No PST populations" = go3, "All populations 100kb" = go4, "N >= 15 populations 100kb" = go5, "No PST populations 100kb" = go6), file = "bayenv_sup_table7.xlsx")
```

**Step 5) Calculate logistic regression enrichment of SNP region, SNP gene consequence, and cCRE type for Bayenv significant SNPs vs background.**

`bayenv_logistic_enrich.R` → `logistic_enrich.rds`

```r
library(data.table)
library(tidyverse)

wd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/results_bayenv/"
all_bfs <- fread(paste0(wd, "data/all_bfs_annotated.txt.gz"))
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
run_logistic_regression <- function(data, element_cols) {
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
  colnames(results) <- c("annotation", "OR", "LCI", "UCI", "p")
  return(results)
}

logistic_enrich <- run_logistic_regression(enrich_data, uniq_annots) %>% as.data.frame()
#logistic_enrich$annotation_category <- rep(c("Variant region", "Variant type", "Immune cCRE"))

saveRDS(logistic_enrich, file = "data/logistic_enrich.rds")
```

**Step 6) Look for enrichment of SNPs putatively under selection from ancient Eurasian DNA in Bayenv significant SNPs from [Akbari et al](https://www.nature.com/articles/s41586-026-10358-1). The annotation of SNPs is in `Selection_Summary_Statistics_01OCT2025.tsv.gz`.**

`bayenv_analysis.Rmd {akbari snps enrich}` → `akbari_enrich_list.rds` and `figs8d_akbari.rds`

```r
library(data.table)
library(ggplot2)
library(scales)
library(tidyverse)

wd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/results_bayenv/"

# Test whether Bayenv SNPs are enriched for ancient DNA selection SNPs
all_bfs <- fread(paste0(wd, "data/all_bfs_annotated.txt.gz"))
colnames(all_bfs)[7] <- "sig1.5"
joins <- c("chr", "pos", "ref", "alt")
sigs <- c("sig0", "sig0.5", "sig1", "sig1.5", "sig2")
all_bfs$rho_cut <- abs(all_bfs$rho_med) > 0.05
all_bfs <- all_bfs %>% mutate(
  sig0 = med > 0 & rho_cut,
  sig0.5 = med > 0.5 & rho_cut,
  sig1 = med > 1 & rho_cut,
  sig2 = med > 2 & rho_cut
)

akbari <- fread("/project/lbarreiro/SHARED/REFERENCES/Selection_Summary_Statistics_01OCT2025.tsv.gz")
akbari$Category <- ifelse(akbari$POSTERIOR >= 0.95, "hc", "lc")
colnames(akbari)[1:4] <- joins

common <- left_join(select(all_bfs, all_of(c(joins, sigs))), select(akbari, all_of(joins), "Category"), by = joins)

akbari_enrich_list <- lapply(sigs, function(x) {
  tbl <- table(select(common, "Category", all_of(x)), useNA = "ifany")
  colnames(tbl) <- c("non-outlier", "outlier")
  rownames(tbl) <- c("high-confidence", "low-confidence", "NA")
  return(tbl)
})
names(akbari_enrich_list) <- sigs
saveRDS(akbari_enrich_list, paste0(wd, "data/akbari_enrich_list.rds"))

# Transform to enrichment prop and plot
enrich_df <- data.frame(
  cut = sub("sig", "", names(akbari_enrich_list)) %>% as.numeric(),
  hvl = sapply(akbari_enrich_list, function(x)
    (x[1,2] * sum(x[2,])) / (sum(x[1,]) * x[2,2])),
  hva = sapply(akbari_enrich_list, function(x)
    (x[1,2] * sum(x[2:3,])) / (sum(x[1,]) * (x[2,2] + x[3,2])))
)
enrich_df <- pivot_longer(enrich_df, cols = c("hvl", "hva"), names_to = "type", values_to = "enrich_ratio")
enrich_df <- enrich_df %>% filter(type == "hvl")

# Get p-values of low confidence
akbari_chisq <- lapply(akbari_enrich_list, function(x) chisq.test(x[1:2,]))
enrich_df$p <- sapply(akbari_chisq, function(x) x$p.value)
enrich_df$plab <- scientific(enrich_df$p)

# Plot enrichment
akbari_plot <- ggplot(enrich_df, aes(x = cut, y = enrich_ratio)) +
  geom_hline(yintercept = 1, linetype = 2, color = "black") +
  geom_line() +
  geom_point() +
  geom_text(aes(label = plab), nudge_y = 0.4, size = 5, size.unit = "pt") +
  labs(x = expression(paste("log"[10], " Bayes factor cutoff")),
       y = "Enrichment ratio", color = "Background type") +
  theme_classic(base_size = 7)

saveRDS(akbari_plot, paste0(wd, "data/figs8d_akbari.rds"))
```

**Step 7) Code not shown, but we ran a test run using all SNPs to see whether common SNPs are more likely to be outliers and plotted that result for panel a. Arrange the plots using data from above to create Supplementary Figure 8.**

`bayenv_analysis.Rmd {arrange sup plots}`

```r
library(data.table)
library(ggplot2)
library(scales)
library(tidyverse)
library(khroma)
library(cowplot)

wd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/results_bayenv/"
muted <- color("muted")

mat_plot <- readRDS(paste0(wd, "data/figs8a_mat.rds"))
akbari_plot <- readRDS(paste0(wd, "data/figs8d_akbari.rds"))

# Generate panel b maf plot
load("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/results_OLD/global_2pops_ecdf.RData")
r1_out <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/results_OLD/r1_outlier_bfs_aaf.txt", select = 6)
r1_out$maf <- ifelse(r1_out$aaf > 0.5, 1-r1_out$aaf, r1_out$aaf)
out_ecdf <- ecdf(r1_out$maf)

maf_plot <- ggplot(data.frame(x = c(0, 0.5)), aes(x = x, color = c("Bayenv significant", "Global"))) +
  geom_abline(slope = 2, intercept = 0, linetype = 2) +
  geom_vline(xintercept = 0.05) +
  stat_function(aes(color = "1"), fun = out_ecdf, geom = "line") +
  stat_function(aes(color = "2"), fun = global_maf_ecdf, geom = "line") +
  labs(x = "MAF", y = "Cumulative proportion of SNPs", color = "") +
  scale_color_manual(values = c("tomato", "deepskyblue"), labels = c("Bayenv outliers", "All SNPs")) +
  theme_classic(base_size = 7) +
  theme(legend.position = "top")

saveRDS(maf_plot, paste0(wd, "data/figs8b_maf.rds")

# Generate panel c logistic plot
logistic_enrich <- readRDS(paste0(wd, "data/logistic_enrich.rds"))
logistic_enrich$annotation[8] <- "no cCRE"
logistic_enrich$annotation <- factor(logistic_enrich$annotation, levels =
                                     c("intergenic", "UTR", "exonic", "intronic",
                                       "noncoding", "nonsynonymous", "synonymous",
                                       "no cCRE", "CTCF", "TF", "H3K4me3", "dELS", "pELS", "PLS"))
logistic_enrich$sig <- factor(ifelse(as.numeric(logistic_enrich$p) < 0.05, 1, 0))
logistic_enrich$annotation_category <- rep(c("Variant region", "Variant type", "Immune cCRE"), times = c(4,3,7))
logistic_enrich[,2:5] <- sapply(logistic_enrich[,2:5], as.numeric)
logistic_enrich$annotation_category <- factor(logistic_enrich$annotation_category, levels = unique(logistic_enrich$annotation_category))

logistic_plot <- logistic_enrich %>%
  filter(!(annotation %in% c("CTCF", "TF", "H3K4me3"))) %>%
  ggplot(aes(x = OR, y = annotation, color = annotation_category, alpha = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2) +
  # Place p-values above the odds ratio points
  #geom_text(aes(label = ifelse(p < 0.05, paste0("p = ", scientific(p, digits = 2), "*"), "")),
            #nudge_y = 0.5, hjust = 0.5, size = 5, size.unit = "pt", color = "black") +  # Position above with nudge_y
  scale_color_manual(values = muted(8)[2:4]) +
  scale_alpha_manual(values = c(0.5, 1)) +
  facet_wrap(~annotation_category, scales = "free_y", ncol = 1) +
  guides(color = "none", alpha = "none") +
  labs(y = "", x = "Odds Ratio", color = "") +
  theme_classic(base_size = 7) +
  theme(legend.position = c(0.8,0.9),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

saveRDS(logistic_plot, paste0(wd, "data/figs8c_logistic.rds")

# Arrange all the plots for Supplementary Figure 8
g <- grid::grid.grabExpr(print(mat_plot))
supp_figures <- plot_grid(
    plot_grid(g, logistic_plot, nrow = 1, labels = c("a", "c"), rel_widths = c(1, 0.5), label_size = 8),
    plot_grid(maf_plot, akbari_plot, labels = c("b", "d"), nrow = 1, label_size = 8),
  rel_heights = c(1, 0.4), nrow = 2)
supp_figures

save_plot(paste0(wd, "figures/supplementary_figure_8.pdf"), supp_figures, base_width = 7, base_height = 7)
```

**Step 8) Plot all tested SNPs with a Manhattan plot of the logBFs.**

`bayenv_analysis.Rmd {manhattan plot}` → `fig5a_manhattan.rds`

```r
library(data.table)
library(ggplot2)
library(scales)
library(tidyverse)
library(ggfastman)
library(ggrepel)
library(ggnewscale)

wd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/results_bayenv/"
gd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/gs_data/"

# Load gene resources
hg19_genes <- readRDS(paste0(gd, "hg19_genes.rds"))
go_kegg_hallmark_term2gene <- readRDS(paste0(gd, "go_kegg_hallmark_term2gene.rds"))
go_kegg_hallmark_term2name <- readRDS(paste0(gd, "go_kegg_hallmark_term2name.rds"))
go_child2 <- fread(paste0(gd, "GO_child2.csv"), header = F)
colnames(go_child2) <- c("term_id", "term_name")

# Manhattan plot of tested BFs
all_bfs <- fread(paste0(wd, "data/all_bfs_annotated.txt.gz"))
manhattan <- all_bfs %>% dplyr::select(chr, pos, med, sig, ref, alt, gene, vartype, varregion) %>% rename(pvalue = med) %>% mutate(color = ifelse(sig, "#DF536B", NA))
manhattan$rows <- 1:nrow(manhattan)
manhattan$chr <- paste0("chr", manhattan$chr)
manhattan$pvalue[manhattan$pvalue > 1.5 & !manhattan$sig] = 1.5

# Annotate genes by GO membership
top25 <- manhattan %>% arrange(desc(pvalue)) %>% filter(sig & !duplicated(gene)) %>% slice_head(n = 30)
top25$gene1 <- sapply(top30$rows, function(x) {
  genes <- strsplit(manhattan$gene[x], ";") %>% unlist()
  if (length(genes) > 1) {
    d1 <- manhattan$pos[x] - hg19_genes[hg19_genes$gene_name == genes[1], 2:3]
    d2 <- manhattan$pos[x] - hg19_genes[hg19_genes$gene_name == genes[2], 2:3]
    if (nrow(d1) == 0) genes[2]
    else if (nrow(d2) == 0) genes[1]
    else ifelse(min(abs(d1)) <= min(abs(d2)), genes[1], genes[2])
  } else {
    manhattan$gene[x]
  }
})
top25 <- top25 %>% filter(!duplicated(gene1)) %>% slice_head(n = 25)
manhattan$gene1 <- NA
manhattan$gene1[top25$rows] <- top25$gene1

go_kegg <- left_join(go_kegg_hallmark_term2gene, go_kegg_hallmark_term2name, by = "term")
top25_go <- go_kegg %>% filter(gene %in% top25$gene1)
top25$annot <- case_when(
  top25$gene1 %in% unique(top25_go$gene[grep("immune", top25_go$name)]) ~ "Immune",
  top25$gene1 %in% unique(top25_go$gene[grep("sperm", top25_go$name)]) ~ "Spermatogenesis",
  top25$gene1 %in% unique(top25_go$gene[grep("metab|lipid|leucine", top25_go$name)]) ~ "Metabolic",
  top25$gene1 %in% unique(top25_go$gene[grep("neur|nerv", top25_go$name)]) ~ "Nervous",
  top25$gene1 %in% unique(top25_go$gene[grep("devel|grow", top25_go$name)]) ~ "Development",
  T ~ "Other"
)

manhattan$annot <- NA
manhattan$annot[top25$rows] <- top25$annot
manhattan$annot <- factor(manhattan$annot, levels = c("Other", "Immune", "Metabolic", "Nervous", "Spermatogenesis"))

manhattan_plot <- fast_manhattan(manhattan, build = "hg19", color1 = "dodgerblue", color2 = "lightskyblue", log10p = F, highlight = "color", alpha = 0.7, pointsize = 3.5, pixels = c(2048,512), dodge_x = T) +
  geom_hline(yintercept = 1.5, linetype = 2, color = "black") +
  geom_text_repel(data = . %>% filter(!is.na(gene1)),
                  aes(label = gene1), color = "black", size = 5/.pt,
                  segment.size = 0.25, segment.linetype = 1, min.segment.length = 0.25,
                  ylim = c(2,NA),
                  force = 0.7) +
  new_scale("shape") +
  geom_point(data = . %>% filter(!is.na(gene1)),
             aes(shape = annot), size = 1, color = "black") +
  scale_y_continuous(limits = c(-1.5, 51), expand = c(0.01, 0), name = expression(paste("log"[10], " Bayes factor"))) +
  guides(shape = guide_legend(title = NULL, nrow = 1, position = "top"), color = "none", x = guide_axis(n = 1)) +
  theme_classic(base_size = 7)

saveRDS(manhattan_plot, paste0(wd, "data/fig5a_manhattan.rds"))
```

**Step 9) Plot the 4 representative significant SNPs with their population AFs, normalized within geographical regions.**

`bayenv_analysis.Rmd {af plot}` → `fig5b_af_plot.rds`

```r
library(data.table)
library(ggplot2)
library(scales)
library(tidyverse)

wd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/results_bayenv/"

# Load in population information
popinfo <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_popinfo.txt")
popinfo$region <- factor(popinfo$region, levels = c("West Central Africa", "East Central Africa", "East Africa", "Southern Africa", "South Asia", "East Asia", "Southeast Asia", "Oceania"))
popinfo$subsistence <- factor(popinfo$subsistence, levels = c("Hunter-gatherer", "Agriculturalist", "Pastoralist"))

# Process the significant SNPs df for plotting
bayenv_outlier <- fread(paste0(wd, "data/outlier_bfs_annotated.txt.gz"))
bayenv_outlier <- pivot_longer(bayenv_outlier, cols = popinfo$pop_code, names_to = "pop_code", values_to = "pop_aaf")
bayenv_outlier <- left_join(bayenv_outlier, popinfo[,c(1,3,5)], by = "pop_code")
bayenv_outlier <- bayenv_outlier %>%
  filter(region != "Oceania" & region != "East Asia") %>% # drop cuz no hg/ag
  mutate(chr_pos = paste0(chr, ":", pos, " ", varregion, " ", ifelse(ccres == ".", "", ccres),
                          " ", gsub(";", ",", gene), ", ", "logBF = ", round(med, 2)))
bayenv_outlier$chr_pos <- factor(bayenv_outlier$chr_pos, levels = unique(bayenv_outlier$chr_pos))
bayenv_outlier <- bayenv_outlier %>% group_by(chr_pos, region) %>% mutate(group_mean = mean(pop_aaf), adj_aaf = scale(pop_aaf, scale = F))
bayenv_outlier$region_abbr <- case_when(
  bayenv_outlier$region == "West Central Africa" ~ "WC Africa",
  bayenv_outlier$region == "East Central Africa" ~ "EC Africa",
  bayenv_outlier$region == "East Africa" ~ "E Africa",
  bayenv_outlier$region == "Southern Africa" ~ "S Africa",
  bayenv_outlier$region == "South Asia" ~ "S Asia",
  bayenv_outlier$region == "Southeast Asia" ~ "SE Asia")
bayenv_outlier$region_abbr <- factor(bayenv_outlier$region_abbr, levels = c("WC Africa", "EC Africa", "E Africa", "S Africa", "S Asia", "SE Asia"))
bayenv_outlier$region <- factor(bayenv_outlier$region, levels = c("West Central Africa", "East Central Africa", "East Africa", "Southern Africa", "South Asia", "Southeast Asia"))

# Select the 4 SNPs we want to plot
top4 <- bayenv_outlier %>% filter(pos %in% c("59920004", "31210215", "33061192", "40939460"))

top4$chr_pos <- rep(
  c("rs1385812 3:59920004 intronic FHIT",
    "rs1401331217 6:31210215 intergenic HLA-C",
    "rs9461833 6:33061192 intergenic HLA-DPB1",
    "rs9916754 17:40939460 exonic WNK4"), each = 27
)

snp_af_plot <- top4 %>%
  filter(region != "Oceania" & region != "East Asia"
         & group_mean != 0 & group_mean != 1
         ) %>%
  ggplot(aes(x = region, y = adj_aaf, fill = subsistence, color = subsistence)) +
    geom_hline(yintercept = 0, color = "black", linetype = 2) +
    stat_summary(fun.y = mean, geom = "crossbar", position = position_dodge(width = 0.5), fatten = 1.5, width = 0.5, show.legend = F) +
    geom_point(position = position_dodge(width = 0.5), size = 0.7, alpha = 1) +
    scale_discrete_manual(aesthetics = c("color", "fill"), values = c(
      "Hunter-gatherer" = "#C21F84",
      "Agriculturalist" = "#FDAE61",
      "Pastoralist" = "#40E0D0")) +
    labs(x = "", y = "Regional AAF deviation", fill = "", color = "") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    facet_wrap(~chr_pos, scales = "free_y", ncol = 1) +
    theme_classic(base_size = 7) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "top",
          strip.placement = "bottom",
          strip.text.x = element_text(angle = 0),
          strip.text.y = element_text(angle = 0, hjust = 0),
          strip.background = element_blank(),
          #panel.spacing.x = unit(0, "cm"),
          axis.line.x = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

saveRDS(snp_af_plot, paste0(wd, "data/fig5b_snp_af.rds"))
```

**Step 10) Plot the gene body GO results classifying gene sets into broader biological functions.**

`bayenv_analysis.Rmd {go_enrich}` → `fig5c_go_enrich.rds`

```r
library(data.table)
library(ggplot2)
library(scales)
library(tidyverse)
library(khroma)

wd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/results_bayenv/"
muted <- color("muted")

go_child2 <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/gs_data/GO_child2.csv", header = F)
colnames(go_child2) <- c("term_id", "term_name")
muted <- color("muted")

# Load in enrichment results and filter for gene sets of interest 
load(paste0(wd, "data/outlier_GO.RData"))
go_outlier_sig <- rbind(sig_sets %>% filter(ID %in% go_child2$term_id | !grepl("GO", ID)) %>%
      filter(grepl("development|growth|morpho|pattern|hedgehog|hippo|neur|nerv|synap|axo|dendrit|perikaryon|action potential", Description, ignore.case = T)),
                        sig_sets %>% filter(
      grepl("hormone|cortisol|androgen|estrogen|oxytocin|aldosterone|GnRH|diabetes|insulin|endocrine", sig_sets$Description, ignore.case = T) |
      grepl("muscle|actin |Z disc|band|sarcomere|costamere|myo", sig_sets$Description, ignore.case = T) |
      grepl("immun|leuk|vir|IL2|Th1|MHC|inflamm|T cell|recognition|phagosome|efferocytosis", sig_sets$Description, ignore.case = T) |
      grepl("metabol|carbohydrate|lipid|amino acid|catabolic|leucine|lysine|digest|salivary", sig_sets$Description, ignore.case = T) |
      !grepl("GO:", ID))
)

# Separate enrichments into groups of immune, muscle, nervous, development, and metabolic functions
go_outlier_sig$enrich_group <- case_when(
  grepl("immun|leuk|vir|IL2|Th1|MHC|inflamm|T cell|recognition|diabetes|phagosome|efferocytosis", go_outlier_sig$Description, ignore.case = T) ~ "Immune",
  grepl("muscle|actin |Z disc|band|sarcomere|costamere|myo", go_outlier_sig$Description, ignore.case = T) ~ "Muscle",
  grepl("neur|nerv|axo|synap|dendrit|behavior|action potential", go_outlier_sig$Description, ignore.case = T) ~ "Nervous",
  grepl("development|growth|morpho|pattern|hedgehog|hippo", go_outlier_sig$Description, ignore.case = T) ~ "Development",
  grepl("metabol|carbohydrate|lipid|amino acid|catabolic|leucine|lysine|digest|salivary", go_outlier_sig$Description, ignore.case = T) ~ "Metabolic",
  T ~ "Other"
)

# Set up plotting df
go_outlier_sig$logp <- -log10(go_outlier_sig$p.adjust)
go_outlier_sig$gene_ratio <- sapply(go_outlier_sig$GeneRatio, function(x) eval(parse(text = x)))
go_outlier_sig$desc <- str_to_sentence(go_outlier_sig$Description)
go_outlier_sig$desc[c(40, 27, 42, 28, 29)] <- c(
  "MHC class II protein complex",
  "MHC protein complex",
  "Type I diabetes mellitus",
  "Positive regulation of T cell activation",
  "Regulation of T cell activation")
go_outlier_sig$enrich_group[42] <- "Immune"

go_outlier_sig$enrich_group <- factor(go_outlier_sig$enrich_group, levels = rev(c("Development", "Metabolic", "Nervous", "Muscle", "Immune", "Other")))
go_outlier_sig <- go_outlier_sig %>% group_by(enrich_group) %>% arrange(gene_ratio, .by_group = T)
go_outlier_sig$enrich_group <- factor(go_outlier_sig$enrich_group, levels = c("Development", "Metabolic", "Nervous", "Muscle", "Immune", "Other"))
go_outlier_sig$desc <- factor(go_outlier_sig$desc, levels = unique(go_outlier_sig$desc))

# Plot as lollipop graph
go_enrich_plot <- go_outlier_sig %>% filter(enrich_group != "Other") %>%
  ggplot(aes(y = desc, x = gene_ratio, color = enrich_group)) +
    geom_segment(aes(x = 0, y = desc, xend = gene_ratio, yend = desc), show.legend = F) +
    geom_point(aes(color = enrich_group, size = logp)) +
    scale_color_manual(values = muted(9)[c(5:9)]) +
    scale_size(breaks = c(2.5, 5, 7.5),
               labels = as.character(c(2.5, 5,7.5)),
               range = c(0.5, 2)) +
    labs(x = "Gene ratio", y = NULL, size = expression(paste("-log"[10], " (p)"))) +
    guides(color = guide_legend(title = NULL, order = 2), size = guide_legend(order = 1)) +
    theme_classic(base_size = 7) +
    theme(plot.margin = unit(c(2.5, 5.5, 2.5, 5.5), "pt"),
          legend.position = c(0.7,0.5),
          legend.box = "horizontal",
          legend.spacing = unit(0, "mm"),
          legend.key.spacing = unit(0, "mm"),
          axis.text.y = element_text(size = 5))

saveRDS(go_enrich_plot, paste0(wd, "data/fig5c_go_enrich.rds"))
```

**Step 11) Arrange the above plots to create Figure 5.**

`bayenv_analysis.Rmd {arrange plots}`

```r
library(data.table)
library(ggplot2)
library(scales)
library(tidyverse)
library(khroma)
library(ggfastman)
library(ggrepel)
library(ggnewscale)
library(cowplot)

wd <- "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/bayenv/results_bayenv/"
muted <- color("muted")

manhattan_plot <- readRDS(paste0(wd, "data/fig5a_manhattan.rds"))
snp_af_plot <- readRDS(paste0(wd, "data/fig5b_snp_af.rds"))
go_enrich_plot <- readRDS(paste0(wd, "data/fig5c_go_enrich.rds"))

bayenv_figures <- plot_grid(
    plot_grid(manhattan_plot, nrow = 1, labels = c("a"), label_size = 8),
    plot_grid(snp_af_plot, go_enrich_plot, labels = c("b", "c"), rel_widths = c(0.7, 1), nrow = 1, label_size = 8),
  rel_heights = c(0.7,1), nrow = 2)

save_plot(paste0(wd, "figures/bayenv_figures_v3.pdf"), bayenv_figures, base_width = 7, base_height = 6)
```
