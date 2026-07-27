# 1) Whole-genome mapping and calling procedures

This page details the procedures taken to generate the VCFs used for further downstream analysis. 

## Preprocessing and calling

Preprocessing was performed using the [Sarek germline variant calling pipeline](https://nf-co.re/sarek/3.3.2) with default parameters, which follows the GATK joint germline short variant discovery workflow. Raw FastQs were aligned to the GRCh37 reference genome (`--genome GATK.GRCh37` in sarek), duplicates were marked, and single sample calling was performed to produce a GVCF for each sample. After all samples reached this stage, joint calling was applied across the entire cohort to generate a single raw callset, with VQSR applied at the 99.9 tranche for SNPs and 99.0 tranche for indels. Please refer to the nf-core documentation for how to run this pipeline. 

Raw VCF file: `germline_PASS_indel_SNP_recalibrated_99.9.selectvariants.vcf.gz`

## Generating coverage statistics

**Step 1) Use mosdepth QC outputs from the Sarek pipeline to get sample coverage genome-wide.**

`mosdepth_average_coverage` → `freeze2_coverage.txt`

```bash
cd /project/lbarreiro/USERS/evanwu/
out="/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/coverage/freeze2_coverage.txt"
readarray -d '' regions < <(find ./ -name *.mosdepth.region.dist.txt -print0 | sort)
readarray -d '' summaries < <(find ./ -name *.mosdepth.summary.txt -print0 | sort)
n=${#regions[@]}

# Get data: sample, average depth genome-wide, proportion of bases covered at 0X, 1X, 10X, 20X, 30X
printf "sample\tavg_depth\t0X\t1X\t10X\t20X\t30X\n" > $out

for (( i = 0; i < $n-1; i++ )); do
        sample=$(basename ${regions[i]})
        sample=${sample%%\.*}

        avg=$(tail -n 1 ${summaries[i]} | awk '{print $4}')

        X30=$(grep 'total' ${regions[i]} | grep '\b30\s' | awk '{print$3}')
        X20=$(grep 'total' ${regions[i]} | grep '\b20\s' | awk '{print$3}')
        X10=$(grep 'total' ${regions[i]} | grep '\b10\s' | awk '{print$3}')
        X1=$(grep 'total' ${regions[i]} | grep '\b1\s' | awk '{print$3}')
        X0=$(grep 'total' ${regions[i]} | grep '\b0\s' | awk '{print$3}')

        printf "$sample\t$avg\t$X0\t$X1\t$X10\t$X20\t$X30\n" >> $out
done
```

## VCF filtering and QC

**Step 1) Normalize VCF, select biallelic SNPs, exclude samples with <10X average genome-wide coverage.**

`bcftools_normalize.sh` → `freeze2_normalized{_biallelic,_indel}.vcf.gz`

```bash
sample_rename="freeze2_vcfID_to_studyID_key.txt"
fasta="/project/lbarreiro/SHARED/REFERENCES/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta"

echo "Normalizing VCF"
bcftools reheader -s $sample_rename --threads 10 germline_PASS_indel_SNP_recalibrated_99.9.selectvariants.vcf.gz | \
  bcftools norm -m - -c x -f $fasta --threads 10 -O z -o freeze2_normalized.vcf.gz # this splits multiallelic sites into multiple biallelic entries but we actually don't want that
bcftools index -t --threads 10 freeze2_normalized.vcf.gz

echo "Selecting Biallelic SNPs"
bcftools view -m2 -M2 -v snps -S samples.txt --threads 10 freeze2_normalized.vcf.gz -O z -o freeze2_normalized_biallelic.vcf.gz
bcftools index -t --threads 10 freeze2_normalized_biallelic.vcf.gz

echo "Selecting Indels"
bcftools view -v indels -S samples.txt --threads 10 freeze2_normalized.vcf.gz -O z -o freeze2_normalized_indel.vcf.gz
bcftools index -t --threads 10 freeze2_normalized_indel.vcf.gz
```

**Step 2) Baseline hard filters / recoding:**
- Within Info field, filter out RMSMappingQuality (MQ) < 40.
- Within Info field, filter out Quality by Depth (QD) < 2.0.
- Within Info field, filter out Fisher’s exact test for Strand bias (FS) > 60.
- Within Format field, GQ < 30 is set as missing.
- Within Format field, a depth (DP) < 10 or (ii) DP > 3X the average sample coverage is set as missing.

`cyvcf_hardfilter.py ${chr}` → `freeze2_biallelic_filtered.${chr}.vcf.gz`

```python
import sys
from cyvcf2 import VCF, Writer
import numpy as np
import pandas as pd

# Initialize vcf parsers
wd = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/raw_vcfs/'
vcf_file = wd + 'freeze2_normalized_biallelic.vcf.gz'
vcf = VCF(vcf_file, strict_gt=True, threads=10)
writer = Writer('baseline_filters/freeze2_filtered.' + sys.argv[1] + '.vcf.gz', vcf)

# Read in per-sample mean depths
avg_dp = pd.read_csv(wd + 'depth_summary.txt', sep='\t')
avg_dp['DP_MAX'] = avg_dp['MEAN_DEPTH'] * 3 # Filter at 3X average sample DP
avg_dp = np.array(avg_dp['DP_MAX']).flatten()

# Read in each VCF site per chromosome
for v in vcf(sys.argv[1]):
    output_record = True

    # Filter based on INFO
    if (v.INFO.get('MQ') is None or v.INFO.get('QD') is None or v.INFO.get('FS') is None):
            output_record = False
    elif (v.INFO.get('MQ') < 40.0 or v.INFO.get('QD') < 2.0 or v.INFO.get('FS') > 60.0):
            output_record = False

    # Set GT to missing based on FORMAT
    else:
        set_miss = np.nonzero(np.any((v.gt_quals < 30, v.gt_depths < 10, np.greater(v.gt_depths, avg_dp)), axis=0))[0]
        if len(set_miss) != 0:
            for i in set_miss:
                v.genotypes[i] = [-1]*v.ploidy + [False]
                v.genotypes = v.genotypes

    # Write passed sites
    if output_record: writer.write_record(v)

vcf.close(); writer.close()
```

**Step 3) Calculate the average missingness within each of the 8 calling cohorts and remove sites where the median cohort missingness is > 10%.**

`cyvcf_cohort_missingness.py ${chr}` → `freeze2_biallelic_filtered_site_0.1_cohort_missing.${chr}.vcf.gz`  and `per_cohort_missingness.${chr}.txt.gz`

```python
import sys
import csv
from cyvcf2 import VCF, Writer
import numpy as np
import pandas as pd

# Initialize vcf parsers
vcf_file = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/raw_vcfs/baseline_filters/freeze2_biallelic_filtered.' + sys.argv[1] + '.vcf.gz'
vcf = VCF(vcf_file, strict_gt=True, threads=10)
writer = Writer('baseline_filters_site_0.1_cohort_missing/freeze2_biallelic_filtered_site_0.1_cohort_missing.' + sys.argv[1] + '.vcf.gz', vcf)

# Read in population metadata
metadata = pd.read_csv('/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v3.txt', sep='\t')
metadata = metadata[['study_ID','cohort']]
cohorts_meta = metadata.groupby('cohort')['study_ID'].apply(list)

for i in cohorts_meta.index:
    cohorts_meta[i] = [s for s, sample in enumerate(vcf.samples) if sample in set(cohorts_meta[i])]

# Open file to save missingness info
pcm = open('per_cohort_missingness/per_cohort_missingness.' + sys.argv[1] + '.txt', "a")
pcm_writer = csv.writer(pcm, delimiter='\t')
pcm_writer.writerow(['CHROM','POS','REF','ALT', 'Median'] + cohorts_meta.index.tolist())

filter_level = 0.1
num_cohorts = len(cohorts_meta)

# Read in each VCF site
for record in vcf:
    output_record = True
    per_cohort_miss = np.zeros(num_cohorts)

    # For each cohort calculate the mean missingness across samples
    for c, cohort in enumerate(cohorts_meta.index):
        per_cohort_miss[c] = np.mean(record.gt_types[np.array(cohorts_meta[cohort]).astype(int)] == 2)

    # Filter site if median/mean > 0.1
    med = np.median(per_cohort_miss)
    if med > filter_level: output_record = False
    # if np.mean(per_cohort_miss) > filter_level: output_record = False
    pcm_writer.writerow([record.CHROM, record.POS, record.REF, record.ALT[0], med] + list(per_cohort_miss))

    if output_record: writer.write_record(record)

vcf.close(); writer.close(); pcm.close()
```

**Step 4) Concatenate the per-chromosome filtered VCFs**

`bcftools_concat.sh` → `freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing.vcf.gz`

```bash
#export TMPDIR=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/raw_vcfs/tmp

bcftools concat freeze2_biallelic_filtered_site_0.1_cohort_missing.{1..11}.vcf.gz -O z --threads 10 -o freeze2_tmp1.vcf.gz
bcftools concat freeze2_biallelic_filtered_site_0.1_cohort_missing.{12..22}.vcf.gz -O z --threads 10 -o freeze2_tmp2.vcf.gz
bcftools concat freeze2_tmp1.vcf.gz freeze2_tmp2.vcf.gz -O u --threads 10 | \
  bcftools sort -o freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing.vcf.gz
bcftools index -t freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing.vcf.gz --threads 10
```

**Step 5) Remove multiallelic sites (> 1 variant observed at a site) and monomorphic sites (introduced when we recoded genotypes as missing based on filters, or later steps when removing samples). This constitutes the Baseline VCF.**  

`bcftools_multiallelic.sh` → `freeze2_biallelic_polymorphic.vcf.gz`

```bash
VCF=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing/freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing.vcf.gz

# Extract and remove multiallelic sites
bcftools norm $VCF -m + --threads 10 -O z -o freeze2_multiallelic.vcf.gz
bcftools index -t --threads 10 freeze2_multiallelic.vcf.gz
bcftools view freeze2_multiallelic.vcf.gz -m3 --threads 10 -O z -o multiallelic_sites.vcf.gz
bcftools index -t --threads 10 multiallelic_sites.vcf.gz
bcftools view -m2 -M2 freeze2_multiallelic.vcf.gz --threads 10 -O z -o freeze2_biallelic.vcf.gz
bcftools index -t --threads 10 freeze2_biallelic.vcf.gz

# Extract and remove monomorphic sites from the biallelic VCF, get list of all sites to exclude by combining the multiallelic and monomorphic VCFs
bcftools view freeze2_biallelic.vcf.gz -c 1:minor --threads 10 -O z -o freeze2_biallelic_polymorphic.vcf.gz
bcftools index -t --threads 10 freeze2_biallelic_polymorphic.vcf.gz

bcftools view freeze2_biallelic.vcf.gz -C 0:minor --threads 10 -O z -o monomorphic_sites.vcf.gz
bcftools index -t --threads 10 monomorphic_sites.vcf.gz

bcftools concat multiallelic_sites.vcf.gz monomorphic_sites.vcf.gz -n --threads 10 | bcftools sort -T tmp -O z -o freeze2_multiallelic_monomorphic_sites.vcf.gz
bcftools index -t --threads 10 freeze2_multiallelic_monomorphic_sites.vcf.gz

bcftools query freeze2_multiallelic_monomorphic_sites.vcf.gz -f '%CHROM\t%POS\n' -o multiallelic_monomorphic_sites.txt
```

**Step 6) Calculate the regional MAF statistics for each site.**

`cyvcf_regional_maf.py ${chr}`  → `per_region_maf.${chr}.txt.gz` 

```python
import sys
import csv
from cyvcf2 import VCF, Writer
import numpy as np
import pandas as pd

# Initialize vcf parsers
vcf_file = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/raw_vcfs/baseline_filters_site_0.1_cohort_missing/freeze2_biallelic_filtered_site_0.1_cohort_missing.' + sys.argv[1] + '.vcf.gz'
vcf = VCF(vcf_file, strict_gt=True, gts012=True, threads=10)

# Read in population metadata
metadata = pd.read_csv('/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v3.txt', sep='\t')
metadata = metadata[['study_ID','region']]
metadata['region'] = [i[-1] for i in metadata['region'].str.split()]
regions_meta = metadata.groupby('region')['study_ID'].apply(list)

for i in regions_meta.index:
    regions_meta[i] = [s for s, sample in enumerate(vcf.samples) if sample in set(regions_meta[i])]

# Open file to save MAF info
prm = open('per_region_maf/per_region_maf.' + sys.argv[1] + '.txt', "a")
prm_writer = csv.writer(prm, delimiter='\t')
prm_writer.writerow(['CHROM','POS','REF','ALT', 'MAF'] + regions_meta.index.tolist())

num_regions = len(regions_meta)

# Read in each VCF site
for record in vcf:
    output_record = True
    per_region_maf = np.zeros(num_regions)

    # For each region calculate the MAF across samples and write out the stats
    for r, region in enumerate(regions_meta.index):
        gts = record.gt_types[np.array(regions_meta[region]).astype(int)]
        per_region_maf[r] = np.mean(gts[gts != 3]) / 2

    prm_writer.writerow([record.CHROM, record.POS, record.REF, record.ALT[0], record.aaf] + list(per_region_maf))

vcf.close(); prm.close()

```

**Step 7) Calculate the population MAF statistics for each site.**

`cyvcf_population_maf.py ${chr}`  → `per_population_maf.${chr}.txt.gz`

```python
import sys
import csv
from cyvcf2 import VCF, Writer
import numpy as np
import pandas as pd

# Initialize vcf parsers
vcf_file = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing/per_chr/freeze2_biallelic_filtered_site_0.1_cohort_missing.' + sys.argv[1] + '.vcf.gz'
vcf = VCF(vcf_file, strict_gt=True, gts012=True, threads=10)
#writer = Writer('chroms/freeze1_filtered_0.1_cohort_missing.' + sys.argv[1] + '.vcf.gz', vcf)

# Read in population metadata
metadata = pd.read_csv('/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v2.txt', sep='\t')
metadata = metadata[['study_ID','pop_code']]
metadata['pop_code'] = [i[-1] for i in metadata['pop_code'].str.split()]
pops_meta = metadata.groupby('pop_code')['study_ID'].apply(list)

for i in pops_meta.index:
    pops_meta[i] = [s for s, sample in enumerate(vcf.samples) if sample in set(pops_meta[i])]
print(pops_meta)

# Open file to save MAF info
ppm = open('per_population_maf/per_population_maf.' + sys.argv[1] + '.txt', "a")
ppm_writer = csv.writer(ppm, delimiter='\t')
ppm_writer.writerow(['CHROM','POS','REF','ALT', 'MAF'] + pops_meta.index.tolist())

#filter_level = 0.1
num_pops = len(pops_meta)

# Read in each VCF site
for record in vcf:
    output_record = True
    per_pop_maf = np.zeros(num_pops)

    # For each pop calculate the MAF across samples
    for r, pop in enumerate(pops_meta.index):
        gts = record.gt_types[np.array(pops_meta[pop]).astype(int)]
        per_pop_maf[r] = np.mean(gts[gts != 3]) / 2

    # Filter site if median/mean > 0.1
    #med = np.median(per_cohort_miss)
    #if med > filter_level: output_record = False
    # if np.mean(per_cohort_miss) > filter_level: output_record = False
    ppm_writer.writerow([record.CHROM, record.POS, record.REF, record.ALT[0], record.aaf] + list(per_pop_maf))

    #if output_record: writer.write_record(record)

vcf.close(); ppm.close()#; writer.close()
```

**Step 8) Generate overall stats on a VCF and compare with dbSNP build 156.** 

`gatk_VCF_metrics.sh` → `dbsnp156_filtered_site_0.1_cohort_missing_VCF_metrics.variant_calling_{detail,summary}_metrics`

```bash
echo "Collect VCF Metrics"
gatk CollectVariantCallingMetrics \
  --DBSNP dbsnp_156_autosomes.b37.vcf.gz \
  -I /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/raw_vcfs/baseline_filters_site_0.1_cohort_missing/freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing.vcf.gz \
  -O dbsnp156_filtered_site_0.1_cohort_missing_VCF_metrics \
  --THREAD_COUNT 10
```

## Variant annotation

**Step 1) Generate annotations using Annovar per-chromosome to get gene, gene region, SNP consequence, clinical significance, etc. annotations**

 `run_table_annovar_per_chr.sh` → `chr${chr}_anno.hg19_multianno.txt`

```bash
#!/bin/bash

table_annovar.pl /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing/per_chr/freeze2_biallelic_filtered_site_0.1_cohort_missing.${chr}.vcf.gz humandb/ -buildver hg19 -out chr${chr}_anno -remove -protocol refGene,cytoBand,exac03,avsnp150,clinvar_20221231 -operation g,r,f,f,f -vcfinput -polish --thread 6
```

## Removing heterozygosity and relatedness outliers

**Step 1) Calculate heterozygosity using PLINK.**

`plink_heterozygosity.sh` → `freeze2_baseline_filters_site_0.1_cohort_missing.het`

```bash
#!/bin/bash

plink --vcf freeze2_biallelic_polymorphic.vcf.gz --make-bed --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing

plink --bfile freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing --het --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing
```


**Step 2) Assess heterozygosity within each sample and designate outliers as samples with heterozygosity rates above or below 3SD of the population mean.**

`heterozygosity.R` → `heterozygosity_outliers_to_exclude.tsv`

```R
######################################
####  load packages and datasets  ####
######################################
library(ggplot2)
library(ggrepel)
library(dplyr)

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/heterozygosity/baseline_filters_site_0.1_cohort_missing")
raw <- read.csv("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/updated/baseline_filters_site_0.1_cohort_missing/freeze2_baseline_filters_site_0.1_cohort_missing.het", sep="")
metadata <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v3.txt")

#####################################################
####  calculate rate of heterozygosity and plot  ####
#####################################################
# (N(NM)-O(HOM))/N(NM)
raw$het_rate <- ((raw$N.NM.)-(raw$O.HOM.))/raw$N.NM.

# merge with metadata
df <- merge(raw, metadata, by.x = "IID", by.y = "study_ID", all.x = TRUE)

# function for number of observations (FID counts per group)
give.n <- function(x){
  return(c(y = median(x)*1.16, label = length(x)))
  # experiment with the multiplier to find the perfect position
}

##################################################
#### create per-population summary statistics ####
##################################################

# create table of summary statistics
population_ave <- df %>% group_by(FID) %>% summarize(mean_het_rate = mean(het_rate, na.rm=TRUE),
                                                     median_het_rate = median(het_rate, na.rm=TRUE),
                                                     sd_het_rate = sd(het_rate, na.rm=TRUE),
                                                     above_3sd = mean_het_rate + 3*sd_het_rate,
                                                     below_3sd = mean_het_rate - 3*sd_het_rate)

population_counts <- df %>% group_by(FID) %>% count()

population_ave <- merge(population_ave, population_counts, by = "FID", all.x = TRUE)
population_ave <- population_ave %>% select(FID, n, mean_het_rate, median_het_rate, sd_het_rate, above_3sd, below_3sd)


##############################################################################
#### identify individuals that are 3SD above or below the population mean ####
##############################################################################
df_with_pop_stats <- df %>% left_join(population_ave, by = "FID")

# samples to potentially exclude
excluded_samples <- df_with_pop_stats %>% mutate(exclude = (het_rate >= above_3sd | het_rate <= below_3sd)) %>% filter(exclude)

# print out their info quickly
excluded_samples %>% select(IID, FID, het_rate, mean_het_rate, above_3sd, below_3sd)

# save list of samples to exclude
write.table(excluded_samples$IID, file = "heterozygosity_outliers_to_exclude.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)

# samples to retain
df_exclude_outliers <- df_with_pop_stats %>%
  mutate(exclude = (het_rate >= above_3sd | het_rate <= below_3sd)) %>% filter(!exclude)

save.image(file="heterozygosity.RData")
```

**Step 3) Remove heterozygosity outlier individuals and check to remove monomorphic variants that may have arisen due to individuals being removed. This constitutes the Level 1 VCF.**

`bcftools_nohet.sh` → `freeze2_biallelic_polymorphic_nohet.vcf.gz`

```bash
HET=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet/heterozygosity_outliers_to_exclude.tsv

# Remove het outliers and then monomorphic
bcftools view freeze2_biallelic_polymorphic.vcf.gz -S "^${HET}" --threads 10 -O u | \
  bcftools view -c 1:minor --threads 10 -O z -o freeze2_biallelic_polymorphic_nohet.vcf.gz
bcftools index -t --threads 10 freeze2_biallelic_polymorphic_nohet.vcf.gz
```

**Step 4) Calculate relatedness between all individuals using KING.**

`KING.sh` → `freeze2_baseline_filters_site_0.1_cohort_missing_nohet.kin`

```bash
#!/bin/bash

input_dir=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/updated/baseline_filters_site_0.1_cohort_missing_nohet

prefix=freeze2_baseline_filters_site_0.1_cohort_missing_nohet

# Load in VCF with plink and use king to calculate relatedness
plink --bfile --vcf freeze2_biallelic_polymorphic_nohet.vcf.gz --make-bed --out ${input_dir}/${prefix}
king -b ${input_dir}/${prefix}.bed --related --degree 2 --prefix ${prefix}

# Filter KING outputs based on Inf type for 2nd degree relatedness or higher
head -n 1 ${prefix.kin} > within_family.txt
awk '($15 == "2nd" || $15 == "FS" || $15 == "PO") {print $0}' ${prefix}.kin >> within_family.txt
```

**Step 5) Prepare KING outputs into NAToRA format.**

`make_natora_inputs.R` → `natora_using_Inf_type.tsv`

```R
setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/relatedness/baseline_filters_site_0.1_cohort_missing_nohet")

# import the within KING outputs (no between relatives detected) that are already filtered for 2nd degree relatives based on Inf type
within <- read.delim("within_family.txt")

# make natora file format
within_natora <- cbind(within$ID1, within$ID2, within$Kinship)
natora_ready <- rbind.data.frame(within_natora)

write.table(natora_ready, file = "natora_using_Inf_type.tsv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep="\t")
min(natora_ready$V3)
```

**Step 6) Run NAToRA to get list of individuals to remove.**

`natora.sh` → `freeze2_baseline_filters_site_0.1_cohort_missing_nohet_natora_using_Inf_type_toRemove.txt`

```bash
!#/bin/bash

# Run NAToRA
python NAToRA_Public.py --input natora_using_Inf_type.tsv -d 2 -o freeze2_baseline_filters_site_0.1_cohort_missing_nohet_natora_using_Inf_type
```

**Step 7) Remove related individuals up to 2nd degree and check to remove monomorphic variants that may have arisen due to individuals being removed. This constitutes the Level 2 VCF.**

`bcftools_nohet_nokin.sh` → `freeze2_biallelic_polymorphic_nohet_nokin.vcf.gz`

```bash
KIN=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet_nokin/freeze2_baseline_filters_site_0.1_cohort_missing_nohet_natora_using_Inf_type_toRemove.txt

# Remove kin outliers and then monomorphic
bcftools view freeze2_biallelic_polymorphic_nohet.vcf.gz -S "^${KIN}" --threads 10 -O u | \
  bcftools view -c 1:minor --threads 10 -O z -o freeze2_biallelic_polymorphic_nohet_nokin.vcf.gz
bcftools index -t --threads 10 freeze2_biallelic_polymorphic_nohet_nokin.vcf.gz
```

## Phasing 

**Step 1) Phase using SHAPEIT5 in reference-free mode. Genetic maps are downloaded from [1kG](https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html).**

`shapeit5_reffree.sbatch ${chr}` → `chr${chr}_shapeit5_noref.bcf`

```bash
#!/bin/bash

input_vcf=freeze2_biallelic_polymorphic_nohet_nokin.vcf.gz

gmap=/project/lbarreiro/USERS/bridget/references/maps/impute2-1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt

source activate shapeit5

# Use SHAPEIT5 to phase reference-free
SHAPEIT5_phase_common --input ${input_vcf} --region ${chr} --map ${gmap} --output chr${chr}_shapeit5_noref.bcf --thread 8
```

**Step 2) Concatenate all phased chromosomes. This constitutes the Level 3 VCF.**

`bcftools_concat.sh` → `freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes.vcf.gz`

```bash
#!/bin/bash

bcftools concat chr{1..22}_shapeit5_noref.bcf -o freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes.vcf.gz -Oz

bcftools index freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes.vcf.gz
```

**Step 3) Get the list of all sites, reformat, and split into chunks for parallel submission**

`extract_positions.sh` → `chr${chr}_shapeit5_noref.bcf`

```bash
#!/bin/bash

mkdir -p chunks 

# Split by chromosome and format as chr1:17452-17452
for chr in {1..22}; do
    bcftools query -f '%CHROM\t%POS\n' "${VCFDIR}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" > allpositions_${chr}.txt

    awk -v chr="chr${chr}" '$1 == chr {print $1":"$2"-"$2}' allpositions_${chr}.txt \
        > anc_state/chr${chr}.txt
    echo "chr${chr}: $(wc -l < chr${chr}.txt) positions"
done

# Split into chunks of 100K lines
for chr in {1..22}; do
    split -l 100000 --additional-suffix=".chunk" \
        chr${chr}.txt \
        chunks/chr${chr}.txt_
done

echo "Total chunks:"
ls chunks/*.chunk | wc -l
```

**Step 4) Extract the ancestral alleles across all chunk files using the [EPO human ancestral reference Fastas](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2).**

`get_ancestral.sh ${chunk}` → `${chunk}.ancestral`

```bash
#!/bin/bash

ANCDIR=/project/lbarreiro/USERS/bridget/references/human_ancestor_GRCh37_e59

chr=$( echo ${chunk} | cut -d"." -f1 | cut -d"r" -f2 )

# Write ancestral allele to file
/software/samtools-1.22.1-el8-x86_64/bin/samtools faidx ${ANCDIR}/homo_sapiens_ancestor_${chr}.fa -r ${chunk}.rename | \
       paste - - > ${chunk}.ancestral
```

**Step 5) Combine the per-chunk outputs and filter for high and low-confidence calls.**

`filter_ancestral.sh ${chunk}` → `${chunk}.ancestral`

```bash
# Create a directory for the combined files (optional)
mkdir -p combined
cd combined

# Loop through each chromosome
for chr_num in {1..22}; do
    # Find all chunked ancestral files for the chromosome, sorted alphabetically
    chunk_files=$(ls chunks/chr${chr_num}.txt_*.chunk.ancestral | sort)

    # Combine them into a single file, preserving order
    cat $chunk_files > chr${chr_num}.txt.ancestral
done

# Retain low- and high-confidence ancestral calls
for chr in {1..22}; do
    awk '$2 ~ /^[ACTGactg]$/' chr${chr}.txt.ancestral > chr${chr}.ancestral.filtered
done

# Add rsID as chrXposXXXXX format
for chr in {1..22}; do
    awk -F'\t' '$2 ~ /^[ACTGactg]$/ {
        gsub(">", "");
        split($1, arr, ":|-");
        printf("%spos%s\t%s\n", arr[1], arr[2], $2);
    }' chr${chr}.ancestral.filtered > tmpfile && mv tmpfile chr${chr}.ancestral.filtered;
done

# Remove duplicate positions
for chr in {1..22}; do 
    awk 'NR==FNR{count[$1]++; next} count[$1]==1' chr${chr}.ancestral.filtered chr${chr}.ancestral.filtered > chr${chr}.ancestral.filtered.norepeats
done

# Convert the lower case letters to all upper cases
for chr in {1..22}; do
    awk '{print $1, toupper($2)}' chr${chr}.ancestral.filtered.norepeats > chr${chr}.ancestral.filtered.norepeats.upper
done
```

**Step 6) Using the list of ancestral sites, recode the VCF with ancestral allele as REF and derived allele as ALT and remove sites without ancestral information. This constitutes the Level 4 VCF.**

`ancestral_recode.sh` → `freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes_anc_inclc_clean.vcf.gz`

```bash
#!/bin/bash

# Annotate VCF with fake rsIDs that matches the outputs above
bcftools annotate --set-id '\chr%CHROM\pos%POS' freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes.vcf.gz -Oz -o freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes_fakeIDs.vcf.gz 

# STEP 1 FILTERING: ONLY RETAIN VARIANTS IN THE VCF THAT ARE PRESENT IN chr${chr}.ancestral.filtered.norepeats.upper

zcat freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes_fakeIDs.vcf.gz | \
awk 'NR==FNR { keys[$1]; next }
     /^#/ { print; next }
     ($3 in keys) { print }' <(for chr in {1..22}; do
          awk '{print $1}' combined/chr${chr}.ancestral.filtered.norepeats.upper
      done) - | \
bgzip -c > freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes_ancestralonly_inclc.vcf.gz

# STEP 2 FILTERING: ONLY RETAIN VARIANTS WHERE THE SECOND COLUMN OF chr${chr}.ancestral.filtered.norepeats MATCHES EITHER COLUMNS 4 OR 5 OF THE VCF

zcat freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes_ancestralonly_inclc.vcf.gz | \
awk 'NR==FNR { 
         ancestral[$1] = $2; 
         next 
     }
     /^#/ { 
         print; 
         next 
     }
     {
         key = $3;
         if ( (key in ancestral) && (ancestral[key] == $4 || ancestral[key] == $5) )
             print
     }' <(for chr in {1..22}; do
              awk '{print $1, $2}' combined/chr${chr}.ancestral.filtered.norepeats.upper
          done) - | \
bgzip -c > freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes_anc_inclc.vcf.gz

# finally use plink to recode ancestral/derived information
cat chr{1..22}.ancestral.filtered.norepeats.upper > ancestral_all.norepeats.upper

plink2 \
  --vcf freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes_anc_inclc.vcf.gz \
  --ref-allele force /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet_nokin_SHAPEIT5noref/anc_state/combined/inc_low_confidence/ancestral_all.norepeats.upper 2 1 \
  --recode vcf \
  --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes_anc_inclc_clean

# compress and index
bgzip -c freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes_anc_inclc_clean.vcf > freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes_anc_inclc_clean.vcf.gz

bcftools index freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_shapeit5_noref_autosomes_anc_inclc_clean.vcf.gz
```
