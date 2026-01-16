# Mapping and genotype calling procedures

# **Preprocessing (FastQ → GVCF)**

Used Sarek germline variant calling pipeline: [https://nf-co.re/sarek](https://nf-co.re/sarek), [https://github.com/nf-core/sarek/tree/3.3.2](https://github.com/nf-core/sarek/tree/3.3.2)

Aligned to GRCh37 reference and subsetted for WGS callable regions from GATK: `/project/lbarreiro/SHARED/REFERENCES/Homo_sapiens/GATK/GRCh37/Annotation/intervals/wgs_calling_regions_Sarek.list`

![Untitled](Untitled.png)

Tool parameters are the defaults provided by Sarek

`BCFTOOLS_STATS:
bcftools: 1.16-90-g2191405
BWAMEM1_MEM:
bwa: 0.7.17-r1188
samtools: '1.16'
CRAM_TO_BAM:
samtools: '1.16'
CRAM_TO_BAM_RECAL:
samtools: '1.16'
CREATE_INTERVALS_BED:
gawk: 4.0.2
CUSTOM_DUMPSOFTWAREVERSIONS:
python: 3.9.7
yaml: '6.0'
FASTP:
fastp: 0.23.2
FASTQC:
fastqc: 0.11.9
GATK4_APPLYBQSR:
gatk4: 4.2.6.1
GATK4_BASERECALIBRATOR:
gatk4: 4.2.6.1
GATK4_GATHERBQSRREPORTS:
gatk4: 4.2.6.1
GATK4_HAPLOTYPECALLER:
gatk4: 4.2.6.1
GATK4_MARKDUPLICATES:
gatk4: 4.2.6.1
samtools: '1.16'
INDEX_CRAM:
samtools: '1.16'
INDEX_MARKDUPLICATES:
samtools: '1.16'
MERGE_CRAM:
samtools: '1.16'
MERGE_HAPLOTYPECALLER:
gatk4: 4.2.6.1
MOSDEPTH:
mosdepth: 0.3.3
SAMTOOLS_STATS:
samtools: '1.16'
TABIX_BGZIPTABIX_INTERVAL_SPLIT:
tabix: '1.14'
TABIX_DBSNP:
tabix: '1.14'
Workflow:
Nextflow: 22.10.7
nf-core/sarek: 3.1.2`

# Joint Calling (GVCF → VCF)

GATK 4.2.3.0: SelectVariants → GenomicsDBImport → GenotypeGVCFs → VariantRecalibrator (SNP and Indel) → ApplyVQSR (99.9 tranche) → SelectVariants PASS 

# VCF processing

Working directory for scripts/files: `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/raw_vcfs` 

Final directory: `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf` 

Raw VCF file: `germline_PASS_indel_SNP_recalibrated_99.9.selectvariants.vcf.gz`

**Step 1) Normalize VCF and split into biallelic version vs. indel version, exclude samples with <10X average genome-wide coverage (based on `samples.txt`, refer to `../freeze2_coverage_corrected.txt`).**

`bcftools_normalize.sh`  → `freeze2_normalized{_biallelic,_indel}.vcf.gz`

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

- Within Info field, filter out RMSMappingQuality (MQ) < 40
- Within Info field, filter out Quality by Depth (QD) < 2.0
- Within Info field, filter out Fisher’s exact test for Strand bias (FS) > 60
- Within Format field, GQ < 30 is set as missing
- Within Format field, a depth (DP) < 10 or (ii) DP > 3X the average sample coverage is set as missing (refer to `depth_summary.txt`)

`cyvcf_submit.sh`  and `cyvcf_hardfilter.py`  → `baseline_filters/per_chr/freeze2_biallelic_filtered.{chr}.vcf.gz`

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

**Step 3) Calculate the average missingness within each of the 8 calling cohorts and remove sites where the median cohort missingness is  >10% (more robust to cohort batch effects)**

`cyvcf_submit.sh`  and `cyvcf_cohort_missingness.py`  → `baseline_filters_site_0.1_cohort_missing/per_chr/freeze2_biallelic_filtered_site_0.1_cohort_missing.{chr}.vcf.gz`  and `per_cohort_missingness/per_cohort_missingness.{chr}.txt.gz`

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
metadata = pd.read_csv('/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v1.txt', sep='\t')
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

**Step 5) Remove multiallelic sites (>1 variant observed at a site) and monomorphic sites (introduced when we recoded genotypes as missing based on filters, or later steps when removing samples)**  

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

**Step 6) Remove individuals that are heterozygosity and relatedness outliers and check to remove monomorphic variants that may have arisen due to variant individuals being removed**

`bcftools_nohet_nokin.sh` → `freeze2_biallelic_polymorphic_nohet.vcf.gz` and `freeze2_biallelic_polymorphic_nohet_nokin.vcf.gz`

```bash
HET=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet/heterozygosity_outliers_to_exclude.tsv
KIN=/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing_nohet_nokin/freeze2_baseline_filters_site_0.1_cohort_missing_nohet_natora_using_Inf_type_toRemove.txt

# Remove het outliers and then monomorphic
bcftools view freeze2_biallelic_polymorphic.vcf.gz -S "^${HET}" --threads 10 -O u | \
  bcftools view -c 1:minor --threads 10 -O z -o freeze2_biallelic_polymorphic_nohet.vcf.gz
bcftools index -t --threads 10 freeze2_biallelic_polymorphic_nohet.vcf.gz

# Remove kin outliers and then monomorphic
bcftools view freeze2_biallelic_polymorphic_nohet.vcf.gz -S "^${KIN}" --threads 10 -O u | \
  bcftools view -c 1:minor --threads 10 -O z -o freeze2_biallelic_polymorphic_nohet_nokin.vcf.gz
bcftools index -t --threads 10 freeze2_biallelic_polymorphic_nohet_nokin.vcf.gz
```

**Step 7) Calculate the regional and population alternate AF statistics for each site**

`cyvcf_submit.sh`  and `cyvcf_regional_maf.py`  → `per_region_maf/per_region_maf.{chr}.txt.gz` 

```python
import sys
import csv
from cyvcf2 import VCF, Writer
import numpy as np
import pandas as pd

# Initialize vcf parsers
vcf_file = '/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/raw_vcfs/baseline_filters_site_0.1_cohort_missing/freeze2_biallelic_filtered_site_0.1_cohort_missing.' + sys.argv[1] + '.vcf.gz'
vcf = VCF(vcf_file, strict_gt=True, gts012=True, threads=10)
#writer = Writer('chroms/freeze1_filtered_0.1_cohort_missing.' + sys.argv[1] + '.vcf.gz', vcf)

# Read in population metadata
metadata = pd.read_csv('/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v1.txt', sep='\t')
metadata = metadata[['study_ID','region']]
metadata['region'] = [i[-1] for i in metadata['region'].str.split()]
regions_meta = metadata.groupby('region')['study_ID'].apply(list)

for i in regions_meta.index:
    regions_meta[i] = [s for s, sample in enumerate(vcf.samples) if sample in set(regions_meta[i])]

# Open file to save MAF info
prm = open('per_region_maf/per_region_maf.' + sys.argv[1] + '.txt', "a")
prm_writer = csv.writer(prm, delimiter='\t')
prm_writer.writerow(['CHROM','POS','REF','ALT', 'MAF'] + regions_meta.index.tolist())

#filter_level = 0.1
num_regions = len(regions_meta)

# Read in each VCF site
for record in vcf:
    output_record = True
    per_region_maf = np.zeros(num_regions)

    # For each region calculate the MAF across samples
    for r, region in enumerate(regions_meta.index):
        gts = record.gt_types[np.array(regions_meta[region]).astype(int)]
        per_region_maf[r] = np.mean(gts[gts != 3]) / 2

    # Filter site if median/mean > 0.1
    #med = np.median(per_cohort_miss)
    #if med > filter_level: output_record = False
    # if np.mean(per_cohort_miss) > filter_level: output_record = False
    prm_writer.writerow([record.CHROM, record.POS, record.REF, record.ALT[0], record.aaf] + list(per_region_maf))

    #if output_record: writer.write_record(record)

vcf.close(); prm.close()#; writer.close()

```

`cyvcf_submit.sh`  and `cyvcf_population_maf.py`  → `per_population_maf/per_population_maf.{chr}.txt.gz`

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

Note: Be careful with these outputs because they are not quite MAF stats— they are the frequencies of the ALT allele, which is not always the minor allele.

**Step 8) Overall stats on autosomal VCF** 

working directory: `../vcf_stats`

`gatk_VCF_metrics.sh` → `filtered_site_0.1_cohort_missing_VCF_metrics.variant_calling_{detail,summary}_metrics` 

This uses dbsnp build 129 to conservatively identify overlapping variants and more accurate Ti/Tv ratio, we also ran with a newer build 156 to see overlap with known variants:

`gatk_VCF_metrics.sh` → `dbsnp156_filtered_site_0.1_cohort_missing_VCF_metrics.variant_calling_{detail,summary}_metrics`

```bash
gatk=/project/lbarreiro/USERS/bridget/programs/gatk-4.2.6.1/gatk
module load java

echo "Collect VCF Metrics"
$gatk CollectVariantCallingMetrics \
  --DBSNP dbsnp_156_autosomes.b37.vcf.gz \
  -I /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/raw_vcfs/baseline_filters_site_0.1_cohort_missing/freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing.vcf.gz \
  -O dbsnp156_filtered_site_0.1_cohort_missing_VCF_metrics \
  --THREAD_COUNT 10
```