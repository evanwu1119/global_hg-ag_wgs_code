# Global PCAs

Using `baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf` filters

**Step 1) Generate PLINK fileset using VCF** `baseline_filters_site_0.1_cohort_missing_nohet_nokin`

```bash
plink --vcf freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.vcf.gz --make-bed --out freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin

918 people (0 males, 0 females, 918 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 918 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.980674.
62908805 variants and 918 people pass filters and QC.

mv freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin.* /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/baseline_filters_site_0.1_cohort_missing_nohet_nokin
```

**Step 2) Update PLINK famID to population labels**

```r
setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/baseline_filters_site_0.1_cohort_missing_nohet_nokin")

ref <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v2.txt")
fam <- read.table("freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin.fam")
fam$V1 <- ref$pop_code[match(fam$V1, ref$study_ID)]
fam$V6 <- 1

write.table(fam, file = "freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin.fam", col.names =F, row.names=F, quote=F)
```

**Step 3) Add MAF filter (1%)** 

```bash
plink --bfile freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin --maf 0.01 --make-bed --out freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf

62908805 variants loaded from .bim file.
918 people (0 males, 0 females, 918 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf.nosex .
918 phenotype values loaded from .fam.
Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
phenotypes to be ignored, use the --allow-no-sex flag.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 918 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.980674.
48718164 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
14190641 variants and 918 people pass filters and QC.
Among remaining phenotypes, 0 are cases and 918 are controls.

plink --bfile freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf --pca --out freeze2_baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf
```

**Step 3) Plot PCA with `freeze2_pcas.R` in `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/pca/`**

Do this using `baseline_filters_site_0.1_cohort_missing_nohet_nokin_0.01_maf` filters

# Global PCAs

### 12/22/24 Updated after re-filtering VCF to include polymorphic variants only

---

Working directory: `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/updated`

**Step 1) Generate PLINK fileset using VCF** `baseline_filters_site_0.1_cohort_missing_nohet_nokin`

```bash
plink --vcf freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.vcf.gz --make-bed --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin

57349118 variants loaded from .bim file.
918 people (0 males, 0 females, 918 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.nosex
.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 918 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.98089.
57349118 variants and 918 people pass filters and QC.

mv freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.* /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/updated/baseline_filters_site_0.1_cohort_missing_nohet_nokin
```

**Step 2) Update PLINK famID to population labels**

```r
setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/plink/updated/baseline_filters_site_0.1_cohort_missing_nohet_nokin")

ref <- read.delim("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/metadata/freeze2_metadata_v2.txt")
fam <- read.table("freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.fam")
fam$V1 <- ref$pop_code[match(fam$V1, ref$study_ID)]
fam$V6 <- 1

write.table(fam, file = "freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin.fam", col.names =F, row.names=F, quote=F)
```

**Step 3) Add MAF filter (1%) and run PCA**

```bash
plink --bfile freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin --maf 0.01 --make-bed --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf

57349118 variants loaded from .bim file.
918 people (0 males, 0 females, 918 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf.nosex
.
918 phenotype values loaded from .fam.
Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
phenotypes to be ignored, use the --allow-no-sex flag.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 918 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.98089.
43424966 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
13924152 variants and 918 people pass filters and QC.

plink --bfile freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf --pca --out freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing_nohet_nokin_0.01_maf

cp *eigen* /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/pca/
```

**Step 4) Plot PCA with `freeze2_pcas_updated.R` in `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/pca/`**