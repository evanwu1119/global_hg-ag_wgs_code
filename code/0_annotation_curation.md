# Annotation curation

## Gene and gene set resources

---

The set of all gene sets tested was a combination of the most up-to-date GO, KEGG, and Hallmark gene set annotations. GO annotations were obtained by downloading GMTs while the others were from querying bioMart. These were then concatenated and formatted into appropriate objects for clusterProfiler enrichment functions. We tested 22,916 gene sets in total. 

```r
# Import Hallmark, KEGG, and GO genesets into TERM2GENE and TERM2NAME tables
go_bp_term2gene <- read.gmt("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/pbs/gs_data/hsapiens.GO:BP.name.gmt") # 15472 genesets
go_bp_term2name <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/pbs/gs_data/GO_BP_names.txt",header = F)
colnames(go_bp_term2name) <- c("term", "name")

go_cc_term2gene <- read.gmt("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/pbs/gs_data/hsapiens.GO:CC.name.gmt") # 1998 genesets
go_cc_term2name <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/pbs/gs_data/GO_CC_names.txt",header = F)
colnames(go_cc_term2name) <- c("term", "name")

go_mf_term2gene <- read.gmt("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/pbs/gs_data/hsapiens.GO:MF.name.gmt") # 5037 genesets
go_mf_term2name <- fread("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/pbs/gs_data/GO_MF_names.txt",header = F)
colnames(go_mf_term2name) <- c("term", "name")

kegg_gs <- getGenesets(org = "hsa", gene.id.type = "SYMBOL", db = "kegg", mode = "biomart") # 359 KEGG pathway genesets
kegg_terms <- strsplit(names(kegg_gs), "_")
kegg_ids <- sapply(kegg_terms, function(x) x[1])
kegg_names <- sapply(kegg_terms, function(x) paste0(x[-1], collapse = " "))
kegg_term2name <- data.frame(kegg_ids, kegg_names)
colnames(kegg_term2name) <- c("term", "name")
kegg_term2gene <- lapply(1:length(kegg_ids), function(x) data.frame(kegg_ids[x], kegg_gs[[x]])) %>% bind_rows()
colnames(kegg_term2gene) <- c("term", "gene")

hallmark_gs <- getGenesets(org = "hsa", gene.id.type = "SYMBOL", db = "msigdb", cat = "H") # 50 genesets
hallmark_terms <- strsplit(names(hallmark_gs), "_")
hallmark_ids <- sapply(hallmark_terms, function(x) x[1])
hallmark_names <- sapply(hallmark_terms, function(x) paste0(x[-1], collapse = " "))
hallmark_term2name <- data.frame(hallmark_ids, hallmark_names)
colnames(hallmark_term2name) <- c("term", "name")
hallmark_term2gene <- lapply(1:length(hallmark_ids), function(x) data.frame(hallmark_ids[x], hallmark_gs[[x]])) %>% bind_rows()
colnames(hallmark_term2gene) <- c("term", "gene")

go_kegg_hallmark_term2gene <- rbind(go_bp_term2gene, go_mf_term2gene, go_cc_term2gene, kegg_term2gene, hallmark_term2gene)
go_kegg_hallmark_term2name <- rbind(go_bp_term2name, go_mf_term2name, go_cc_term2name, kegg_term2name, hallmark_term2name) # 22916 genesets tested total
saveRDS(go_kegg_hallmark_term2gene, "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/pbs/gs_data/go_kegg_hallmark_term2gene.rds")
saveRDS(go_kegg_hallmark_term2name, "/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/pbs/gs_data/go_kegg_hallmark_term2name.rds")
go_kegg_hallmark_term2gene <- readRDS("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/pbs/gs_data/go_kegg_hallmark_term2gene.rds")
go_kegg_hallmark_term2name <- readRDS("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/evan/pbs/gs_data/go_kegg_hallmark_term2name.rds")

```

Gene locations and annotations were obtained from Gencode V46 lifted over to GRCh37. We calculated PBS for all testable genes: genes in autosomal regions that were contained in at least one tested gene set. We therefore evaluated 23,317 genes and used this as our statistical background for ORA. 

```r
# Get testable genes based on all the genesets we are testing and intersect with Gencode V46, we get a total of 21,317 genes as our background
hg19_genes <- import("gs_data/gencode.v46lift37.basic.annotation.gtf.gz")
hg19_genes <- hg19_genes %>% filter(type == "gene" & gene_name %in% unique(go_kegg_hallmark_term2gene$gene)) %>% dplyr::select(gene_id, gene_name)
seqlevelsStyle(hg19_genes) <- "NCBI"
hg19_genes <- hg19_genes %>% filter(seqnames %in% as.character(seq(1:22)) & !duplicated(gene_name)) # remove non-autosomal genes and take the first sequence of duplicated genes
saveRDS(hg19_genes, "hg19_genes.rds")
hg19_genes <- readRDS("hg19_genes.rds")
background <- unique(sort(hg19_genes$gene_name)) # statistical background for GO tests
```

Code to get the 100kb around gene centers, negative/greater than seqlength values shouldn’t affect the findOverlaps function (as I tested) so we just take the calculated ranges as-is. 

```bash
hg19_genes <- readRDS("hg19_genes.rds")
starts <- start(hg19_genes)
ends <- end(hg19_genes)

# Calculate windows
centers <- round((starts + ends) / 2)
window_size <- 100000
window_starts <- centers - (window_size / 2)
window_ends <- centers + (window_size / 2)
windows <- IRanges(window_starts, window_ends)

# Save object
hg19_genes_100kb <- hg19_genes
ranges(hg19_genes_100kb) <- windows
saveRDS(hg19_genes_100kb, "hg19_genes_100kb.rds")
```

## Subsetting second child GO terms

---

GO API commands: [https://www.ebi.ac.uk/QuickGO/api/index.html#!/gene_ontology/findChildrenUsingGET_1](https://www.ebi.ac.uk/QuickGO/api/index.html#!/gene_ontology/findChildrenUsingGET_1) 

Goal is to extract terms 2 children down from the root terms:

- GO:0008150 biological process
- GO:0005575 cellular component
- GO:0003674 molecular function

**Step 1) Get term IDs of first children of the root terms (60 terms)**

```bash
curl -X GET -H 'Accept:application/json' 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A0008150%2CGO%3A0005575%2CGO%3A0003674/children' | \
	grep -oh "\w*GO:\w*" | \
	grep -v "GO:0008150\|GO:0005575\|GO:0003674" > GO_child1.txt
```

**Step 2) Get term ID and names of second children (1,413 terms), some first children (4 terms) are also second children but this will exclude those as well**

```bash
QUERY=$(awk '(ORS="%2C"){sub(":","%3A",$1)}1' GO_child1.txt)
IGNORE=$(awk '(FS="\\|")1' GO_child1.txt)

curl -X GET -H 'Accept:application/json' "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/${QUERY}/children" | \
	tee >(grep -oP '(?<=\"id\":)\".*?,' > tmp1.txt) | grep -oP '(?<=\"name\":)\".*?\"' > tmp2.txt
	
paste -d "" tmp1.txt tmp2.txt | \
	sort | uniq | \
	grep -v "${IGNORE}" > GO_child2.csv 
	
rm tmp1.txt tmp2.txt
```

# KEGG gene set curation (for diet + metabolism)

## Background

- Human pathway map for reference: [https://www.genome.jp/kegg-bin/show_organism?menu_type=genome_info&org=hsa](https://www.genome.jp/kegg-bin/show_organism?menu_type=genome_info&org=hsa)
- **KEGG MODULE**
    - The **KEGG MODULE** database consists of KEGG modules identified by M numbers and KEGG reaction modules identified by RM numbers, which are manually defined functional units of gene sets and reaction sets, respectively. KEGG modules are further divided into pathway modules and signature modules as shown below.
    - pathway modules – functional units of gene sets in metabolic pathways, including molecular complexes
    - signature modules – functional units of gene sets that characterize phenotypic features
    - reaction modules – functional units of successive reaction steps in metabolic pathways
    - KEGG modules are associated with dynamically generated graphical diagrams. For example, [M00002](https://www.genome.jp/module/M00002) represents glycolysis core module involving three-carbon compounds and its organism specific module takes the form of [hsa_M00002](https://www.genome.jp/module/hsa_M00002)

***Goal: Extract all KEGG genes from the following Metabolism (09100) module***

- One level down is:
    - Carbohydrate metabolism
        - E.g. for Glycolysis / Gluconeogenesis pathway genes, do [`https://rest.kegg.jp/get/hsa00010`](https://rest.kegg.jp/get/hsa00010)
        - If you do [`https://rest.kegg.jp/link/hsa/hsa00010`](https://rest.kegg.jp/link/hsa/hsa00010) you will just fetch a tab-delimited file with two columns (1st column indicates path, 2nd column indicates the gene with the hsa IDs, e.g. hsa:10327)
    - Energy metabolism
    - Lipid metabolism
    - Nucleotide metabolism
    - Amino acid metabolism
    - Metabolism of other amino acids
    - Glycan biosynthesis and metabolism
    - Metabolism of cofactors and vitamins
    - Metabolism of terpenoids and polyketides
    - Biosynthesis of other secondary metabolites
    - Xenobiotics biodegradation and metabolism

![Screenshot 2024-03-19 at 16.54.45.png](Screenshot_2024-03-19_at_16.54.45.png)

## List of metabolic pathways to extract genes from and their KEGG codes (hsaXXXX)

---

***09101 Carbohydrate metabolism***

```
hsa00010
hsa00020
hsa00030
hsa00040
hsa00051
hsa00052
hsa00053
hsa00500
hsa00520
hsa00620
hsa00630
hsa00640
hsa00650
hsa00562
```

***09102 Energy metabolism***

```bash
hsa00190
hsa00910
hsa00920
```

***09103 Lipid metabolism***

```bash
hsa00061
hsa00062
hsa00071
hsa00100
hsa00120
hsa00140
hsa00561
hsa00564
hsa00565
hsa00600
hsa00590
hsa00591
hsa00592
hsa01040
```

***09104 Nucleotide metabolism***

```bash
hsa00230
hsa00240
```

***09105 Amino acid metabolism***

```bash
hsa00250
hsa00260
hsa00270
hsa00280
hsa00290
hsa00310
hsa00220
hsa00330
hsa00340
hsa00350
hsa00360
hsa00380
hsa00400
```

***09106 Metabolism of other amino acids***

```bash
hsa00410
hsa00430
hsa00440
hsa00450
hsa00470
hsa00480
```

***09107 Glycan biosynthesis and metabolism***

```bash
hsa00510
hsa00513
hsa00512
hsa00515
hsa00514
hsa00532
hsa00534
hsa00533
hsa00531
hsa00563
hsa00601
hsa00603
hsa00604
hsa00511
```

***09108 Metabolism of cofactors and vitamins***

```bash
hsa00730
hsa00740
hsa00750
hsa00760
hsa00770
hsa00780
hsa00785
hsa00790
hsa00670
hsa00830
hsa00860
hsa00130
```

***09109 Metabolism of terpenoids and polyketides***

```
hsa00900

Add the following genes (sub-categories were so small there wasn't a pathway):

7086 TKT; transketolase
84076 TKTL2; transketolase like 2
8277 TKTL1; transketolase like 1
23483 TGDS; TDP-glucose 4,6-dehydratase
23483 TGDS; TDP-glucose 4,6-dehydratase
```

***09110 Biosynthesis of other secondary metabolites***

```
hsa00232
hsa00524

Add the following genes (sub-categories were so small there wasn't a pathway):

1610 DAO; D-amino acid oxidase
9060 PAPSS2; 3'-phosphoadenosine 5'-phosphosulfate synthase 2
9061 PAPSS1; 3'-phosphoadenosine 5'-phosphosulfate synthase 1
3101 HK3; hexokinase 3
3098 HK1; hexokinase 1
3099 HK2; hexokinase 2
80201 HKDC1; hexokinase domain containing 1
2645 GCK; glucokinase
27349 MCAT; malonyl-CoA-acyl carrier protein transacylase
31 ACACA; acetyl-CoA carboxylase alpha
32 ACACB; acetyl-CoA carboxylase beta
```

***09111 Xenobiotics biodegradation and metabolism***

```
hsa00980
hsa00982
hsa00983
```

Note: All of the above information is stored in `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/curate_genes/KEGG_09100_Metabolism` Each *_pathways.txt file is comprised of the codes for pathways within each category. This is used to fetch all the gene-level information from the API.

**Step 2: Run script to fetch data from KEGG API** `fetch_KEGG_genes.sh`

- Remember to do `chmod +x fetch_KEGG_genes.sh`
- This script will iterate over all **`_pathways.txt`** files in the current directory, extract the prefix from each filename, fetch the data for each ID in the file, and save the results to separate files with the pattern **`${prefix}_${id}_genes.txt`**, waiting 1 second between each query as per the API's rate limit

```bash
#!/bin/bash

# Function to fetch data for each ID
fetch_data() {
    local pathway_file="$1"
    local prefix="${pathway_file%_*}"  # Extract prefix from the filename
    while IFS= read -r id; do
        curl -s "https://rest.kegg.jp/link/hsa/$id" > "${prefix}_${id}_genes.txt"
        sleep 1  # Wait for 1 second
    done < "$pathway_file"
}

# Iterate over each _pathways.txt file in the current directory
for pathway_file in *_pathways.txt; do
    if [ -f "$pathway_file" ]; then
        fetch_data "$pathway_file"
    fi
done
```

**Step 3: Concatenate each pathway’s gene list under the broader category (e.g. Carbohydrate metabolism, Energy metabolism)**

```bash
for genes_file in *_genes.txt; do
    if [ -f "$genes_file" ]; then
        prefix="${genes_file%_*_*}"  # Extract prefix including the pathway information from the filename
        cat "$genes_file" >> "${prefix}_ALLGENES.txt"
    fi
done
```

**Step 4: Concatenate all of them into a broader metabolism gene list**

```r
cat *ALLGENES.txt > 09100_METABOLISM_ALLGENES.txt
```

**Step 5: Add the individual genes that weren’t grouped under a pathway for *09109 Metabolism of terpenoids and polyketides and 09110 Biosynthesis of other secondary metabolites (see above)***

Once I have the list of hsa gene IDs (which are NCBI IDs), I can then convert to ENSEMBL IDs

- The hsa gene ID to NCBI gene ID conversion table is downloaded [here](https://rest.kegg.jp/conv/ncbi-geneid/hsa)

# Immune regulatory annotations

> Many uses of cCREs are based on the regulatory role associated with their biochemical signatures. Thus, we putatively defined cCREs in one of the following annotation groups based on each element’s dominant biochemical signals across all available biosamples. Analogous to GENCODE's catalog of genes, which are defined irrespective of their varying expression levels and alternative transcripts across different cell types, we provide a general, cell type-agnostic classification of cCREs based on the max-Zs as well as its proximity to the nearest annotated TSS:
> 
> 1. *cCREs with promoter-like signatures (cCRE-PLS)* fall within 200 bp (center to center) of an annotated GENCODE TSS and have high DNase and H3K4me3 signals (evaluated as DNase and H3K4me3 max-Z scores, defined as the maximal DNase or H3K4me3 Z scores across all biosamples with data; see Methods).
> 2. *cCREs with enhancer-like signatures (cCRE-ELS)* have high DNase and H3K27ac max-Z scores and must additionally have a low H3K4me3 max-Z score if they are within 200 bp of an annotated TSS. The subset of cCREs-ELS within 2 kb of a TSS is denoted proximal (cCRE-pELS), while the remaining subset is denoted distal (cCRE-dELS).
> 3. *DNase-H3K4me3 cCREs* have high H3K4me3 max-Z scores but low H3K27ac max-Z scores and do not fall within 200 bp of a TSS.
> 4. *CTCF-only cCREs* have high DNase and CTCF max-Z scores and low H3K4me3 and H3K27ac max-Z scores.

## Notes

---

- **Working directory:** `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/functional_annotations`
- **Annotations from Zhipeng Weng:**
    
    `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/i-cCREs.tsv`
    

## Methods

---

- Liftover annotations from Weng Lab from build 38 to 37
- Set up simple logistic regression model using a single predictor variable (the cCRE in question), binarizing iHS scores ( > 99th percentile in the population == “high iHS”)
    - e.g. high_iHS ~ proximal enhancer-like signature

## Liftover from GRCh38 to GRCh37

---

**Step 1) Liftover .bed file from build 38 to 37**

- use [https://genome.ucsc.edu/cgi-bin/hgLiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) with `i-cCREs_GRCh38.bed` as the input
- the resulting file `lifted_i-cCREs_GRCh37.bed` has a few issues that need to be cleaned up —
    - some ranges converted from one chromosome to another, those just need to be dropped
    - it doesn’t have the header and info columns re. which cell-types the cCREs are present in

**Step 2) Clean dataframes up with** `i-cCREs_liftover_cleanup.R`

**Important:** Final cleaned up .tsv file with proper annotations and GRCh37 coordinates are in 
`i-cCREs_GRCh37.tsv`

## Results from 100kb analysis

---

![image.png](image%2027.png)

![image.png](image%2028.png)

![image.png](image%2029.png)

![image.png](image%2030.png)

![image.png](image%2031.png)

- 100kb iHS scores are in `freeze2_iHS_per_window_100kb_df.RData`
- Script for plotting + running logistic regression models is `immune_cCREs_100kb.R` and workspace is saved in `immune_cCREs_100kb.RData`
    - The object `filtered_combined_results` has the pvals for populations > 15 individuals for each element type

### *Getting a combined p-value for each element type*

```r
filtered_combined_results %>%
  group_by(element_type) %>%
  summarise(Combined_P_Value = pchisq(-2 * sum(log(P_Value)), df = 2 * n(), lower.tail = FALSE))
  
 
 ## A tibble: 4 × 2
  element_type Combined_P_Value
  <chr>                   <dbl>
1 CA_H3K4me3           6.28e- 1
2 PLS                  3.27e-11
3 dELS                 2.88e- 2
4 pELS                 2.23e-18
```

![image.png](image%2032.png)

## Results from 50kb analysis

---

![image.png](image%2033.png)

![image.png](image%2034.png)

![image.png](image%2035.png)

![image.png](image%2036.png)

![image.png](image%2037.png)

```bash
# Number of windows where each cCRE appears
pELS         CA    CA-CTCF CA-H3K4me3       dELS        PLS      CA-TF        TF 
18114       8741       8643       5507      38034      14141       2977      246
```

- 50kb iHS scores are in `freeze2_iHS_per_window_50kb_df.RData`
- Script for plotting + running logistic regression models is `immune_cCREs_50kb.R` and workspace is saved in `immune_cCREs_50kb.RData`

![image.png](image%2038.png)

## Overall Takeaways from 10, 50, 100kb windowed analysis

---

- Of the 4 tested cCREs, pELS, dELS, and PLS presence seem to be predictive (to varying degrees) of high iHS scores (i.e. show evidence of selection)
- The statistical significance of modeled odds ratios seem correlated with the number of elements that exist + the number of windows (observations) we can test

## Cell-type-specific analysis on 100kb windows

---

Refer to `immune_cCREs_100kb_CTanalysis_streamlined.R`  for script that annotates, models, and plots each cell type condition (-specific vs. -present) 

![image.png](image%2039.png)

```
# for PLS
1: glm.fit: fitted probabilities numerically 0 or 1 occurred
2: glm.fit: fitted probabilities numerically 0 or 1 occurred
3: glm.fit: fitted probabilities numerically 0 or 1 occurred
4: glm.fit: fitted probabilities numerically 0 or 1 occurred

# This is happening for the PLS and H3K models, where lots of probabilities are close to 0, probably because we have very sparse data for those elements 
```

All cell-type specific results (logistic regression odds ratios and plots) are in 
`/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/functional_annotations/immune_100kb_CTspecific_results`

![image.png](image%2040.png)

![class_level_FE_meta_analysis_CA_H3K4me3.png](class_level_FE_meta_analysis_CA_H3K4me3.png)

![class_level_FE_meta_analysis_dELS.png](class_level_FE_meta_analysis_dELS.png)

![class_level_FE_meta_analysis_PLS.png](class_level_FE_meta_analysis_PLS.png)

![class_level_FE_meta_analysis_pELS.png](class_level_FE_meta_analysis_pELS.png)

## Fixed effect Meta-analysis

---

- Analysis uses `library(metafor)` package
- Refer to this paper for reference on fixed vs. random effects meta-analysis
    
    [fixed vs random effect models.pdf](fixed_vs_random_effect_models.pdf)
    

***Immune analysis***

- `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/functional_annotations/immune_cCREs_100kb_fixed_effect_meta.R` and
- `immune_cCREs_100kb_fixed_effect_meta.RData`
- Detailed figures (expand to view)
    
    [immune-cCREs_100kb_fixed_effect_meta.pdf](immune-cCREs_100kb_fixed_effect_meta.pdf)
    

***Brain analysis***

- `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS/functional_annotations/brain_cCREs_100kb_fixed_effect_meta.R`
- `brain_cCREs_100kb_fixed_effect_meta.RData`
- Detailed figures (expand to view)
    
    ![image.png](image%2041.png)
    
    ![image.png](image%2042.png)
    
    ![image.png](image%2043.png)
    
    ![image.png](image%2044.png)
    

***Comparison between estimates for immune-cCREs and brain-cCREs***

Analysis is done in `brain_cCREs_100kb_fixed_effect_meta.R`

![image.png](image%2045.png)

![image.png](image%2046.png)

[Liftover annotations from Weng Lab from build 38 to 37 (see the initial iHS analysis page for details)](Liftover%20annotations%20from%20Weng%20Lab%20from%20build%2038%20t%201ad749e8fe65805ba514d7033220b7bd.md)
