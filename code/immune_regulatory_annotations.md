# Immune regulatory annotations

![image.png](image%2026.png)

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