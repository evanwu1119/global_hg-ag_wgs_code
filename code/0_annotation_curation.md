# Annotation curation

This page details how all of the annotations used for downstream analyses were processed if they were not used as available publicly or obtained from providers. Look out for usage of these output files in the the following code pages. 

## Gene and gene set resources

**Step 1) The set of all gene sets tested was a combination of the most up-to-date GO, KEGG, and Hallmark gene set annotations. GO annotations were obtained by [downloading GMTs from the GOSt browser tool](https://biit.cs.ut.ee/gprofiler/gost), while the others were from querying `bioMart`. These were then concatenated and formatted into appropriate objects for `clusterProfiler` enrichment functions. We tested 22,916 gene sets in total.**

`geneset_annot.R` → `go_kegg_hallmark_term2gene.rds` and `go_kegg_hallmark_term2name.rds`

```r
library(clusterProfiler)
library(data.table)
library(EnrichmentBrowser)
library(tidyverse)

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
```

**Step 2) Gene locations and annotations were obtained from Gencode V46 lifted over to GRCh37. We subsetted for all testable genes: genes in autosomal regions that were contained in at least one tested gene set. We therefore evaluated 23,317 genes and used this as our statistical background for ORA.**

`geneset_test.R` → `hg19_genes.rds`

```r
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(rtracklayer)

# Get testable genes based on all the genesets we are testing and intersect with Gencode V46, we get a total of 21,317 genes as our background
hg19_genes <- import("gs_data/gencode.v46lift37.basic.annotation.gtf.gz")
hg19_genes <- hg19_genes %>% filter(type == "gene" & gene_name %in% unique(go_kegg_hallmark_term2gene$gene)) %>% dplyr::select(gene_id, gene_name)
seqlevelsStyle(hg19_genes) <- "NCBI"
hg19_genes <- hg19_genes %>% filter(seqnames %in% as.character(seq(1:22)) & !duplicated(gene_name)) # remove non-autosomal genes and take the first sequence of duplicated genes
saveRDS(hg19_genes, "hg19_genes.rds")
hg19_genes <- readRDS("hg19_genes.rds")
background <- unique(sort(hg19_genes$gene_name)) # statistical background for GO tests
```

## Retrieving second child GO terms

**Step 1) Some of our enrichment results were too large and redundant due to GO's hierarchical structure, so we wanted to grab more meaningful functions to intrepret. Get term IDs of first children of the root terms (60 terms).**

`get_go_child1.sh` → `GO_child1.txt`

```bash
curl -X GET -H 'Accept:application/json' 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A0008150%2CGO%3A0005575%2CGO%3A0003674/children' | \
	grep -oh "\w*GO:\w*" | \
	grep -v "GO:0008150\|GO:0005575\|GO:0003674" > GO_child1.txt
```

**Step 2) Get term ID and names of second children (1,413 terms), some first children (4 terms) are also second children but this will exclude those as well.**

`get_go_child2.sh` → `GO_child2.txt`

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

## KEGG gene set curation (for metabolism and endocrine gene sets)

**Step 1) Our goal is to extract all KEGG genes from the Metabolism (09100) module and Endocrine (01952) submodule. All pathway information in `KEGG_09100_Metabolism/${subcategory}_pathways.txt` and `endocrine/endocrine_system_reference.txt` is comprised of the codes for pathways within each subcategory. We use the following as examples of getting one pathway from one subcategory:**
- Metabolism (09100)
    - Carbohydrate metabolism (09101) 
        - E.g. for Glycolysis / Gluconeogenesis pathway (ko00010) genes, do [`https://rest.kegg.jp/get/hsa00010`](https://rest.kegg.jp/get/hsa00010)

	- Endocrine system (09152)
		- E.g. for Estrogen signaling pathway (ko04915) genes, do [`https://rest.kegg.jp/get/hsa04915`](https://rest.kegg.jp/get/hsa04915)

**Step 2) Run script to fetch data from KEGG API using `${subcategory}_pathways.txt` files for Metabolism, and the `endocrine_system_reference.txt` file for Endocrine system.** 

`fetch_KEGG_genes.sh ${subcategory}` → `${subcategory}_${pathway}_genes.txt`

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

**Step 3) Concatenate each pathway’s gene list under the broader category (e.g. Carbohydrate metabolism, Energy metabolism).**

`concat_KEGG_genes.sh` → `${category}_ALLGENES.txt`

```bash
for genes_file in *_genes.txt; do
    if [ -f "$genes_file" ]; then
        prefix="${genes_file%_*_*}"  # Extract prefix including the pathway information from the filename
        cat "$genes_file" >> "${prefix}_ALLGENES.txt"
    fi
done

cat *ALLGENES.txt > 09100_METABOLISM_ALLGENES.txt
```

**Step 4) Add the individual genes manually that weren’t grouped under a pathway for subcategories that are too small.**

## cCRE annotations

**Step 1) Liftover annotation file from providers from build hg38 to hg37 using [UCSC LiftOver in-browser](https://genome.ucsc.edu/cgi-bin/hgLiftOver) with `${tissue}-cCREs_GRCh38.bed` as the input. The resulting file `lifted_${tissue}-cCREs_GRCh37.bed` has a few issues that need to be cleaned up.**

**Step 2) Clean up dataframes. Replace names and inputs with other tissue files for same outputs.** 

`i-cCREs_liftover_cleanup.R` → `${tissue}-cCREs_GRCh37.tsv`

```r
#####################################
#####         READ ME           #####
#####################################

# Extracting the lifted-over immune cCREs to redo analysis
# Lifted from GRCh38 to GRCh37

setwd("/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/analysis/iHS")
library(data.table)
library(gtools)

cCREs_37 <- fread("lifted_i-cCREs_GRCh37.bed", header=FALSE)
cCREs_38 <- fread("i-cCREs_GRCh38.bed", header = FALSE)
colnames(cCREs_37) <- c("chr", "start", "end", "internal_ID", "SCREEN_ID")
colnames(cCREs_38) <- c("chr", "start", "end", "internal_ID", "SCREEN_ID")

# Filter rows where SCREEN_ID matches between the two dataframes
matching_rows <- merge(cCREs_37, cCREs_38, by = "SCREEN_ID", suffixes = c("_37", "_38"))

# Identify rows where chr values are different
mismatched_rows <- matching_rows[matching_rows$chr_37 != matching_rows$chr_38, ]

# Identify SCREEN_IDs where chr values are different
mismatched_SCREEN_IDs <- matching_rows$SCREEN_ID[matching_rows$chr_37 != matching_rows$chr_38]

# Retain only rows in cCREs_37 where SCREEN_ID is not in mismatched_SCREEN_IDs
cCREs_37_filtered <- cCREs_37[!cCREs_37$SCREEN_ID %in% mismatched_SCREEN_IDs, ]


# Add back in annotation columns + original header
cCREs_38_annotations <- fread("i-cCREs_GRCh38.tsv")

# Select only the relevant columns from cCREs_38_annotations
annotations_to_merge <- cCREs_38_annotations[, 5:19]

cCREs_37_final <- merge(cCREs_37_filtered, annotations_to_merge,
                        by.x = "SCREEN_ID", by.y = "cCRE",
                        all.x = TRUE)

# Reorder columns and rows
cCREs_37_final$chr <- factor(cCREs_37_final$chr, levels = mixedsort(unique(cCREs_37_final$chr)))
# Order the dataframe by chr and start
cCREs_37_final <- cCREs_37_final[order(cCREs_37_final$chr, cCREs_37_final$start), ]

# Rearrange columns: SCREEN_ID right after internal_ID
col_order <- c(
  "chr", "start", "end", "internal_ID", "SCREEN_ID",
  setdiff(colnames(cCREs_37_final), c("chr", "start", "end", "internal_ID", "SCREEN_ID"))
)

# Reorder the dataframe
cCREs_37_final <- cCREs_37_final[, ..col_order]

# Rename columns
setnames(cCREs_37_final, old = c("internal_ID", "SCREEN_ID"), new = c("rDHS", "cCRE"))

# Save cleaned up object
fwrite(cCREs_37_final, "i-cCREs_GRCh37.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```
