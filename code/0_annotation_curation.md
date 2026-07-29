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
