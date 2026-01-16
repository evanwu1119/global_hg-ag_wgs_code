# Gene enrichment resources

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