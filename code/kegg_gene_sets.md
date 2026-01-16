# KEGG gene set curation (for diet + metabolism)

## Background

---

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