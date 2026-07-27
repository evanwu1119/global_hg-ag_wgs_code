# 2) Phasing and variant annotation

blah

## Variant annotation

**Step 1) Generate annotations using Annovar per-chromosome to get gene region, SNP consequence, clinical significance, etc. annotations**

 `run_table_annovar_per_chr.sh` → `chr${chr}_anno.hg19_multianno.txt`

```bash
#!/bin/bash

table_annovar.pl /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing/per_chr/freeze2_biallelic_filtered_site_0.1_cohort_missing.${chr}.vcf.gz humandb/ -buildver hg19 -out chr${chr}_anno -remove -protocol refGene,cytoBand,exac03,avsnp150,clinvar_20221231 -operation g,r,f,f,f -vcfinput -polish --thread 6
```

**Step 2) Parse annotations to create a summary**

```bash
#!/bin/bash

# To create a summary for every chromosome, Loop through all files matching the pattern "chr*_multianno.txt"

for file in chr*_multianno.txt; do
    output="${file/_multianno.txt/_anno_summary.txt}"
    awk -F'\t' 'NR>1 {print $6}' "$file" | sort | uniq -c > "$output"
done
```

**Step 2) Collapse some categories of the per-chr summaries**

- When you see exonic;splicing, that means that variants both overlapping a coding region + within 2-bp of a splicing junction
- In the above example, we’d want to add 3 entries to both the exonic and splicing categories respectively, and remove the exonic;splicing category

```bash
for file in chr*anno_summary.txt; do
    # Initialize associative array
    declare -A counts

    while read -r line; do
        # Extract the count and the category/categories from the line
        count=$(echo "$line" | awk '{print $1}')
        category=$(echo "$line" | awk '{$1=""; print $0}' | sed 's/^ //')

        if [[ $category == *";"* ]]; then
            # If the category contains ";", split it and update counts for each
            IFS=';' read -ra split_categories <<< "$category"
            for split_category in "${split_categories[@]}"; do
                counts[$split_category]=$((counts[$split_category]+count))
            done
        else
            # Otherwise, just update the count for this category
            counts[$category]=$((counts[$category]+count))
        fi
    done < "$file"

    output="${file/_summary.txt/_summary_clean.txt}"

    # Write the counts to the new file
    for category in "${!counts[@]}"; do
        echo "${counts[$category]} $category" >> "$output"
    done

    unset counts
done
```

