# Variant Annotation

Use ANNOVAR: [https://annovar.openbioinformatics.org/en/latest/user-guide/startup/](https://annovar.openbioinformatics.org/en/latest/user-guide/startup/)

Use this vcf: `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing/freeze2_autosomes_biallelic_filtered_site_0.1_cohort_missing.vcf.gz`

## Generating annotations by chromosome

---

The file conversion step is very slow so don’t try to do this using the whole genome vcf

Script is `/project/lbarreiro/USERS/bridget/programs/annovar/run_table_annovar_per_chr.sh`

```bash
for chr in {1..22}; do export chr; sbatch run_table_annovar_per_chr.sh ; done

#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --account=pi-lbarreiro
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=70G
#SBATCH --ntasks=6
#SBATCH --job-name=annovar

table_annovar.pl /project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing/per_chr/freeze2_biallelic_filtered_site_0.1_cohort_missing.${chr}.vcf.gz humandb/ -buildver hg19 -out chr${chr}_anno -remove -protocol refGene,cytoBand,exac03,avsnp150,clinvar_20221231 -operation g,r,f,f,f -vcfinput -polish --thread 6
```

All outputs are in `/project/lbarreiro/USERS/bridget/huntergatherer/freeze2/vcf/baseline_filters_site_0.1_cohort_missing/annotation`

## Parsing annotations

---

```bash
# To get a count summary of how many entries are in each category under column X of the *multianno.txt outputs, do the following. Example below is for column 6, variant type. 
awk -F'\t' 'NR>1 {print $6}' chr22_anno.hg19_multianno.txt | sort | uniq -c

   9571 downstream
  15988 exonic
      3 exonic;splicing
 319693 intergenic
 378255 intronic
   7753 ncRNA_exonic
      6 ncRNA_exonic;splicing
  76207 ncRNA_intronic
     42 ncRNA_splicing
     78 splicing
   9164 upstream
    497 upstream;downstream
  15510 UTR3
   3101 UTR5
     19 UTR5;UTR3
     
# To create a summary for every chromosome, Loop through all files matching the pattern "chr*_multianno.txt"

for file in chr*_multianno.txt; do
    output="${file/_multianno.txt/_anno_summary.txt}"
    awk -F'\t' 'NR>1 {print $6}' "$file" | sort | uniq -c > "$output"
done
```

**Step 2) Collapse some categories of the per-chr summaries**

Reference: [https://annovar.openbioinformatics.org/en/latest/user-guide/gene/#output-file-1-refseq-gene-annotation](https://annovar.openbioinformatics.org/en/latest/user-guide/gene/#output-file-1-refseq-gene-annotation)

- When you see exonic;splicing, that means that variants both overlapping a coding region + within 2-bp of a splicing junction
- In the above example, we’d want to add 3 entires to both the exonic and splicing categories respectively, and remove the exonic;splicing category

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

The all-chrs output is here `all_autosomes_hg19_anno_summary_clean.txt`

```bash
650868 exonic
695568 UTR3
137059 UTR5
4022414 ncRNA_intronic
256911 ncRNA_exonic
3811 splicing
438547 downstream
413931 upstream
32415182 intergenic
1438 ncRNA_splicing
23890795 intronic
```

**Percentage breakdown of variant type**

![Screenshot 2024-04-09 at 18.16.32.png](Screenshot_2024-04-09_at_18.16.32.png)