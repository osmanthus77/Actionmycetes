# 对antismash结果进行分析

```shell

# filter into each product one C domain DNA sequence file
for product in corbomycin GP6738 rimomycin misaugamycin; do
    cat domain_Cdna.txt | grep -F -A 1 -f ${product}.txt | sed '/^--$/d' > domain_Cdna_${product}.txt
    wc -l domain_Cdna_${product}.txt
done


```

```shell
# extract all domains from all antismash results of 1559 complete genomes
cd ~/project/Actionmycetes/antismash/antismash_summary

cat aa2/aa_all_1559_filter.tsv | grep -v '^$' |
    while IFS=$'\t' read -r strain cluster family level product _ _ _; do
        echo "${strain}\t${cluster}"
        js=$(ls ../antismash_result/${family}/${level}/${strain}/regions.js);
        type=$(python aa_antismash_pp_for_complete_uncomplete.py "${js}" "${cluster}" "${strain}");
        echo ${type} | sed "s/]/]\n/g" >> aa2/domain_aa/domain_aa_all_1559.txt;
    done

# correct the format
cat aa2/domain_aa/domain_aa_all_1559.txt |
    sed "s/\[//g" | sed "s/\]//g"| sed "s/'//g" | sed "s/,/\n/g" |
    sed "s/ >/>/g" | sed "s/+/\t/g" |
    sed "s/\t/\n/g" |
    sed '/^$/d' \
    >> aa2/domain_aa/ok_domain_aa_all_1559.txt

# filter the C domain
cat aa2/domain_aa/ok_domain_aa_all_1559.txt |
        grep -A 1 -E "Condensation|Cglyc" |
        sed '/^--$/d' \
        >> aa2/domain_aa/domain_Caa_1559.txt

cat aa2/domain_aa/domain_Caa_1559.txt | grep '>' | sed 's/>//g' >> aa2/domain_aa/domain_Caa_1559_annotate.tsv

```






## 比对构树

```shell
# building a tree

cd ~/project/Actionmycetes/
mkdir -p analysis2/sequence_dna analysis2/sequence_aa analysis2/msa analysis2/trim analysis2/fasttree
cd analysis

# Multiple Sequence Alignment
for file in sequence_dna/*; do
    seq=$(basename "$file" .txt)
    echo sequence_dna/${seq}
    mafft --auto sequence_dna/${seq}.txt > msa/${seq}.aln.fa
done
    
# trim
for file in sequence_dna/*; do
    seq=$(basename "$file" .txt)
    trimal -in msa/${seq}.aln.fa -out trim/${seq}.trim.fa -automated1
done

# fasttree v2.2 installation in macOS
# cd ~/biosoft/fasttree2
# wget https://morgannprice.github.io/fasttree/FastTree.c
# gcc -O3 -fopenmp-simd -funsafe-math-optimizations -march=native -o FastTree FastTree.c -lm
# echo 'export PATH=$HOME/biosoft/fasttree2/FastTree:$PATH' >> ~/.bashrc
# source ~/.bashrc

# tree
for file in sequence_dna/*; do
    seq=$(basename "$file" .txt)
    echo ${seq}
    FastTree -nt trim/${seq}.trim.fa > fasttree/${seq}.nwk
done


```

```shell
# to visualize

for file in sequence_dna/*; do
    seq=$(basename "$file" .txt)
    cat $file |
    grep ">" |
    sed "s/^>//g" > annotate/${seq}_annotate.txt
done

# then annotate the sequence file by hand
# strain	Location	Gene	AminoAcid
# Streptomy_albon_ATCC_12461_GCF_008704395_1-Condensation_Starter-2-304-5.1-active	Location1	gene1	Hpg
# Streptomy_albon_ATCC_12461_GCF_008704395_1-Condensation_LCL-1062-1365-5.4-active	Location2	gene1	D-Dpg
# Streptomy_albon_ATCC_12461_GCF_008704395_1-Condensation_LCL-2135-2440-5.7-active	Location3	gene1	D-Trp
# Streptomy_albon_ATCC_12461_GCF_008704395_1-Cglyc-19-317-6.1-active	Location4	gene2	D-Hpg

# annotate with iTOL script
script=fasttree/table2itol/table2itol.R
Rscript ${script} -a -D fasttree/domain_Cdna3 -i Strain -l AminoAcid -w 0.5 annotate/domain_Cdna_annotate2.txt
Rscript ${script} -a -D fasttree/domain_Cdna3 -C fasttree/domain_Cdna3/colours.yml -i Strain -l AminoAcid -w 0.5 annotate/domain_Cdna_annotate3.txt
Rscript ${script} -a -D fasttree/domain_Cdna3 -i Strain -l Cdomain5 -w 0.5 annotate/domain_Cdna_annotate2.txt
Rscript ${script} -a -D fasttree/domain_Cdna3 -i Strain -l GPA -w 0.5 annotate/domain_Cdna_annotate2.txt
Rscript ${script} -a -D fasttree/domain_Cdna3 -i Strain -l Type1 -w 0.5 annotate/domain_Cdna_annotate2.txt
Rscript ${script} -a -D fasttree/domain_Cdna3 -i Strain -l Type2 -w 0.5 annotate/domain_Cdna_annotate2.txt

Rscript ${script} -a -D fasttree/domain_Cdna3 -i Strain -l Type -w 0.5 annotate/domain_Cdna_annotate2.txt


for name in corbomycin GP6738 misaugamycin rimomycin; do
    for circle in Location Gene AminoAcid; do
        script=fasttree/table2itol/table2itol.R
        Rscript ${script} -a -D fasttree/${name}2 -i strain -l ${circle} -w 0.5 sequence_dna/2domain_Cdna_${name}_annotate.txt
    done
done

```



```shell
# 提取 NPRS DNA序列

cd ~/project/Actionmycetes/antismash/antismash_summary

for level in complete_taxon complete_untaxon; do
        product_file="product/cluster_gpa_${level}.tsv"

        identifier_file="aa1/nrps_dna/nrps_dna_identifier.tsv"

        dna_file="aa1/nrps_dna/nrps_dna.tsv"

        cat "$product_file"| grep -v '^$' |
        while IFS=$'\t' read -r strain cluster product _ _ _ _; do
            region=$(echo ${cluster} | sed -E 's/c[0-9]+//g' | sed -E 's/r//g')
            clu=$(echo ${cluster} | sed -E 's/r[0-9]+//g' | sed -E 's/c//g')
			json_path=$(ls product/antismash_${level}/${strain}/*genomic.json 2>/dev/null | head -n 1)
            cds=$(python antimash_cds.py ${json_path} ${region} ${clu})
            identifier=$(echo $cds | sed -E "s/\[//g "| sed "s/\]//g" | sed "s/'//g"  | sed 's/, /\n/g')
            echo ${identifier} | sed "s#^#${strain}\tr${region}c${clu}\t${product}\t#g" \
            >> "$identifier_file"
        done
done

# cds格式
# ['SXIN_RS02445', 'SXIN_RS02440', 'SXIN_RS02435', 'SXIN_RS02430']

for level in complete_taxon complete_untaxon; do
    identifier_file="aa1/nrps_dna/nrps_dna_identifier.tsv"

    dna_file="aa1/nrps_dna/nrps_gene_dna.txt"

    cat "$identifier_file" | grep -v '^$' |
    while IFS=$'\t' read -r strain cluster product identifier _ _; do
        js_file=product/antismash_${level}/${strain}/regions.js
        python antismash_js_cds.py ${js_file} ${cluster} ${identifier} ${strain} \
        >> "$dna_file"
    done
done


# numbers of sequence in one strain merge into one sequence
awk '
/^>/ {
    split(substr($0, 2), a, "+")
    h = a[1]
    if (!(h in seen)) {
        seen[h] = 1
        order[++count] = h
    }
    current = h
    next
}
{
    seq[current] = seq[current] $0
}
END {
    for (i = 1; i <= count; i++) {
        h = order[i]
        print ">" h
        print seq[h]
    }
}
' nrps_gene_dna.txt >> nrps_gene.txt

```

```shell
# split sequence into GPA one file
cd ~/project/Actionmycetes/antismash/antismash_summary

for i in 1-4 corbomycin GP6738 misaugamycin rimomycin; do
    cat domain_C_GPA.txt | grep "${i}" | cut -f 1 > C${i}.txt
done

for i in 1-4 corbomycin GP6738 misaugamycin rimomycin; do
    cat domain_Caa.txt | grep -A 1 -F -i -f C${i}.txt | sed '/^--$/d' > domain_Caa_${i}.txt
done

for i in 1-4 corbomycin GP6738 misaugamycin rimomycin; do
    mafft --auto sequence_aa/domain_Caa_${i}.txt > msa/domain_Caa_${i}.aln.fa
    trimal -in msa/domain_Caa_${i}.aln.fa -out trim/domain_Caa_${i}.trim.fa -automated1
    FastTree trim/domain_Caa_${i}.trim.fa > fasttree/domain_Caa_${i}.nwk
done

```

```shell
# extract P450 locus_tag and translation from region.js
cd ~/project/Actionmycetes/antismash/antismash_summary

for level in complete_taxon complete_untaxon; do
    product_file="product/cluster_gpa_${level}.tsv"

    P450_file="aa1/p450_aa/p450_gene_aa_raw.txt"

    cat "$product_file" | grep -v '^$' |
    while IFS=$'\t' read -r strain cluster product _ _ _ _; do
        echo "${strain}\t${cluster}"
        js_file=product/antismash_${level}/${strain}/regions.js
        python antismash_js_p450.py ${js_file} ${cluster} ${strain} \
        >> "$P450_file"
    done
done

cat aa1/p450_aa/p450_gene_aa.txt | grep '>' | sed 's/>//g' | sed -E "s/-/\t/g; s/\+/\t/g" >> aa1/p450_aa/p450_identifier_raw.tsv
```

```shell
# extract BGC dna sequence from .gbk file
cd ~/project/Actionmycetes/antismash/antismash_summary

for level in complete_taxon complete_untaxon; do
    product_file="product/cluster_gpa_${level}.tsv"
    bgc_file="aa1/bgc_dna/bgc_dna.txt"
    
    cat "$product_file" | grep -v '^$' |
    while IFS=$'\t' read -r strain cluster product _ _ _; do
        clu=$(echo $cluster | awk 'match($0, /r[0-9]*c([0-9]+)/) {
            num = substr($0, RSTART + 3, RLENGTH - 3);
            printf "region%03d\n", num;
            }')
        echo "${strain}_${clu}"
        cluster_gbk_file=$(ls product/antismash_${level}/${strain}/*${clu}.gbk | head -n 1)
        echo "${cluster_gbk_file}"
        sequence=$(awk '/^ORIGIN/{flag=1} flag' ${cluster_gbk_file} | 
            sed -E "s/[0-9]+//g" | sed "s/ *//g" | tr -d '\n' | sed 's/ORIGIN//g' | sed "s/\/\///g")
        echo ">${strain}_${cluster}\n${sequence}" >> ${bgc_file}
    done
done
# There are 'three' regions need to replace correct sequence by hand

awk '/^ORIGIN/{flag=1} flag' Streptomyces_sp_DSM_11171_feglymycin_bgc.gbk | sed -E "s/[0-9]+//g" | sed "s/ *//g" | tr -d '\n' | sed "s/ORIGIN/>Streptomyces_sp_DSM_11171_feglymycin\n/" >> Streptomyces_sp_DSM_11171_feglymycin.txt

```



## BGC clinker

```shell
# extract BGC and flanking sequence about 10k

cd ~/project/Actionmycetes/antismash/antismash_summary

# extract BGC location in genome from html file
for level in complete_taxon complete_untaxon; do
    product_file="product/cluster_gpa_${level}.tsv"
    output_file="aa1/bgc_flanking_10k/bgc_10k_locus_location.tsv"
    
    cat "$product_file" | grep -v '^$' |
    while IFS=$'\t' read -r strain cluster product _ _ _; do
        locus=$(cat product/antismash_${level}/${strain}/index.html |
        pup "div.page#${cluster} div.description-container div.heading text{}" |
        grep 'NZ_' |
        sed 's/\..*//')

        location=$(cat product/antismash_${level}/${strain}/index.html |
        pup "div.page#${cluster} div.description-text text{}" |
        grep 'Location:' |
        sed -E "s/.*Location: ([0-9,]+) - ([0-9,]+) nt\. \(total: ([0-9,]+) nt\).*/\1\t\2\t\3/" )

        echo "${strain}\t${level}\t${cluster}\t${locus}\t${location}" \
        >> "$output_file"
        
    done
done

# extract BGC and flanking sequence with python script
python antismash_flanking.py aa1/bgc_flanking_10k/bgc_10k_locus_location.tsv product aa1/bgc_flanking_10k/sequence
cat aa1/bgc_flanking_10k/sequence/* >> aa1/bgc_flanking_10k/bgc_flanking_all.txt

```


```shell
# use BGC gbk file to clinker
cd ~/project/Actionmycetes/antismash/antismash_summary

for level in complete_taxon complete_untaxon; do
    product_file="product/cluster_gpa_${level}.tsv"

    # copy .gbk file to new dictionary    just region1
    cat "$product_file" | grep -v '^$' |
    while IFS=$'\t' read -r strain cluster product _ _ _; do
        clu=$(echo $cluster | awk 'match($0, /r[0-9]*c([0-9]+)/) {
            num = substr($0, RSTART + 3, RLENGTH - 3);
            printf "region%03d\n", num;
            }')
        echo "${strain}_${clu}"
        cluster_gbk_file=$(ls product/antismash_${level}/${strain}/*${clu}.gbk | head -n 1)
        echo "${cluster_gbk_file}"

        cp ${cluster_gbk_file} ~/project/Actionmycetes/cluster_clinker/gbk_file/${strain}_${cluster}.gbk
        
    done
done

# There are 'three' regions need to replace correct gbk file by hand

cd ~/project/Actionmycetes/antismash/cluster_clinker
clinker gbk_file/*gbk -s bgc_clinker_68.json -p bgc_cliner_68_raw.svg

```


```shell
# known GPA genome run antismash
cd ~/project/Actionmycetes/antismash/known_gpa_antismash
source ~/miniconda3/bin/activate
conda activate antismash

cat known_gpa_strain_list.txt |
parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
    echo {1};
    antismash --taxon bacteria -c 4 --cb-general --cc-mibig --cb-knownclusters  --pfam2go --asf --genefinding-tool prodigal genome/{1}/*genomic.gbff.gz --output-dir antismash_result/{1}
    "


```


```shell
# extract domain fron known GPA in mibig file

cd ~/project/Actionmycetes/antismash/antismash_summary

for i in $(cat product/knownGPA_mibig_list.txt | sed "s/\t/,/g"); do
    echo ${i};
    sample=$(echo ${i} | cut -d "," -f 1);
    num=$(echo ${i} | cut -d "," -f 2);
    gpa=$(echo ${i} | cut -d "," -f 3)
    js="product/knownGPA_mibig_antismash/${sample}/regions.js";
    type=$(python aa_antismash_pp_for_complete_uncomplete.py "${js}" "${num}" "${sample}" "${gpa}");
    echo ${type} | sed "s/]/]\n/g" |
    sed "s/\[//g" | sed "s/\]//g"| sed "s/'//g" | sed "s/,/\n/g" |
    sed "s/ >/>/g" | sed "s/+/\t/g" |
    sed "s/\t/\n/g" |
    sed '/^$/d' \
    >> aa1/domain_aa/knownGPA_mibig/domain_aa_knownGPA_mibig.txt;
done

# extract domain C\A\T
cat aa1/domain_aa/knownGPA_mibig/domain_aa_knownGPA_mibig.txt |
grep -A 1 -E "Condensation|Cglyc" |
sed '/^--$/d' \
>> aa1/domain_aa/knownGPA_mibig/domain_Caa.txt

cat aa1/domain_aa/knownGPA_mibig/domain_aa_knownGPA_mibig.txt |
grep -A 1 -E "AMP-binding" |
sed '/^--$/d' \
>> aa1/domain_aa/knownGPA_mibig/domain_Aaa.txt

cat aa1/domain_aa/knownGPA_mibig/domain_aa_knownGPA_mibig.txt |
grep -A 1 -E "\-PCP\-" |
sed '/^--$/d' \
>> aa1/domain_aa/knownGPA_mibig/domain_Taa.txt

# combine with previous C/A/T domain aminoacid sequence file in complete
for i in C A T; do
    cat aa1/domain_aa/domain_${i}aa_68.txt aa1/domain_aa/knownGPA_mibig/domain_${i}aa.txt >> aa1/domain_aa/domain_${i}aa_79.txt
    cp aa1/domain_aa/domain_${i}aa_79.txt ~/project/Actionmycetes/analysis_79/sequence_aa/domain_${i}aa_79.txt
done


cd ~/project/Actionmycetes/analysis_79
for i in C A T; do
    cat sequence_aa/domain_${i}aa_79.txt | grep '>' | sed 's/>//g' >> annotate/${i}domain_79_annotate.txt
done

for i in C A T; do
    mafft --auto sequence_aa/domain_${i}aa_79.txt > msa/domain_${i}aa_79.aln.fa
    trimal -in msa/domain_${i}aa_79.aln.fa -out trim/domain_${i}aa_79.trim.fa -automated1
    Fasttree trim/domain_${i}aa_79.trim.fa > fasttree/domain_${i}aa_79.nwk
done

```


```shell
# extract P450 locus_tag and translation of knwon bgc result from region.js
cd ~/project/Actionmycetes/antismash/antismash_summary

    product_file="product/knownGPA_mibig_list.txt"

    P450_file="aa1/p450_aa/p450_gene_aa_raw_11.txt"

    cat "$product_file" | grep -v '^$' |
    while IFS=$'\t' read -r strain cluster _ _ _ _; do
        echo "${strain}\t${cluster}"
        js_file=product/knownGPA_mibig_antismash/${strain}/regions.js
        python antismash_js_p450_11.py ${js_file} ${cluster} ${strain} \
        >> "$P450_file"
    done


cat aa1/p450_aa/p450_gene_aa_raw_11.txt | grep '>' | 
    sed 's/>//g' | sed -E "s/-/\t/g; s/\+/\t/g" \
    >> aa1/p450_aa/p450_identifier_raw_11.tsv

# filter the true P450 gene identifier by hand

cat aa1/p450_aa/p450_identifier_11.tsv | tsv-select -f 3 >> aa1/p450_aa/identifier_11.tsv

# then filter the raw fasta file
cat aa1/p450_aa/p450_gene_aa_raw_11.txt | grep -A 1 -F -i -f aa1/p450_aa/identifier_11.tsv | sed '/--/d' \
    >> aa1/p450_aa/p450_gene_aa_11.txt

# tree
cd ~/project/Actionmycetes/analysis_p450
mafft --auto sequence_aa/p450_gene_aa_79.txt > msa/p450_aa_79.aln.fa
trimal -in msa/p450_aa_79.aln.fa -out trim/p450_aa_79.trim.fa -automated1
Fasttree trim/p450_aa_79.trim.fa  > fasttree/p450_aa_79.nwk
```

```shell
# delete 4 abnormal BGC and fasttree again
cd ~/project/Actionmycetes/analysis_75


for i in C A T; do
    mafft --auto sequence_aa/domain_${i}aa_75.txt > msa/domain_${i}aa_75.aln.fa
    trimal -in msa/domain_${i}aa_75.aln.fa -out trim/domain_${i}aa_75.trim.fa -automated1
    Fasttree trim/domain_${i}aa_75.trim.fa > fasttree/domain_${i}aa_75.nwk
done

for i in C A T; do
    script=fasttree/table2itol/table2itol.R
    Rscript ${script} -a -D fasttree/domain_${i}aa_75  -i Strain -l Strain -w 0.5 annotate/${i}domain_75_annotate.xlsx
done

```

```shell
# filter the C domain and A domain of the central Hpg
cd ~/project/Actionmycetes/analysis_75

cat annotate/Adomain_75_annotate.xlsx | grep 'A_H0' | tsv-select -f 1 >> sequence_aa/Adomain_75_Hpg_filter.txt
cat annotate/Cdomain_75_annotate.xlsx | grep 'C_H0' | tsv-select -f 1 >> sequence_aa/Cdomain_75_Hpg_filter.txt
#there is only strain colum in files


# add the C domain and A domain of fegmycin's Hpg

cat sequence_aa/domain_Caa_75.txt | grep -A 1 -F -i -f sequence_aa/Cdomain_75_Hpg_filter.txt | sed '/--/d' >> sequence_aa/domain_Caa_75_Hpg.txt
cat sequence_aa/domain_Aaa_75.txt | grep -A 1 -F -i -f sequence_aa/Adomain_75_Hpg_filter.txt | sed '/--/d' >> sequence_aa/domain_Aaa_75_Hpg.txt

for i in C A; do
    mafft --auto sequence_aa/domain_${i}aa_75_Hpg.txt > msa/domain_${i}aa_75_Hpg.aln.fa
    trimal -in msa/domain_${i}aa_75_Hpg.aln.fa -out trim/domain_${i}aa_75_Hpg.trim.fa -automated1
    Fasttree trim/domain_${i}aa_75_Hpg.trim.fa > fasttree/domain_${i}aa_75_Hpg.nwk
done

for i in C A; do
    script=fasttree/table2itol/table2itol.R
    Rscript ${script} -a -D fasttree/domain_${i}aa_75_Hpg  -i Strain -l Strain -w 0.5 annotate/${i}domain_75_Hpg_annotate.xlsx
done

Rscript ${script} -a -D fasttree/domain_Aaa_75_Hpg_2  -i Strain -l Strain -w 0.5 annotate/Adomain_75_Hpg_annotate.xlsx
```