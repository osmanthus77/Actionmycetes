
# 使用antismash对基因组数据进行预测

## 1 将基因组数据进行分类

在`genomes_download.md`中下载完成基因组数据，并进行检验、minhash及计数，共得到9个科的共计12302菌株的基因组数据，另有一个菌株（Allon_opa_DSM_45601_GCF_003002095_1）科水平分类不确定   

| item    | count |
| ------- | ----: |
| strain  | 12302 |
| species |  1783 |
| genus   |   140 |
| family  |    10 |
| order   |     5 |
| class   |     1 |

9个科包括：
| family                |
| --------------------- |
| Actinomycetaceae      |
| Carbonactinosporaceae |
| Micromonosporaceae    |
| Nocardiopsidaceae     |
| Pseudonocardiaceae    |
| Streptomycetaceae     |
| Streptosporangiaceae  |
| Thermomonosporaceae   |
| Treboniaceae          |


### 1.1 根据组装水平以及是否是未知物种分类

根据`Actionmycetes.assembly.tsv`（原始基因组）及`collect.pass.tsv`（经过n50pass的、质量较高的），先分别生成4个tsv，complete_taxon、complete_untaxon、uncomplete_taxon、uncomplete_untaxon，后续运行antismash根据这个分类，其中以pass的进行antismash分析

> complete：Complete Genome、Chromosome
> 
> uncomplete：Contig、Scaffold
> 
> taxon：菌株有物种分类
> 
> untaxon：菌株无物种分类，名称上包含`sp.`

```shell
cd ~/project/nwr/Actionmycetes/summary

cat Actionmycetes.assembly.tsv |
    grep -E "Complete\ Genome|Chromosome" |
    grep -v "sp\." |
    tsv-select -f 1,4 \
    > Act_complete_taxon.tsv

cat Actionmycetes.assembly.tsv |
    grep -E "Complete\ Genome|Chromosome" |
    grep "sp\." |
    tsv-select -f 1,4 \
    > Act_complete_untaxon.tsv

cat Actionmycetes.assembly.tsv |
    grep -E -v "Complete\ Genome|Chromosome" |
    grep -v "sp\." |
    tsv-select -f 1,4 \
    > Act_uncomplete_taxon.tsv

cat Actionmycetes.assembly.tsv |
    grep -E -v "Complete\ Genome|Chromosome" |
    grep "sp\." |
    tsv-select -f 1,4 \
    > Act_uncomplete_untaxon.tsv

# 以下四个文件为后续需要使用的
cat collect.pass.tsv |
    grep -E "Complete\ Genome|Chromosome" |
    grep -v "sp\." |
    tsv-select -f 1,2 \
    > Act_pass_complete_taxon.tsv

cat collect.pass.tsv |
    grep -E "Complete\ Genome|Chromosome" |
    grep "sp\." |
    tsv-select -f 1,2 \
    > Act_pass_complete_untaxon.tsv

cat collect.pass.tsv |
    grep -E -v "Complete\ Genome|Chromosome" |
    grep -v "sp\." |
    tsv-select -f 1,2 \
    > Act_pass_uncomplete_taxon.tsv

cat collect.pass.tsv |
    grep -E -v "Complete\ Genome|Chromosome" |
    grep "sp\." |
    tsv-select -f 1,2 \
    > Act_pass_uncomplete_untaxon.tsv

wc -l Act_pass*
#     850 Act_complete_taxon.tsv_pass.tsv
#     709 Act_complete_untaxon.tsv_pass.tsv
#    6060 Act_uncomplete_taxon.tsv_pass.tsv
#    4912 Act_uncomplete_untaxon.tsv_pass.tsv
#   12531 total

# 10个文件+ genome.taxon.tsv 全部复制到 ~/project/Actionmycetes/antismash/group
```

### 1.2 根据科分类进行再次分类


```shell
cd ~/project/Actionmycetes/antismash

for family in $(cat ../genomes/Count/family.lst); do
    mkdir -p group/family_GCF
    mkdir -p group/family_result
    cat group/genome.taxon.tsv | grep ${family} | cut -f 1 > group/family_GCF/${family}_GCF.tsv
    cat group/Act_complete_taxon_pass.tsv | tsv-join -f group/family_GCF/${family}_GCF.tsv -k 1 > group/family_result/${family}_pass_complete_taxon.tsv
    cat group/Act_complete_untaxon_pass.tsv | tsv-join -f group/family_GCF/${family}_GCF.tsv -k 1 > group/family_result/${family}_pass_complete_untaxon.tsv
    cat group/Act_uncomplete_taxon_pass.tsv | tsv-join -f group/family_GCF/${family}_GCF.tsv -k 1 > group/family_result/${family}_pass_uncomplete_taxon.tsv
    cat group/Act_uncomplete_untaxon_pass.tsv | tsv-join -f group/family_GCF/${family}_GCF.tsv -k 1 > group/family_result/${family}_pass_uncomplete_untaxon.tsv
done

wc -l group/genome.taxon.tsv
# 12302 group/genome.taxon.tsv

wc -l group/family_result/*
```

```txt
wc -l group/family_result/*
     105 group/family_result/Actinomycetaceae_pass_complete_taxon.tsv
      18 group/family_result/Actinomycetaceae_pass_complete_untaxon.tsv
     393 group/family_result/Actinomycetaceae_pass_uncomplete_taxon.tsv
      99 group/family_result/Actinomycetaceae_pass_uncomplete_untaxon.tsv
       0 group/family_result/Carbonactinosporaceae_pass_complete_taxon.tsv
       0 group/family_result/Carbonactinosporaceae_pass_complete_untaxon.tsv
       3 group/family_result/Carbonactinosporaceae_pass_uncomplete_taxon.tsv
       0 group/family_result/Carbonactinosporaceae_pass_uncomplete_untaxon.tsv
      66 group/family_result/Micromonosporaceae_pass_complete_taxon.tsv
      65 group/family_result/Micromonosporaceae_pass_complete_untaxon.tsv
     558 group/family_result/Micromonosporaceae_pass_uncomplete_taxon.tsv
     331 group/family_result/Micromonosporaceae_pass_uncomplete_untaxon.tsv
      18 group/family_result/Nocardiopsidaceae_pass_complete_taxon.tsv
       4 group/family_result/Nocardiopsidaceae_pass_complete_untaxon.tsv
      95 group/family_result/Nocardiopsidaceae_pass_uncomplete_taxon.tsv
      30 group/family_result/Nocardiopsidaceae_pass_uncomplete_untaxon.tsv
      61 group/family_result/Pseudonocardiaceae_pass_complete_taxon.tsv
      49 group/family_result/Pseudonocardiaceae_pass_complete_untaxon.tsv
     417 group/family_result/Pseudonocardiaceae_pass_uncomplete_taxon.tsv
     173 group/family_result/Pseudonocardiaceae_pass_uncomplete_untaxon.tsv
     537 group/family_result/Streptomycetaceae_pass_complete_taxon.tsv
     549 group/family_result/Streptomycetaceae_pass_complete_untaxon.tsv
    4029 group/family_result/Streptomycetaceae_pass_uncomplete_taxon.tsv
    4043 group/family_result/Streptomycetaceae_pass_uncomplete_untaxon.tsv
      13 group/family_result/Streptosporangiaceae_pass_complete_taxon.tsv
       9 group/family_result/Streptosporangiaceae_pass_complete_untaxon.tsv
     289 group/family_result/Streptosporangiaceae_pass_uncomplete_taxon.tsv
     173 group/family_result/Streptosporangiaceae_pass_uncomplete_untaxon.tsv
       8 group/family_result/Thermomonosporaceae_pass_complete_taxon.tsv
       6 group/family_result/Thermomonosporaceae_pass_complete_untaxon.tsv
     114 group/family_result/Thermomonosporaceae_pass_uncomplete_taxon.tsv
      45 group/family_result/Thermomonosporaceae_pass_uncomplete_untaxon.tsv
       0 group/family_result/Treboniaceae_pass_complete_taxon.tsv
       0 group/family_result/Treboniaceae_pass_complete_untaxon.tsv
       1 group/family_result/Treboniaceae_pass_uncomplete_taxon.tsv
       0 group/family_result/Treboniaceae_pass_uncomplete_untaxon.tsv
   12301 total
```

> 组装水平及物种分类是以`collect.pass.tsv`为基础的，该文件只经过n50检验，未完成minhash和重新count，因此该文件中还包含`minhash`筛选掉的异常菌株。这里结果为12531个
> 
> 科分类是以`genome.taxon.tsv`为基础的，该文件经过n50检验、minhash、重新count之后生成。这里结果为12301

### 1.3 补全缺失信息

```shell
cd ~/project/Actionmycetes/antismash
mkdir -p group/omission

# 12301
for i in $(ls group/family_result/*); do
    cat ${i} |
    cut -f 1 \
    >> group/omission/all_12301_name.tsv;
done
# 12531
for i in $(ls group/Act_pass*.tsv); do
    cat ${i} |
    cut -f 1 \
    >> group/omission/all_12531_name.tsv;
done

# 排序
sort group/omission/all_12301_name.tsv > group/omission/all_12301_name_sort.tsv
sort group/omission/all_12531_name.tsv > group/omission/all_12531_name_sort.tsv

# 筛选缺少的 只在file1出现，即 12531 中
comm -23 group/omission/all_12531_name_sort.tsv group/omission/all_12301_name_sort.tsv > group/omission/omission_name.tsv

head -n 5 group/omission/omission_name.tsv
# Actinoma_sp_DLS_62_GCF_041463655_1
# Actinomyces_nae_GCF_901873715_1
# Actinomyces_oris_CCUG_33920_GCF_001937675_1
# Actinomyces_oris_ORNL_0101_GCF_034648895_1
# Actinomyces_oris_P6N_GCF_001937665_1

# 补充完整信息
cat group/Actionmycetes.assembly.tsv | tsv-join -f group/omission/omission_name.tsv -k 1 > group/omission/omission_collect_pass.tsv
```

```shell
# 查看omission中组装水平
cat group/omission/omission_collect_pass.tsv | tsv-select -f 5 | tsv-uniq
# Scaffold
# Contig
# Complete Genome
# Chromosome

# 检查物种分类情况
cat group/omission/omission_collect_pass.tsv | tsv-select -f 4 | cut -d " " -f 2 | tsv-uniq
# 输出包括 sp. 和 物种名
```
根据这些输出，按照组装水平和物种分类，跟之前的一致，omission中仍然按照这个方法分类

下面，检查分类学上属于什么genus

```shell
# 检查分类学上属于什么genus
cat group/omission/omission_collect_pass.tsv | tsv-select -f 4 | cut -d " " -f 1 | tsv-uniq

Actinomadura
Actinomyces
Allonocardiopsis 
Amycolatopsis
Buchananella
Gleimia
Kitasatospora
Longispora
Microbispora
Micromonospora
Nonomuraea
Pseudonocardia
Saccharopolyspora
Schaalia
Streptacidiphilus
Streptomyces
Trueperella
Varibaculum
Winkia
```

再确定genus属于什么family

> Thermomonosporaceae: Actinomadura
> 
> Actinomycetaceae: Actinomyces、Buchananella、Gleimia、Schaalia、Trueperella、Varibaculum、Winkia
> 
> Pseudonocardiaceae: Amycolatopsis、Pseudonocardia、Saccharopolyspora
> 
> Streptomycetaceae: Kitasatospora、Streptacidiphilus、Streptomyces
> 
> Streptosporangiaceae: Microbispora、Nonomuraea
> 
> Micromonosporaceae: Longispora、Micromonospora
> 
> no rank: Allonocardiopsis, 搜索`collect.pass.tsv`后对应菌株只有一个（Allon_opa_DSM_45601_GCF_003002095_1）

再将这些分到各自的family tsv中

```bash
declare -A family_genera_map
family_genera_map["Thermomonosporaceae"]="Actinomadura"
family_genera_map["Actinomycetaceae"]="Actinomyces|Buchananella|Gleimia|Schaalia|Trueperella|Varibaculum|Winkia"
family_genera_map["Pseudonocardiaceae"]="Amycolatopsis|Pseudonocardia|Saccharopolyspora"
family_genera_map["Streptomycetaceae"]="Kitasatospora|Streptacidiphilus|Streptomyces"
family_genera_map["Streptosporangiaceae"]="Microbispora|Nonomuraea"
family_genera_map["Micromonosporaceae"]="Longispora|Micromonospora"

input_file="group/omission/omission_collect_pass.tsv"
output_dir="group/omission"

for family in "${!family_genera_map[@]}"; do
    echo "处理: ${family}"
    genera_pattern="${family_genera_map[$family]}"
    temp_file=$(mktemp)

    grep -E "${genera_pattern}" "${input_file}" > "${temp_file}"

    grep -E "Complete Genome|Chromosome" "${temp_file}" | grep -v "sp\." > "${output_dir}/o_${family}_complete_taxon.tsv"
    grep -E "Complete Genome|Chromosome" "${temp_file}" | grep "sp\." > "${output_dir}/o_${family}_complete_untaxon.tsv"
    grep -E "Scaffold|Contig" "${temp_file}" | grep -v "sp\." > "${output_dir}/o_${family}_uncomplete_taxon.tsv"
    grep -E "Scaffold|Contig" "${temp_file}" | grep "sp\." > "${output_dir}/o_${family}_uncomplete_untaxon.tsv"

    rm "${temp_file}" # 删除临时文件
done
```

然后再补充回去
```shell
for family in Thermomonosporaceae Actinomycetaceae Pseudonocardiaceae Streptomycetaceae Streptosporangiaceae Micromonosporaceae; do
    for level in complete uncomplete; do
        for species in taxon untaxon; do
            cat group/omission/o_${family}_${level}_${species}.tsv | cut -f 1,4 \
            >> group/family_result/${family}_pass_${level}_${species}.tsv
        done
    done
done
```

## 2 antismash（complete水平）

### 2.1 检查gbff文件

```shell
cd ~/project/Actionmycetes/antismash
mkdir -p run/check
ln -s ~/project/Actionmycetes/genomes/ASSEMBLY ASSEMBLY

for family in $(cat ../genomes/Count/family.lst); do
    cat group/family_result/${family}_*.tsv | tr " " "_" |
    sed "s/sp\./sp/g" |
    sed "s/(//g" |
    sed "s/)//g" |
    sed "s/\.//g" |
    sed "s/\[//g" |
    sed "s/\]//g" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 '
    if [ ! -e ASSEMBLY/{2}/{1}/*genomic.gbff.gz ]; then
        echo ${family}  ASSEMBLY/{2}/{1} >> run/check/error.tsv;
    else 
        echo ${family} ASSEMBLY/{2}/{1}  >> run/check/success.tsv;
    fi
    '
done

    cat group/family_result/${family}_*.tsv | tr " " "_" |
    sed "s/sp\./sp/g" |
    sed "s/(//g" |
    sed "s/)//g" |
    sed "s/\.//g" |
    sed "s/\[//g" |
    sed "s/\]//g" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 '
    if [ ! -e ASSEMBLY/{2}/{1}/*genomic.gbff.gz ]; then
        echo ${family}  ASSEMBLY/{2}/{1} >> run/check/error.tsv;
    else 
        echo ${family} ASSEMBLY/{2}/{1}  >> run/check/success.tsv;
    fi
    '
done


wc -l run/check/success.tsv
#   12530 run/check/success.tsv
# 缺少的一个为没有family分类的菌株
```

### 2.2 运行antismash

```shell
# 选取几个作为测试
for family in Actinomycetaceae; do
    for level in complete_taxon complete_untaxon; do
        mkdir -p test/${family}/${level};
        cat test/${family}_pass_${level}.tsv | tr " " "_" |
        sed "s/sp\./sp/g" |
        sed "s/(//g" |
        sed "s/)//g" |
        sed "s/\.//g" |
        sed "s/\[//g" |
        sed "s/\]//g" |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 2 "
            echo "${family}_${level}_{2}_{1}";
            mkdir -p test/${family}/${level}/{2};
            antismash --taxon bacteria -c 4 --cb-general --cc-mibig --cb-knownclusters --pfam2go ASSEMBLY/{2}/{1}/*genomic.gbff.gz --output-dir test/${family}/${level}/{2}/{1}
        "
    done
done


# 正式运行
for family in $(cat ../genomes/Count/family.lst); do
    for level in complete_taxon complete_untaxon; do
        mkdir -p antismash_result/${family}/${level};
        cat group/family_result/${family}_pass_${level}.tsv | tr " " "_" |
        sed "s/sp\./sp/g" |
        sed "s/(//g" |
        sed "s/)//g" |
        sed "s/\.//g" |
        sed "s/\[//g" |
        sed "s/\]//g" |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            echo "${family}_${level}_{2}_{1}";
            mkdir -p antismash_result/${family}/${level}/{2};
            antismash --taxon bacteria -c 4 --cb-general --cc-mibig --cb-knownclusters --pfam2go ASSEMBLY/{2}/{1}/*genomic.gbff.gz --output-dir antismash_result/${family}/${level}/{2}/{1}
        "
    done
done
```
