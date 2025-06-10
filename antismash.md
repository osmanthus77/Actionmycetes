
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
# 激活环境
cd ~/project/Actionmycetes/antismash
source ~/miniconda3/bin/activate
conda antivate antismash

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
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            echo "${family}_${level}_{2}_{1}";
            time antismash --taxon bacteria -c 6 --cb-general --cc-mibig --cb-knownclusters --pfam2go ASSEMBLY/{2}/{1}/*genomic.gbff.gz --output-dir test/${family}/${level}/{1}
        "
    done
done


# 正式运行
for family in $(cat family.lst); do
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
            antismash --taxon bacteria -c 4 --cb-general --cc-mibig --cb-knownclusters --pfam2go ASSEMBLY/{2}/{1}/*genomic.gbff.gz --output-dir antismash_result/${family}/${level}/{1}
        "
    done
done
```

### 2.3 检查并补齐

```shell
cd ~/project/Actionmycetes/antismash/run
mkdir completion

# 检查生成结果目录的数量
for family in $(cat ../family.lst); do
    for level in complete_taxon complete_untaxon; do
        ls ../antismash_result/${family}/${level} |
        sed "s/$/\t${family}\t${level}/g" >> completion/all_end_names.tsv
    done
done

wc -l completion/all_end_names.tsv 
#    1556 completion/all_end_names.tsv

# 检查空目录
for family in $(cat ../family.lst); do
    for level in complete_taxon complete_untaxon; do
        find ../antismash_result/${family}/${level} -mindepth 1 -maxdepth 1 -type d -empty |
        cut -f 1 |
        cut -d "/" -f 5 |
        sed "s/$/\t${family}\t${level}/g" >> completion/all_empty_names.tsv
    done
done

wc -l completion/all_empty_names.tsv 
#      24 completion/all_empty_names.tsv
```

`1.1`中统计得到complete共1559个条目，antismash得到1556个目录，数量不匹配，并且存在24个空目录。根据运行过程中stdout，可知是部分菌株gbff文件内容出现问题，无法成功注释。可使用fna文件注释

```shell
# 寻找1559条中无注释部分

# 1559条菌株信息
cd ~/project/Actionmycetes/antismash/run

for family in $(cat ../family.lst); do
    for level in complete_taxon complete_untaxon; do
        cat ../group/family_result/${family}_pass_${level}.tsv |
        cut -f 1 |
        sed "s/$/\t${family}\t${level}/g" >> completion/all_complete_names.tsv
    done
done
wc -l all_complete_names.tsv 
#    1559 all_complete_names.tsv

# 找出没跑的
cd ~/project/Actionmycetes/antismash/run/completion
sort all_end_names.tsv > all_end_names_sort.tsv
sort all_complete_names.tsv > all_complete_names_sort.tsv

comm -23 all_complete_names_sort.tsv all_end_names_sort.tsv > no_run.tsv
wc -l no_run.tsv 
#       4 no_run.tsv

# 合并空目录和没跑的
cat all_empty_names.tsv | cut -f 1 >> run_wrong.tsv
cat no_run.tsv | cut -f 1 >> run_wrong.tsv

wc -l run_wrong.tsv 
#      28 run_wrong.tsv

# 获取genus-species信息
for family in $(cat ../../family.lst); do
    for level in complete_taxon complete_untaxon; do
        cat  ../../group/family_result/${family}_pass_${level}.tsv |
        tr "" "_" |
        sed "s/sp\./sp/g" |
        sed "s/(//g" |
        sed "s/)//g" |
        sed "s/\.//g" |
        sed "s/\[//g" |
        sed "s/\]//g" >> all_complete_names_taxon.tsv
    done 
done

# 找出需要重跑 antismash 的，并整合信息
sort run_wrong.tsv >  run_wrong_sort.tsv
sort all_complete_names_taxon.tsv > all_complete_names_taxon_sort.tsv

cat all_complete_names_taxon_sort.tsv | tsv-join -f run_wrong_sort.tsv -k 1 > run_wrong_taxon_sort.tsv
wc -l run_wrong_taxon_sort.tsv 
#      28 run_wrong_taxon_sort.tsv

# 并补充 antismash 运行过程中出现 error 和 warning 的条目
wc -l run_wrong_taxon_sort.tsv 
#      36 run_wrong_taxon_sort.tsv
```

一共统计得到28条菌株需要重新运行 antismash

```shell
# 利用fna基因组文件进行antismash
cd ~/project/Actionmycetes/antismash/run/completion
source ~/miniconda3/bin/activate
conda antivate antismash
mkdir fna result

cat run_wrong_taxon_sort.tsv | tr " " "_" |
    parallel --colsep "\t" --no-run-if-empty --linebuffer -k -j 1 '
    echo "{2}\t{1}";
    mkdir fna/{1};
    cp ../../ASSEMBLY/{2}/{1}/*genomic.fna.gz fna/{1}/ ;
    rm fna/{1}/*cds_from_genomic.fna.gz;
    rm fna/{1}/*rna_from_genomic.fna.gz;
    gzip -df fna/{1}/*.gz;
    antismash --taxon bacteria -c 4 --cb-general --cc-mibig --cb-knownclusters \
        --pfam2go --asf --genefinding-tool prodigal \
        fna/{1}/*genomic.fna --output-dir result/{1}
    '
# 手动移动到对应目录下
```

### 2.4 检查补齐结果

```shell
cd ~/project/Actionmycetes/antismash
mkdir -p antismash_summary/strains_raw

# 统计 strain 信息
for family in $(cat family.lst); do
    for level in complete_taxon complete_untaxon; do
        find antismash_result/${family}/${level}  -mindepth 1 -maxdepth 1 -type d |
            xargs -I {} basename {} > antismash_summary/strains_raw/strains_${family}_${level}_8.lst ;
        wc -l antismash_summary/strains_raw/strains_${family}_${level}_8.lst
    done
done

wc -l antismash_summary/strains_raw/*       
#      22 antismash_summary/strains_raw/strains_Actinomycetaceae_complete_untaxon_8.lst
#       0 antismash_summary/strains_raw/strains_Carbonactinosporaceae_complete_taxon_8.lst
#       0 antismash_summary/strains_raw/strains_Carbonactinosporaceae_complete_untaxon_8.lst
#      66 antismash_summary/strains_raw/strains_Micromonosporaceae_complete_taxon_8.lst
#      65 antismash_summary/strains_raw/strains_Micromonosporaceae_complete_untaxon_8.lst
#      18 antismash_summary/strains_raw/strains_Nocardiopsidaceae_complete_taxon_8.lst
#       4 antismash_summary/strains_raw/strains_Nocardiopsidaceae_complete_untaxon_8.lst
#      61 antismash_summary/strains_raw/strains_Pseudonocardiaceae_complete_taxon_8.lst
#      49 antismash_summary/strains_raw/strains_Pseudonocardiaceae_complete_untaxon_8.lst
#     576 antismash_summary/strains_raw/strains_Streptomycetaceae_complete_taxon_8.lst
#     553 antismash_summary/strains_raw/strains_Streptomycetaceae_complete_untaxon_8.lst
#      13 antismash_summary/strains_raw/strains_Streptosporangiaceae_complete_taxon_8.lst
#      10 antismash_summary/strains_raw/strains_Streptosporangiaceae_complete_untaxon_8.lst
#       8 antismash_summary/strains_raw/strains_Thermomonosporaceae_complete_taxon_8.lst
#       6 antismash_summary/strains_raw/strains_Thermomonosporaceae_complete_untaxon_8.lst
#       0 antismash_summary/strains_raw/strains_Treboniaceae_complete_taxon_8.lst
#       0 antismash_summary/strains_raw/strains_Treboniaceae_complete_untaxon_8.lst
#    1559 total

# 检查空目录
for family in $(cat family.lst); do
    for level in complete_taxon complete_untaxon; do
        find antismash_result/${family}/${level}  -mindepth 1 -maxdepth 1 -type d -empty |
        cut -f 1 |
        cut -d "/" -f 4 |
        sed "s/$/\t${family}\t${level}/g" >> antismash_summary/empty_dir.tsv
    done
done

wc -l antismash_summary/empty_dir.tsv
#       0 antismash_summary/empty_dir.tsv
```

空目录数目为0，说明结果补齐成功，complete的全部菌株antismash都运行成功

## 3 统计结果信息（complete水平）

### 3.1 统计 overview-table 信息

```shell
cd  ~/project/Actionmycetes/antismash/antismash_summary
mkdir -p table/overview table/mibig

for family in $(cat ../family.lst); do
    for level in complete_taxon complete_untaxon; do
        for strains in $(cat strains_raw/strains_${family}_${level}_8.lst); do
            mkdir -p table/overview/raw;
            cat ../antismash_result/${family}/${level}/${strains}/index.html |
            pup 'table.region-table tbody tr td text{}' |
            sed 's/Region/|Region/g' |
            grep '\S' > table/overview/raw/${strains}_overview_raw.tsv
        done
    done
done

# 修改格式
for family in $(cat ../family.lst); do
    for level in complete_taxon complete_untaxon; do
        for strains in $(cat strains_raw/strains_${family}_${level}_8.lst); do
            perl html.pl table/overview/raw/${strains}_overview_raw.tsv |
            sed "s/^/${strains}_/g; s/Region /cluster/g" |
            sed -E 's/cluster1\./cluster/g'
            sort | uniq \
            >> table/overview/overview_whole_${level}_8.tsv
        done
    done
done

wc -l table/overview/*tsv
#   23870 table/overview/overview_whole_complete_taxon_8.tsv
#   21812 table/overview/overview_whole_complete_untaxon_8.tsv
#   45682 total


```
### 3.2 统计 mibig 产物信息

```shell
cd  ~/project/Actionmycetes/antismash/antismash_summary

for family in $(cat ../family.lst); do
    for level in complete_taxon complete_untaxon; do
        for strains in $(cat strains_raw/strains_${family}_${level}_8.lst); do
            cat ../antismash_result/${family}/${level}/${strains}/index.html |
            pup 'table.region-table tbody tr.linked-row a attr{href}' |
            grep -v "https://docs.antismash.secondarymetabolites.org/" |    #排除 antismash 自身代谢产物数据库
            grep -B 1 "https://mibig.secondarymetabolites.org/" |    #筛选出根据 mibig 预测的信息
            grep -v "\-\-" |
            sed -E 's/https:\/\/mibig.secondarymetabolites.org\/go\///g; s/#//g; s/\/[0-9]//g' |
            paste - - |
            sed "s/^/${strains}_/g" |
            sed -E 's/r1c/cluster/g; s/r([0-9]+)c([0-9]+)/cluster\1.\2/g' |
            sort | uniq \
            >> table/mibig/mibig_whole_${level}_8.tsv
        done 
    done
done

# 统计预测产物MiBIG参考信息的个数
wc -l table/mibig/*tsv
#   13269 table/mibig/mibig_whole_complete_taxon_8.tsv
#   11593 table/mibig/mibig_whole_complete_untaxon_8.tsv
#   24862 total
```


### 3.3 筛选产物

创建GPA名称的list文件，根据该文件对antismash结果进行筛选

```shell
cd  ~/project/Actionmycetes/antismash/antismash_summary
mkdir product
for level in complete_taxon complete_untaxon; do
    cat table/overview/overview_whole_${level}_8.tsv |
        tsv-join -d 1 -f table/mibig/mibig_whole_${level}_8.tsv -k 1 --append-fields 2 | # 合并product和MiBIG
        tsv-filter --or --str-in-fld 7:High --str-in-fld 7:Medium  | # 筛选：相似度Medium和High,可选
        sed 's/_cluster/\tcluster/g' |
        grep -F -f "gpa.tsv" > product/gpa_${level}_8.tsv 
        wc -l product/gpa_${level}_8.tsv
done

wc -l product/*
#      43 product/gpa_complete_taxon_8.tsv
#      23 product/gpa_complete_untaxon_8.tsv

```

## 4 提取domain序列

### 4.1 提取全部domain序列

#### 准备产物文件

```shell
# 修改第二列格式
cd ~/project/Actionmycetes/antismash/antismash_summary

for level in complete_taxon complete_untaxon; do
    cat product/gpa_${level}_8.tsv | 
    sed 's/cluster/r1c/g' |
    cut -f 1,2 > product/cluster_gpa_${level}.tsv
done

# 注意：可能存在cluster2.1的情况，手动修改为 r2c1

# 拷贝相关antismash结果到新文件夹
for family in $(cat ../family.lst); do
    for level in complete_taxon complete_untaxon; do
        cat product/cluster_gpa_${level}.tsv |
            cut -f 1 |
            parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 "
                echo {};
                mkdir -p product/antismash_${level}/{};
                cp -r ../antismash_result/${family}/${level}/{} product/antismash_${level}/
                "
    done
done

```

#### 提取全部domain序列

手动将`antismash_pp_for_complete_uncomplete.py`脚本中`dna`替换为`sequence`并保存为另一文件，针对dna和aa的`py`脚本内容除`dna和aa`外无区别

```shell
cd ~/project/Actionmycetes/antismash/antismash_summary

for level in complete_taxon complete_untaxon; do
    for se in dna aa; do
        mkdir domain_${se};

        for i in $(cat product/cluster_gpa_${level}.tsv | sed "s/\t/,/g"); do
            echo ${i};
            sample=$(echo ${i} | cut -d "," -f 1);
            num=$(echo ${i} | cut -d "," -f 2);
            echo ${num};
            js="product/antismash_${level}/${sample}/regions.js";
            type=$(python ${se}_antismash_pp_for_complete_uncomplete.py "${js}" "${num}" "${sample}");
            echo ${type} | sed "s/]/]\n/g" >> domain_${se}/domain_${se}_all_${level}.txt;
        done
    done
done

# 修改格式
for level in complete_taxon complete_untaxon; do
    for se in dna aa; do
        cat domain_${se}/domain_${se}_all_${level}.txt |
        sed "s/\[//g" | sed "s/\]//g"|sed "s/'//g" | sed "s/,/\n/g"|
        sed "s/ >/>/g" | sed "s/+/\t/g" |
        sed "s/\t/\n/g" |
        sed '/^$/d' \
        > domain_${se}/ok_domain_${se}_all_${level}.txt
    done
done

```

### 4.2 提取特定domain序列

**注意**：筛选特定 domain 时，要用 domain 全称，如下表
| domain 简称 | domain 全称 |
|:----------:|--------------|
| A domain | AMP-binding |
| C domain | Condensation |
| T domain | PCP |

ps: 糖肽类抗生素中 `C domain`还有一种类型，叫`Cglyc`，筛选时需要与`Condensation`一起筛选，才是全部Cdomain

```shell
cd ~/project/Actionmycetes/antismash/antismash_summary

for level in complete_taxon complete_untaxon; do
    for se in dna aa; do
        cat domain_${se}/ok_domain_${se}_all_${level}.txt |
        grep -E "Condensation|Cglyc" \
        > domain_${se}/donmain_C${se}.txt

        cat domain_${se}/ok_domain_${se}_all_${level}.txt |
        grep -E "AMP-binding" \
        > domain_${se}/donmain_A${se}.txt

        cat domain_${se}/ok_domain_${se}_all_${level}.txt |
        grep -E "PCP" \
        > domain_${se}/donmain_T${se}.txt
    done
done


```