# 下载 *Actionmycetes* 基因组数据
分歧杆菌目( *Mycobacteriales* )通常是病原菌，下载基因组数据时排除掉
[Glycopeptide Antibiotics](https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1521-3773(19990802)38:15%3C2096::AID-ANIE2096%3E3.0.CO;2-F)列举了1999年之前发现的糖肽抗生素及其来源细菌，分别在Streptosporangiales、Pseudonocardiales、Kitasatosporales、Micromonosporales、Actinomycetales五个目中

## 1 统计分类信息

```shell
nwr member Actinomycetes -r order |
    keep-header -- tsv-sort -k2,2 |
    rgr md stdin --num

nwr member \
    Streptosporangiales Pseudonocardiales Kitasatosporales Micromonosporales Actinomycetales |
    tsv-summarize -H -g 3 --count |
    rgr md stdin --fmt

nwr member \
    Streptosporangiales Pseudonocardiales Kitasatosporales Micromonosporales Actinomycetales \
    -r "species group" -r "species subgroup" |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3 -k2,2 |
    rgr md stdin --num

```

| #tax_id | sci_name                   | rank  | division |
| ------: | -------------------------- | ----- | -------- |
| 1643683 | Acidothermales             | order | Bacteria |
|    2037 | Actinomycetales            | order | Bacteria |
|  622450 | Actinopolysporales         | order | Bacteria |
|   85004 | Bifidobacteriales          | order | Bacteria |
| 1389450 | Candidatus Actinomarinales | order | Bacteria |
| 2039638 | Candidatus Nanopelagicales | order | Bacteria |
|  414714 | Catenulisporales           | order | Bacteria |
| 2495577 | Cryptosporangiales         | order | Bacteria |
|   85013 | Frankiales                 | order | Bacteria |
| 1643682 | Geodermatophilales         | order | Bacteria |
|   85014 | Glycomycetales             | order | Bacteria |
| 2805415 | Jatrophihabitantales       | order | Bacteria |
| 1217098 | Jiangellales               | order | Bacteria |
|  622452 | Kineosporiales             | order | Bacteria |
|   85011 | Kitasatosporales           | order | Bacteria |
|   85006 | Micrococcales              | order | Bacteria |
|   85008 | Micromonosporales          | order | Bacteria |
| 2793120 | Motilibacterales           | order | Bacteria |
|   85007 | Mycobacteriales            | order | Bacteria |
| 1643684 | Nakamurellales             | order | Bacteria |
|   85009 | Propionibacteriales        | order | Bacteria |
|   85010 | Pseudonocardiales          | order | Bacteria |
| 2495578 | Sporichthyales             | order | Bacteria |
|   85012 | Streptosporangiales        | order | Bacteria |

| rank             |  count |
| ---------------- | -----: |
| order            |      5 |
| family           |      9 |
| no rank          |    173 |
| species          | 42,209 |
| genus            |    158 |
| strain           |    407 |
| subspecies       |    165 |
| clade            |      4 |
| species group    |     19 |
| species subgroup |     10 |

| #tax_id | sci_name                              | rank             |
| ------: | ------------------------------------- | ---------------- |
| 2893673 | Amycolatopsis japonica group          | species group    |
| 2893674 | Amycolatopsis methanolica group       | species group    |
| 2893671 | Prauserella coralliicola group        | species group    |
| 2893672 | Prauserella salsuginis group          | species group    |
| 1477431 | Streptomyces albidoflavus group       | species group    |
| 2867120 | Streptomyces albogriseolus group      | species group    |
| 1852274 | Streptomyces albus group              | species group    |
| 2867194 | Streptomyces althioticus group        | species group    |
| 2838335 | Streptomyces aurantiacus group        | species group    |
| 1481185 | Streptomyces cinnamoneus group        | species group    |
| 2849069 | Streptomyces diastaticus group        | species group    |
| 2867193 | Streptomyces griseoincarnatus group   | species group    |
|  629295 | Streptomyces griseus group            | species group    |
|  661169 | Streptomyces melanosporofaciens group | species group    |
| 2838332 | Streptomyces phaeochromogenes group   | species group    |
| 2894086 | Streptomyces pseudogriseolus group    | species group    |
| 2867164 | Streptomyces rochei group             | species group    |
| 2867121 | Streptomyces violaceoruber group      | species group    |
| 2839105 | Streptomyces violaceusniger group     | species group    |
| 1482558 | Streptomyces albovinaceus subgroup    | species subgroup |
| 1482561 | Streptomyces anulatus subgroup        | species subgroup |
| 1482564 | Streptomyces atroolivaceus subgroup   | species subgroup |
| 1482566 | Streptomyces bacillaris subgroup      | species subgroup |
| 1482592 | Streptomyces fimicarius subgroup      | species subgroup |
| 1482593 | Streptomyces flavovirens subgroup     | species subgroup |
| 1482596 | Streptomyces griseus subgroup         | species subgroup |
| 1482599 | Streptomyces halstedii subgroup       | species subgroup |
| 1482601 | Streptomyces microflavus subgroup     | species subgroup |
| 1482603 | Streptomyces puniceus subgroup        | species subgroup |

## 2 统计组装信息

组装信息分为两种：`RefSeq`（参考基因组）、`Genbank`（所有测序数据）   
这一步需要从 refseq 和 genbank 中生成得到包含下载信息的两个文件

```shell
mkdir -p ~/project/nwr/Actionmycetes/summary
cd ~/project/nwr/Actionmycetes/summary

nwr member \
    Streptosporangiales Pseudonocardiales Kitasatosporales Micromonosporales Actinomycetales \
    -r genus |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv 
# 158 genus.list.tsv

#refseq信息
cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -separator $'\t' ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS1.tsv

#genbank信息
cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -separator $'\t' ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB1.tsv

wc -l RS*.tsv GB*.tsv
# 8075 RS1.tsv
# 8975 GB1.tsv

#统计总数
for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1	13792
#GB1	16373
```

## 3 下载数据

### 3.1 生成下载的数据表格

- 统计模式物种信息

从genus.list.tsv中得到属名的标识号，再用SQL筛选模式物种信息存入reference.tsv

```shell
cd ~/project/nwr/Actionmycetes/summary

GENUS=$(
    cat genus.list.tsv|
        cut -f 1 |
        tr "\n" "," | 
        sed 's/,$//'
)

echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -separator $'\t' ~/.nwr/ar_refseq.sqlite \
    > reference.tsv

wc -l reference.tsv 
# 1732 reference.tsv

cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv
```

- refseq 中筛选上一步骤中得到的species信息

```shell
cd ~/project/nwr/Actionmycetes/summary

# refseq
SPECIES=$(
    cat RS1.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

# 排除未知种sp.和杂交种x
echo "
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% x %'  
    " |
    sqlite3 -separator $'\t' ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv
wc -l raw.tsv 
#7515 raw.tsv

# 添加未知种
echo "
    SELECT
        genus || ' sp. ' || infraspecific_name || ' ' || assembly_accession AS name,
        genus || ' sp.', genus, ftp_path, biosample, assembly_level,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species LIKE '% sp.%'
    " |
    sqlite3 -separator $'\t' ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv
wc -l raw.tsv
#13793 raw.tsv

# 单独保存 refseq 中 assembly_accession
cat raw.tsv |
    tsv-select -H -f "assembly_accession" \
    > rs.acc.tsv
```

- genbank 中筛选上一步骤中得到的species信息

```shell
cd ~/project/nwr/Actionmycetes/summary
SPECIES=$(
    cat GB1.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

# 排除未知种和杂交种，并排除refseq中出现的
echo "
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% x %'
    " |
    sqlite3 -separator $'\t' ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

# 添加未知种
echo "
    SELECT
        genus || ' sp. ' || infraspecific_name || ' ' || assembly_accession AS name,
        genus || ' sp.', genus, ftp_path, biosample, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species LIKE '% sp.%'
    " |
    sqlite3 -separator $'\t' ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

# 统计、检查，保证每行字段相同
cat raw.tsv |
    rgr dedup stdin |
    datamash check
#16374 lines, 7 fields
```

- 创建基于目标基因组的表格

```shell
cd ~/project/nwr/Actionmycetes/summary
# 物种名称标准化及缩写
cat raw.tsv |
    grep -v '^#' |
    rgr dedup stdin |
    tsv-select -f 1-6 |
    perl ../abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\tbiosample\tspecies\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[6]}++; # abbr_name
        $seen{$F[6]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\t%s\n}, $F[6], $F[3], $F[4], $F[1], $F[5];
        ' |
    tsv-filter --or --str-in-fld 2:ftp --str-in-fld 2:http |
    keep-header -- tsv-sort -k4,4 -k1,1 \
    > Actionmycetes.assembly.tsv

# 统计检查
datamash check < Actionmycetes.assembly.tsv
# 16356 lines, 5 fields

# 去重
cat Actionmycetes.assembly.tsv |
    tsv-uniq -f 1 --repeated

# 删除无ftp链接的
cat Actionmycetes.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp
```

### 3.2 正式下载之前检查计数

在下载基因组数据refseq和genbank之前，进行检查计数。分类水平及统计，genus属水平数量统计

```shell
cd ~/project/nwr/Actionmycetes

# 准备下载用到的脚本，生成三个sh脚本和一个species.tsv
nwr template summary/Actionmycetes.assembly.tsv \
    --count \
    --rank genus

# 运行strains.sh，生成菌株的分类tsv文件、数量统计tsv文件
bash Count/strains.sh

# 数量统计tsv文件转化为markdown格式
cat Count/taxa.tsv |
    rgr md stdin --fmt

# 运行rank.sh，生成genus count文件和genus list文件
bash Count/rank.sh
mv Count/genus.count.tsv Count/genus.before.tsv

# 转化为markdown格式
cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:100 |
    rgr md stdin --num
```
`Count/taxa.tsv`文件内容：
| item    |  count |
| ------- | -----: |
| strain  | 16,349 |
| species |  2,032 |
| genus   |    145 |
| family  |     10 |
| order   |      5 |
| class   |      1 |

`Count/genus.before.tsv`部分内容：
| genus             | #species | #strains |
| ----------------- | -------: | -------: |
| Actinomadura      |       72 |      171 |
| Actinomyces       |       48 |      487 |
| Actinoplanes      |       36 |      152 |
| Amycolatopsis     |       88 |      293 |
| Kitasatospora     |       35 |      395 |
| Micromonospora    |      123 |      846 |
| Nocardiopsis      |       47 |      156 |
| Nonomuraea        |       66 |      287 |
| Pseudonocardia    |       53 |      209 |
| Schaalia          |       13 |      107 |
| Streptomyces      |      816 |    11289 |
| Streptosporangium |       28 |      110 |
| Trueperella       |        7 |      142 |

### 3.3 下载并检查

```shell
cd ~/project/nwr/Actionmycetes

# 生成下载所需文件/脚本
nwr template summary/Actionmycetes.assembly.tsv \
    --ass

# 下载
bash ASSEMBLY/rsync.sh
# 由于后续antismash需要用到gbff注释文件，将rsync.sh中--exclude="\*\_genomic.gbff.gz"删除掉

# 检查
bash ASSEMBLY/check.sh

# 对n50等信息进行检查。生成两个文件（n50.tsv and n50.pass.tsv）
bash ASSEMBLY/n50.sh 50000 500 500000
# n50.sh需要提前安装hnsm


# Actinot_timonense_UMB9819B_GCF_030227995_1没有基因组数据，NCBI上只有report.txt、stats.txt、md5.txt三个文件，无法计算n50，故在n50.tsv中删除这些
# Kib_sp_MJ126_NF4_GCF_000826545_1
# Microm_radicis_AZ1_13_GCF_003583405_1
# Streptos_coryd_CCUG_49571_GCF_042658485_1
# Streptos_feng_CCUG_58335_GCF_042658525_1
# Streptos_jiao_CGMCC_4_1884_GCF_042654945_1
# Streptos_kro_CCUG_49563_GCF_042654805_1
# Streptos_sheng_CGMCC_1_16303_GCF_042655985_1
# Streptos_tarax_CGMCC_4_1896_GCF_042663065_1

cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50" --max "C" --min "S"
# N50_min	C_max	S_min
# 5021	1959	1335511

# Actinomyces_isr_NCTC12972_GCA_900637225_1同上，无法计算n50
# Gl_europaea_UMB10119A_GCA_030230005_1
# Streptomy_sp_Soil_LB_1_GCA_041879315_1

# 然后再手动筛选，生成n50.pass.tsv
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --ge "N50:100000" |
    tsv-filter -H --ge "S:1000000" |
    tsv-filter -H --le "C:1000" \
    > ASSEMBLY/n50.pass.tsv

cat ASSEMBLY/n50.tsv |         
    tsv-summarize -H --quantile "N50:0.1,0.5" --quantile "C:0.5,0.9" --quantile "S:0.1,0.5" |
    datamash transpose
# N50_pct10	30608.2
# N50_pct50	264239
# C_pct50	76
# C_pct90	435
# S_pct10	5847056.8
# S_pct50	8284047

# collect收集每个sample相关信息。生成collect.tsv
bash ASSEMBLY/collect.sh

# finish检验筛选并保存。生成
bash ASSEMBLY/finish.sh


cp ASSEMBLY/collect.pass.tsv summary/
cat ASSEMBLY/counts.tsv |
    rgr md stdin --fmt
```

| #item            | fields |  lines |
| ---------------- | -----: | -----: |
| url.tsv          |      3 | 16,355 |
| check.lst        |      1 | 16,355 |
| collect.tsv      |     20 | 16,356 |
| n50.tsv          |      4 | 16,344 |
| n50.pass.tsv     |      4 | 12,532 |
| collect.pass.tsv |     23 | 12,532 |
| pass.lst         |      1 | 12,531 |
| omit.lst         |      1 |    512 |
| rep.lst          |      1 |  1,601 |
| sp.lst           |      1 |  5,581 |

ps：
> `rsync.sh`从NCBI上进行下载数据
>
> `check.sh`检查下载情况，对下载的文件进行检验`md5`码
>
> `n50.sh`进行n50等组装质量的检验。 后面的三个数字分别是 LEN_N50 表示 N50 值 ， N_CONTIG 表示 contig 的数量 ， LEN_SUM 表示序列总长度，即所有序列长度的总和，just在计算之前作为一个参考
> 
> `collect.sh`对每个sample的相关信息进行收集
>
> `finish.sh`对下载好的各sample文件的质量的检验和筛选，并生成相关文件来保存这些结果。生成counts.tsv、collect.pass.tsv和4个lst文件，`omit.lst`无蛋白序列/编码区序列、`collect.pss.tsv`筛选colletc.tsv中n50值pass的菌株并给无蛋白序列的菌株注释栏标记为No、`rep.lst`筛选collect.pss.tsv中有refseq的菌株、`sp.lst`筛选collect.pss.tsv中没有species taxon的菌株、`counts.tsv`文件行数（菌株数）计数


## 4 Minhash

```shell
cd ~/project/Actionmycetes/genomes

# 生成 minhash 文件和脚本，一个 species.tsv 和 五个sh脚本
nwr template summary/Actionmycetes.assembly.tsv \
    --mh \
    --parallel 8 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.005

# 生成序列的sketch
bash MinHash/compute.sh

# 去重
bash MinHash/nr.sh

find MinHash -name "NR.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/NR.lst
find MinHash -name "redundant.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/redundant.lst
wc -l summary/NR.lst summary/redundant.lst
# 8896 summary/NR.lst
# 2737 summary/redundant.lst

# 异常菌株
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
# 244

# 排除untaxon的
cat MinHash/abnormal.lst |
    tsv-join -e -f ASSEMBLY/sp.lst \
    > MinHash/tmp.lst
mv MinHash/tmp.lst MinHash/abnormal.lst
cat MinHash/abnormal.lst | wc -l
# 224

# 按照新条件重新 minhash
nwr template summary/Actionmycetes.assembly.tsv \
    --mh \
    --parallel 8 \
    --in ASSEMBLY/rep.lst \
    --not-in ASSEMBLY/sp.lst \
    --not-in MinHash/abnormal.lst \
    --not-in summary/redundant.lst \
    --height 0.4
bash MinHash/dist.sh
```

ps：
> `compute.sh`：将fasta序列按照kmer 21用一致的哈希函数转换为hash值，筛选最小的10000生成sketch来代表这个序列
>
> `nr.sh`：序列去重
>
> `abnormal.sh`：筛选出基因组序列的相似度不正常的部分


## 5 重新count

```shell
cd ~/project/Actionmycetes/genomes

# 按照 pass 和 minhash 之后的新条件重新 count
nwr template summary/Actionmycetes.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --rank family --rank genus \
    --lineage family --lineage genus

# 菌株统计，strains.taxon.tsv and taxa.tsv
bash Count/strains.sh
cat Count/taxa.tsv |
    rgr md stdin --num

# rank 分类
bash Count/rank.sh
cat Count/family.count.tsv |
    tsv-filter -H --ge "3:20" |
    rgr md stdin --num

cp Count/strains.taxon.tsv summary/genome.taxon.tsv
```
ps: Allon_opa_DSM_45601_GCF_003002095_1-Allonocardiopsis属，无科分类，Streptosporangiales目

| item    | count |
| ------- | ----: |
| strain  | 12531 |
| species |  1783 |
| genus   |   140 |
| family  |    10 |
| order   |     5 |
| class   |     1 |

| family               | #species | #strains |
| -------------------- | -------: | -------: |
| Actinomycetaceae     |      121 |      633 |
| Micromonosporaceae   |      258 |     1024 |
| Nocardiopsidaceae    |       58 |      147 |
| Pseudonocardiaceae   |      322 |      705 |
| Streptomycetaceae    |      793 |     9355 |
| Streptosporangiaceae |      142 |      489 |
| Thermomonosporaceae  |       86 |      174 |
| Carbonactinosporaceae|        1 |        3 |
| Thermomonosporaceae  |        1 |        1 |
