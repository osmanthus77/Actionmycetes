

## tools installation

```shell
# install by homebrew
brew install cd-hit
brew install hmmer
brew info blast
brew install easel
brew info miniprot

```

## represent sequence and pHMMs of domains

### representative sequences obtained by clustering

```shell
cd ~/project/Actionmycetes/search
mkdir -p 2refer/1domain/1raw 2refer/1domain/2cdhit

# sequence
for domain in A C T E TE X nMT; do
    cp ../analysis_75/sequence_aa/domain_${domain}aa_75.fa 2refer/1domain/1raw
done

# sequence identity by diamond
diamond makedb --in 2refer/1domain/1raw/domain_nMTaa_75.fa -d 2refer/1domain/diamond/nMTdomain_db
diamond blastp -q 2refer/1domain/1raw/domain_nMTaa_75.fa -d 2refer/1domain/diamond/nMTdomain_db.dmnd \
    -o 2refer/1domain/diamond/nMT_identity.tsv \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --sensitive
python 2refer/1domain/diamond/plot_identity.py 2refer/1domain/diamond/nMT_identity.tsv


# clustering
cd-hit -i 2refer/1domain/1raw/domain_Aaa_75.fa -o 2refer/1domain/2cdhit/domain_Aaa_75.fa -c 0.95 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/domain_Caa_75.fa -o 2refer/1domain/2cdhit/domain_Caa_75.fa -c 0.96 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/domain_Taa_75.fa -o 2refer/1domain/2cdhit/domain_Taa_75.fa -c 0.96 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/domain_Eaa_75.fa -o 2refer/1domain/2cdhit/domain_Eaa_75.fa -c 0.93 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/domain_TEaa_75.fa -o 2refer/1domain/2cdhit/domain_TEaa_75.fa -c 0.95 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/domain_Xaa_75.fa -o 2refer/1domain/2cdhit/domain_Xaa_75.fa -c 0.96 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/domain_nMTaa_75.fa -o 2refer/1domain/2cdhit/domain_nMTaa_75.fa -c 0.98 -n 5 -M 10000 -T 8 -d 0

for domain in A C T E TE X nMT; do
    cat 2refer/1domain/2cdhit/domain_${domain}aa_75.fa >> 2refer/represent_domain.all.fa
done

```

### build pHMMs models

```shell
cd ~/project/Actionmycetes/search
mkdir -p 2refer/2hmm/1raw

# sequence
for domain in A C T E TE X nMT; do
    cp 2refer/1domain/1raw/domain_${domain}aa_75.fa 2refer/2hmm/1raw
done

cd 2refer/2hmm/1raw
# copy the annotate file here and modify its format like below（format.tsv for each domain）
#C	Cglyc	C_N2	Actinok_aurantic_DSM_44650_GCF_041897805_1-r1c3-actinoidinB-Cglyc-593-886-2.3-active
#C	Cglyc	C_N1	Actinok_aurantic_DSM_44650_GCF_041897805_1-r1c3-actinoidinB-Cglyc-12-309-3.1-active

# modify the format of fasta sequence name
for domain in A C T E TE X nMT; do
    awk '
    NR==FNR {
        map[$4] = $1"."$2"."$3"."
        next
    }
    /^>/ {
        name = substr($0, 2)
        if (name in map) {
            print ">" map[name] name
        } else {
            print "NOT FOUND: " name > "/dev/stderr"
            print $0
        }
        next
    }
    {print}
    ' ${domain}domain_format_75.tsv domain_${domain}aa_75.fa > fixed_domain_${domain}aa_75.fa
done

cd ~/project/Actionmycetes/search
# deduplicate
for domain in A C T E TE X nMT ; do
    mkdir -p 2refer/2hmm/2class/${domain}domain

    # Obtain domain type list
    cat 2refer/2hmm/1raw/fixed_domain_${domain}aa_75.fa |
        grep '>' |
        sed 's/>//g' |
        cut -d'.' -f 1,2,3 |
        sort | uniq > 2refer/2hmm/2class/${domain}domain.tsv

    # Classify domain sequences
    cat 2refer/2hmm/2class/${domain}domain.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 8 "
            cat 2refer/2hmm/1raw/fixed_domain_${domain}aa_75.fa |
            grep -A 1 {} |
            sed 's/^--$//g' |
            grep -v '^$'  > 2refer/2hmm/2class/${domain}domain/{}.fa
        "

    # deduplicate
    cat 2refer/2hmm/2class/${domain}domain.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 8 "
            cd-hit -i 2refer/2hmm/2class/${domain}domain/{}.fa \
            -o 2refer/2hmm/2class/${domain}domain/{}.de.fa -c 1.00 -n 5 -M 10000 -T 8 -d 0
        "
done

# build pHMM model
for domain in A C T E TE X nMT ; do
    mkdir -p 2refer/2hmm/3phmm/${domain}domain

    cat 2refer/2hmm/2class/${domain}domain.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 8 "
            mafft --auto 2refer/2hmm/2class/${domain}domain/{}.de.fa >2refer/2hmm/3phmm/${domain}domain/{}.mafft.fa
            trimal -in 2refer/2hmm/3phmm/${domain}domain/{}.mafft.fa -out 2refer/2hmm/3phmm/${domain}domain/{}.trim.fa
            esl-reformat stockholm 2refer/2hmm/3phmm/${domain}domain/{}.trim.fa >2refer/2hmm/3phmm/${domain}domain/{}.sto
            hmmbuild 2refer/2hmm/3phmm/${domain}domain/{}.hmm 2refer/2hmm/3phmm/${domain}domain/{}.sto
    "
done




```

## representative sequences search

### search by miniprot

```shell
cd ~/project/Actionmycetes/antismash/group
cat Act*pass.tsv | tr " " "_" | 
    sed "s/sp\./sp/g" |
    sed "s/(//g" |
    sed "s/)//g" |
    sed "s/\.//g" |
    sed "s/\[//g" |
    sed "s/\]//g" \
    >> strains_all.tsv
cp strains_all.tsv ~/project/Actionmycetes/search

cd ~/project/Actionmycetes/genomes/summary 
awk -F'\t' '{split($2, a, "/"); OFS="\t"; print $1, a[10]}' Actionmycetes.assembly.tsv > tmp.tsv
cp tmp.tsv ~/project/Actionmycetes/search
rm tmp.tsv


cd ~/project/Actionmycetes/search
ln -s ~/project/Actionmycetes/genomes/ASSEMBLY ASSEMBLY

tsv-join --filter-file tmp.tsv --key-fields 1 --append-fields 2 strains_all.tsv > strains.all.tsv


cat group/Act_uncomplete_untaxon_pass.tsv | tsv-join -f group/family_GCF/${family}_GCF.tsv -k 1 > group/family_result/${family}_pass_uncomplete_untaxon.tsv

# miniprot index and search
mkdir -p 3search/1miniprot/1index 3search/1miniprot/2search
mkdir strain_split
split -l 20 strains.all.tsv strain_split/strain_
ls strain_split/* | sed -E 's!strain_split/!!g' > strain_split.tsv

for i in $(cat strain_split.tsv); do
    echo "==> ${i}"
    cat strain_split/${i} | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo "=> genome {1} is in index."
        miniprot -t 2 -d 3search/1miniprot/1index/{1}.mpi ASSEMBLY/{2}/{1}/{3}_genomic.fna.gz 
        echo "=> genome {1} is in search."
        miniprot -t 2 --outs=0.85 --gff --aln --trans 3search/1miniprot/1index/{1}.mpi \
            2refer/represent_domain.all.fa > 3search/1miniprot/2search/{1}.gff3
'
    rm -f 3search/1miniprot/1index/*.mpi
done

# Extract miniprot search results
mkdir -p 3search/1miniprot/3result
cat strains.all.tsv | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 12 '
    cat 3search/1miniprot/2search/{1}.gff3 | 
        grep "##PAF" | 
        cut -f 2,3,6,7,8,9,10,15 |
        sed "s/ms:i://g" |
        sort -k4,4 -k6,6n -k7,7n -k8,8nr |
        uniq > 3search/1miniprot/3result/{1}.tsv
'

# modify the format of domain_id(first field)
cat 2refer/2hmm/1raw/*.tsv >> format.tsv
cat strains.all.tsv | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 12 '
    perl 0script/domain_name_modify.pl format.tsv 3search/1miniprot/3result/{1}.tsv > 3search/1miniprot/3result/{1}_tmp.tsv \
        && mv 3search/1miniprot/3result/{1}_tmp.tsv 3search/1miniprot/3result/{1}.tsv
'


# find the threshold for outlier ('0script/remove_outliers.py')(hit length/domain_length)
mkdir -p search_threshold/miniprot search_threshold/blast search_threshold/threshold_b search_threshold/threshold_m
touch search_threshold/strain_75.tsv
for i in $(cat search_threshold/strain_75.tsv); do
    cp 3search/1miniprot/3result/${i}.tsv search_threshold/miniprot
    cp 3search/2blast/3result/${i}.tsv search_threshold/blast
done

for i in $(cat search_threshold/strain_75_new.tsv); do
    python 0script/calc_ratio_outlier.py -i search_threshold/miniprot/${i} -o search_threshold/threshold_m/${i}
    python 0script/calc_ratio_outlier.py -i search_threshold/blast/${i} -o search_threshold/threshold_b/${i}
done





# Handle miniprot results: remove outliers and merge identical matches
mkdir -p 3search/1miniprot/4filter 3search/1miniprot/5site
cat strains.all.tsv | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 --tagstring "TASK {1}" '
    echo "=> remove outliers."
    python 0script/remove_outliers.py -m miniprot -i 3search/1miniprot/3result/{1}.tsv -o 3search/1miniprot/4filter/{1}.outliers.tsv

    echo "=> combine match."
    python 0script/combine_match.py -t 0.95 -i 3search/1miniprot/4filter/{1}.outliers.tsv -o 3search/1miniprot/5site/{1}.tsv
'

```

### search by blast

```shell
cd ~/project/Actionmycetes/search

# create database and search
mkdir -p 3search/2blast/2search

for i in $(cat strain_split_new.tsv); do
    echo "==> ${i}"
    cat strain_split/${i} | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 --tagstring "TASK {1}" '
        mkdir -p 3search/2blast/1db/{1}

        echo "=>decompress step."
        gzip -dk ASSEMBLY/{2}/{1}/{3}_genomic.fna.gz

        echo "=>makeblastdb step."
        makeblastdb -in ASSEMBLY/{2}/{1}/{3}_genomic.fna -dbtype nucl -out 3search/2blast/1db/{1}/{1} 1>/dev/null

        echo "=>remove fna file."
        rm -f ASSEMBLY/{2}/{1}/{3}_genomic.fna

        echo "=>search step."
        tblastn -query 2refer/represent_domain.all.fa -db 3search/2blast/1db/{1}/{1} -out 3search/2blast/2search/{1}.tsv \
            -num_threads 4 \
            -outfmt "6 qseqid qlen sstrand sseqid slen sstart send bitscore pident evalue qcovs"

        echo "=>modify blast outfile format."
        sed -i bak -E "s/plus/+/g;s/minus/-/g" 3search/2blast/2search/{1}.tsv
        rm -f 3search/2blast/2search/{1}.tsvbak
'
    rm -rf 3search/2blast/1db/
done


# Extract blast search results and Format results
mkdir -p 3search/2blast/3result
cat strains.all.tsv | parallel --colsep '\t' parallel --no-run-if-empty --linebuffer -k -j 12 '
    perl 0script/blast2miniprot.pl 3search/2blast/2search/{1}.tsv > 3search/2blast/3result/{1}.tsv
'
# modify the format of domain_id(first field)
cat strains.all.tsv | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 12 '
    perl 0script/domain_name_modify.pl format.tsv 3search/2blast/3result/{1}.tsv > 3search/2blast/3result/{1}_tmp.tsv \
        && mv 3search/2blast/3result/{1}_tmp.tsv 3search/2blast/3result/{1}.tsv
'


# find threshold for quality control ('0script/quality_control.py')(bitscore/domain_length)
for i in $(cat search_threshold/strain_75_new.tsv); do
    python 0script/calc_ratio_qc.py -i search_threshold/blast/${i}.tsv -o search_threshold/threshold_qc/${i}
done
for i in $(cat search_threshold/strain_75_new.tsv); do
    cat search_threshold/threshold_qc/${i}.ratios.tsv | grep "${i}" >> search_threshold/threshold_qc_2/${i}.ratios.tsv
done
for i in $(cat search_threshold/strain_75_new.tsv); do
    cat search_threshold/threshold_qc_2/${i}.ratios.tsv | grep 'nMT.nMT' | awk -F'\t' '!seen[$2]++' >> search_threshold/threshold_qc_2/nMT.ratios.tsv
done


# Handle blast results: 
# remove outliers, matches with low confidence, matches with low quality, and the combine matches
mkdir -p 3search/2blast/4filter 3search/2blast/5site
cat strains.all.tsv | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 --tagstring "TASK {1}" '
    echo "=> remove outliers."
    python 0script/remove_outliers.py -m blast -i 3search/2blast/3result/{1}.tsv -o 3search/2blast/4filter/{1}.outliers.tsv

    echo "=> remove matches with low confidence."
    cat 3search/2blast/4filter/{1}.outliers.tsv |
        tsv-filter --ge 9:70 --le 10:1.0e-20 --ge 11:50 |
        sort -k4,4 -k6,6n -k7,7n -k8,8nr |
        uniq > 3search/2blast/4filter/{1}.confidence.tsv

    echo "=> remove matches with low quality."
    python 0script/quality_control.py 3search/2blast/4filter/{1}.confidence.tsv 3search/2blast/4filter/{1}.quality.tsv

    echo "=> combine matches."
    python 0script/combine_match.py -t 0.95 -i 3search/2blast/4filter/{1}.quality.tsv -o 3search/2blast/5site/{1}.tsv
'


```

```shell
# 

cd ~/project/Actionmycetes/search
find 3search/2blast/5site -name "*.tsv" -type f -size +0c | sed -E 's!3search/2blast/5site/!!g; s!.tsv!!g' >> 3search/2blast/site.tsv
find 3search/1miniprot/5site -name "*.tsv" -type f -size +0c | sed -E 's!3search/1miniprot/5site/!!g; s!.tsv!!g' >> 3search/1miniprot/site.tsv

grep -vxFf 3search/2blast/site.tsv  3search/1miniprot/site.tsv > 3search/miniprot_blast_site_diff.tsv
shuf -n 100 3search/miniprot_blast_site_diff.tsv > 3search/antismash_test.tsv
tsv-join --filter-file strains.all.tsv --key-fields 1 --append-fields 2,3 3search/antismash_test.tsv > 3search/antismash.tsv

source ~/miniconda3/bin/activate
conda activate antismash

cat 3search/antismash.tsv |
parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
    echo {1};
    antismash --taxon bacteria -c 4 --cb-general --cc-mibig --cb-knownclusters --pfam2go --asf --genefinding-tool prodigal ASSEMBLY/{2}/{1}/{3}_genomic.fna.gz --output-dir 3search/antismash/{1}
    "


```
