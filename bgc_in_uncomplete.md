

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

# Build pHMM database
cat 2refer/2hmm/3phmm/*/*.hmm > 2refer/all.hmm
hmmpress 2refer/all.hmm

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

for i in $(cat search_threshold/strain_75_new.tsv); do
    cat search_threshold/threshold_m/${i}.tsv.ratios.tsv | grep "${i}" >> search_threshold/threshold_m/${i}_tmp.tsv
    for domain in A C T E TE X nMT; do
        cat search_threshold/threshold_m/${i}_tmp.tsv | tsv-filter --str-eq 1:"${domain}" >> search_threshold/threshold_m/${domain}_ratio.tsv
    done
    rm search_threshold/threshold_m/${i}_tmp.tsv
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
    cat search_threshold/threshold_qc/${i}.tsv.ratios.tsv | grep "${i}" >> search_threshold/threshold_qc/${i}_tmp.tsv
    for domain in A C T E TE X nMT; do
        cat search_threshold/threshold_qc/${i}_tmp.tsv | tsv-filter --str-eq 1:"${domain}" >> search_threshold/threshold_qc/${domain}_ratio.tsv
    done
    rm search_threshold/threshold_qc/${i}_tmp.tsv
done


# Handle blast results: 
# remove outliers, matches with low confidence, matches with low quality, and the combine matches
mkdir -p 3search/2blast/4filter 3search/2blast/5site
cat strains.all.tsv | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 --tagstring "TASK {1}" '
    echo "=> remove outliers."
    python 0script/remove_outliers.py -m blast -i 3search/2blast/3result/{1}.tsv -o 3search/2blast/4filter/{1}.outliers.tsv

    echo "=> remove matches with low confidence."
    cat 3search/2blast/4filter/{1}.outliers.tsv |
        tsv-filter --ge 9:65 --le 10:1.0e-20 --ge 11:50 |
        sort -k4,4 -k6,6n -k7,7n -k8,8nr |
        uniq > 3search/2blast/4filter/{1}.confidence.tsv

    echo "=> remove matches with low quality."
    python 0script/quality_control.py 3search/2blast/4filter/{1}.confidence.tsv 3search/2blast/4filter/{1}.quality.tsv

    echo "=> combine matches."
    python 0script/combine_match.py -t 0.95 -i 3search/2blast/4filter/{1}.quality.tsv -o 3search/2blast/5site/{1}.tsv
'

```


## Extract Matched Domains Sequences

```shell
cd ~/project/Actionmycetes/search

# Merge search results of miniprot and blast
mkdir -p 3search/3site/1raw
cat strains.all.tsv | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 '
    cat 3search/1miniprot/5site/{1}.tsv 3search/2blast/5site/{1}.tsv |
        awk -F"\t" "{ if (\$5 == 0) { \$5 = 1; } print \$0; }" OFS="\t" |
        sort -k3,3 -k5,5n -k6,6n |
        uniq > 3search/3site/1raw/{1}.tsv
'

# Merge site information and divide the cluster
mkdir -p 3search/3site/2site 3search/3site/3cluster
cat strains.all.tsv | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 '
    python 0script/recombine_site.py -t 0.95 -i 3search/3site/1raw/{1}.tsv -o 3search/3site/2site/{1}.tsv
    python 0script/divide_cluster.py 3search/3site/2site/{1}.tsv 3search/3site/3cluster/{1}.tsv
'

# Extract site information（contig and location） of matched Domains
mkdir -p 3search/4subseq/1region
cat strains.all.tsv | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 '
    cat 3search/3site/3cluster/{1}.tsv |
        grep -v ">" |
        cut -f 7,9,10 |
        perl -nla -e "
            my (\$contig, \$start, \$end) = split /\t/;
            print qq{\$contig:\$start-\$end};
        " >3search/4subseq/1region/{1}.tsv
'

# Extract sequences of matched Domains
cat strains.all.tsv | while IFS=$'\t' read -r strain genus id _; do
    mkdir -p 3search/4subseq/2nucl/${strain}
    gzip -dk ASSEMBLY/${genus}/${strain}/${id}_genomic.fna.gz
    cat 3search/4subseq/1region/${strain}.tsv | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 8 "
        faops region -l 0 ASSEMBLY/${genus}/${strain}/${id}_genomic.fna <(echo {}) 3search/4subseq/2nucl/${strain}/{}.fa
    "
    rm -f ASSEMBLY/${genus}/${strain}/${id}_genomic.fna
done

# Six frame translation (dna->aa) 
cat strains.all.tsv | while IFS=$'\t' read -r strain genus id _; do
    mkdir -p 3search/4subseq/3pro/${strain}
    cat 3search/4subseq/1region/${strain}.tsv | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 8 "
        perl 0script/six_frame_translation.pl 3search/4subseq/2nucl/${strain}/{}.fa 3search/4subseq/3pro/${strain}/{}.fa
    "
done

```

## Determination of Domain Types

```shell
cd ~/project/Actionmycetes/search

# hmmscan based on pHMMs
cat strains.all.tsv | while IFS=$'\t' read -r strain genus id _; do
    mkdir -p 3search/5hmmscan/1scan/${strain}
    cat 3search/4subseq/1region/${strain}.tsv | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 2 "
            hmmscan --cpu 4 --tblout 3search/5hmmscan/1scan/${strain}/{}.txt 2refer/all.hmm 3search/4subseq/3pro/${strain}/{}.fa
        "
done

# Extract hmmscan results
cat strains.all.tsv | while IFS=$'\t' read -r strain genus id _; do
    mkdir -p 3search/5hmmscan/2result/${strain}
    cat 3search/4subseq/1region/${strain}.tsv | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 8 "
            python 0script/extract_hmm.py 3search/5hmmscan/1scan/${strain}/{}.txt 3search/5hmmscan/2result/${strain}/{}.tsv
        "
done

# Merge different results reflecting the same site
cat strains.all.tsv | while IFS=$'\t' read -r strain genus id _; do
    mkdir -p 3search/5hmmscan/3merge/${strain}
    cat 3search/4subseq/1region/${strain}.tsv | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 8 "
        if [ -s 3search/5hmmscan/2result/${strain}/{}.tsv ]; then
            python 0script/merge_lines.py -i 3search/5hmmscan/2result/${strain}/{}.tsv -o 3search/5hmmscan/3merge/${strain}/{}.tsv
        fi
    "
done 

# Merge hmmscan results
mkdir -p 3search/5hmmscan/4combine
cat strains.all.tsv | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
    if [ "$(ls -A 3search/5hmmscan/3merge/{1}/ 2>/dev/null)" ]; then
        cat 3search/5hmmscan/3merge/{1}/*.tsv > 3search/5hmmscan/4combine/{1}.tsv
    else
        touch 3search/5hmmscan/4combine/{1}.tsv
    fi
'

# Merge search and hmmscan results
mkdir -p 3search/6all
cat strains.all.tsv | parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 '
    python 0script/combine_result.py -s 3search/3site/3cluster/{1}.tsv -m 3search/5hmmscan/4combine/{1}.tsv -o 3search/6all/{1}.tsv
'

```



```shell
# find threshold for 'extract_hmm.py'
cd project/Actionmycetes/search/


faops split-name -l 0 2refer/fixed_represent_domain.all.fa 2refer/represent_domain
cat 2refer/fixed_represent_domain.all.fa | grep ">" | sed 's/>//g' > hmm_threshold/fasta.lst
cp 2refer/represent_domain hmm_threshold/

cat hmm_threshold/fasta.lst | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 4 "
    hmmscan --cpu 4 --tblout hmm_threshold/result/{}.txt 2refer/all.hmm hmm_threshold/represent_domain/{}.fa
"

for i in $(cat hmm_threshold/fasta.lst); do
    cat hmm_threshold/result/${i}.txt | 
        perl -nla -e '
            /#/ and next;
            my @p = split /\./, $F[2];
            print if $F[0] eq join(".", @p[0..2]);
        ' \
    >> hmm_threshold/result.tsv
done

```

```shell
# find threshold for 'divide_cluster.py'
cat ok_domain_aa_75.txt | grep -A 1 -E "Condensation|Cglyc|AMP-binding|-PCP|nMT|-X-|Epimerization|Thioesterase" | sed '/^--$/d' \
    >> domain_75.txt


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
    ' format.tsv domain_75.txt > fixed_domain_75.fa

awk '
BEGIN{
    while((getline < "remove.lst")>0) del[$1]=1
}
/^>/{
    name=substr($0,2)
    keep = !(name in del)
}
keep
' fixed_domain_75.fa > domain_75.fa

cat domain_75.fa | grep '>' | sed 's/>//g' > domain_all.tsv

awk '{
    split($0, a, ".")
    domain = a[1]

    split($0, b, "-")
    n = length(b)
    start = b[n-3]
    end   = b[n-2]

    print $0 "\t" domain "\t" start "\t" end
}' domain_all.tsv > distance.tsv

cd antismash/antismash_summary/product/antismash_all


python ../../gene_location_from_gbk.py ../../aa1/nrps_dna/nrps_dna_identifier.tsv gene_location.tsv


        cat ../knownGPA_mibig_list.txt | grep -v '^$' |
        while IFS=$'\t' read -r strain cluster _ _ _ _; do
            region=$(echo ${cluster} | sed -E 's/c[0-9]+//g' | sed -E 's/r//g')
            clu=$(echo ${cluster} | sed -E 's/r[0-9]+//g' | sed -E 's/c//g')
			json_path=$(ls ${strain}/*.json 2>/dev/null | head -n 1)
            cds=$(python ../../antismash_json_nrps_identifier.py ${json_path} ${region} ${clu})
            identifier=$(echo $cds | sed -E "s/\[//g "| sed "s/\]//g" | sed "s/'//g"  | sed 's/, /\n/g')
            echo ${identifier} | sed "s#^#${strain}\tr${region}c${clu}\t#g" \
            >> nrps_dna_identifier.tsv
        done

python ../../gene_location_from_gbk.py nrps_dna_identifier.tsv gene_location.tsv
```




