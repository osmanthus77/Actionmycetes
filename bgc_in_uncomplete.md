

## tools installation

```shell
# install by homebrew
brew install cd-hit
brew install hmmer
brew info blast
brew install easel

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
cd-hit -i 2refer/1domain/1raw/domain_Aaa_75.fa -o 2refer/1domain/2cdhit/domain_Aaa_75_95.fa -c 0.95 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/domain_Caa_75.fa -o 2refer/1domain/2cdhit/domain_Caa_75.fa -c 0.96 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/domain_Taa_75.fa -o 2refer/1domain/2cdhit/domain_Taa_75.fa -c 0.96 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/domain_Eaa_75.fa -o 2refer/1domain/2cdhit/domain_Eaa_75.fa -c 0.93 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/domain_TEaa_75.fa -o 2refer/1domain/2cdhit/domain_TEaa_75.fa -c 0.95 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/domain_Xaa_75.fa -o 2refer/1domain/2cdhit/domain_Xaa_75.fa -c 0.96 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/domain_nMTaa_75.fa -o 2refer/1domain/2cdhit/domain_nMTaa_75.fa -c 0.94 -n 5 -M 10000 -T 8 -d 0

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

