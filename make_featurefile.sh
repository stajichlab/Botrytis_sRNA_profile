#!/usr/bin/bash -l
module load bcftools
module load genometools
DATFILE=BcinereaB05-10.features.gff3
FINAL=BcinereaB05-10.features.sort.gff3
cat genomes/botrytis/FungiDB-60_BcinereaB05-10.gff TE/BcinereaB05-10.RM.gff3 | grep -v -P "\t(gene|mRNA|CDS|three_prime_UTR|five_prime_UTR|protein_coding_gene|ncRNA_gene)\t"  > $DATFILE
gt gff3 -addintrons genomes/botrytis/FungiDB-60_BcinereaB05-10.gff  | grep intron >> $DATFILE
(grep ^"#" $DATFILE; grep -v ^"#" $DATFILE | grep -v "^$" | grep "\t" | sort -k1,1 -k4,4n) > $FINAL

bgzip -f $FINAL
tabix $FINAL.gz
