#!/usr/bin/bash -l
module load bcftools
module load genometools
DIR=genomes/botrytis
IN=$DIR/FungiDB-60_BcinereaB05-10.gff
RPT=TE/BcinereaB05-10.RM.gff3
DATFILE=$DIR/BcinereaB05-10.features.gff3
FINAL=$DIR/BcinereaB05-10.features.sort.gff3
cat $IN $RPT | grep -v -P "\t(gene|mRNA|CDS|three_prime_UTR|five_prime_UTR|protein_coding_gene|ncRNA_gene)\t"  > $DATFILE
gt gff3 -addintrons $IN | grep -P "\tintron\t" >> $DATFILE
(grep ^"#" $DATFILE; grep -v ^"#" $DATFILE | grep -v "^$" | grep "\t" | sort -k1,1 -k4,4n) > $FINAL
bgzip -f $FINAL
tabix $FINAL.gz
