#!/usr/bin/bash -l
tail -n +4 BcinereaB05-10.RM_rmblast/FungiDB-60_BcinereaB05-10_Genome.fasta.out | \
	awk 'BEGIN{OFS="\t"} {print $5,"RepeatMasker","match",$6,$7,".",$9,".","type="$11";family="$10}' | \
	perl -p -e 's/\tC\t/\t-\t/' > BcinereaB05-10.RM.gff3
