#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 8 --mem 32gb
module load fasta
pushd results
tail -n +2 Tomato_matchAligned.min100.tsv | cut -f1,2 | perl -p -e 's/(\S+)\s+(\d+)/>$1.$2\n$1/' > Tomato.min100.fa
tail -n +2 Botrytis_matchAligned.min100.tsv  | cut -f1,2 | perl -p -e 's/(\S+)\s+(\d+)/>$1.$2\n$1/' > Botrytis.min100.fa

ssearch36 -E 0.05 -m 8c Tomato.min100.fa ../mirDB/sly_mature  > Tomato.min100.sly_search.SSEARCH.tab
ssearch36 -E 0.05 -m 8c Tomato.min100.fa ../mirDB/mature.fa > Tomato.min100.mirbase.SSEARCH.tab
ssearch36 -E 0.05 -m 8c Botrytis.min500.fa ../mirDB/sly_mature  > Botrytis.min100.sly_search.SSEARCH.tab
ssearch36 -E 0.05 -m 8c Botrytis.min500.fa ../mirDB/mature.fa  > Botrytis.min100.mirbase.SSEARCH.tab
