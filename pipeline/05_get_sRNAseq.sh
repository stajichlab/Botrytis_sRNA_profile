#!/usr/bin/bash -l
#SBATCH -p short -C xeon -N 1 -n 1 -c 4 --mem 24gb --out logs/seq_sRNA_BAM.log

module load biopython
module load parallel

parallel -j 2 ./scripts/get_BAM_sRNA_unique.py -b {} --sra2name sra2name.tsv ::: $(ls results/*.bam)
