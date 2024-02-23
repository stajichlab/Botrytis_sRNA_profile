#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 4 --mem 24gb --out logs/size_profile_BAM.log

module load biopython
module load parallel

parallel -j 2 ./scripts/process_BAM_sRNA.py -b {} --sra2name sra2name.tsv ::: $(ls results/*.bam)
