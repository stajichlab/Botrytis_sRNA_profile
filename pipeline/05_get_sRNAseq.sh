#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 1 --mem 128gb --out logs/seq_sRNA_BAM.%a.log -a 1,2

module load biopython
module load workspace/scratch
N=$SLURM_ARRAY_ID
if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
	N=1
    fi
fi
BAMFILE=$(ls results/*.bam | sed -n ${N}p)

time ./scripts/get_BAM_sRNA_unique.py -b $BAMFILE --sra2name sra2name.tsv 

