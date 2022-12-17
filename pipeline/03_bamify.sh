#!/usr/bin/bash -l
#SBATCH -p short -N 1 -c 16 -n 1 --mem 96gb --out logs/bamify.log

module load samtools
CPU=8
#if [ $SLURM_CPUS_ON_NODE ]; then
#  CPU=$SLURM_CPUS_ON_NODE
#fi

parallel -j 2 samtools sort -O bam --threads $CPU --write-index -T $SCRATCH/{/.} -o {.}.bam {} ::: $(ls results/*.sam)

parallel -j 2 samtools sort -O bam --threads $CPU --write-index -T $SCRATCH/{/.} -o {.}.bam {} ::: $(ls results_bylib/*/*.sam)

rm -f results_bylib/*/*.sam
rm -f results/*.sam
