#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 128 --mem 64gb

module load RepeatMasker
CPU=64
RepeatMasker -e hmmer -pa $CPU -species fungi -gff -dir BcinereaB05-10.RM_dfam FungiDB-60_BcinereaB05-10_Genome.fasta > FungiDB-60_BcinereaB05-10_Genome.RM_dfam.log
