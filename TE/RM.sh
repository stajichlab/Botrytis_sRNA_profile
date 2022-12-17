#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 128 --mem 64gb 

module load RepeatMasker
CPU=64
RepeatMasker -pa $CPU -species fungi -e rmblast -gff -dir BcinereaB05-10.RM_rmblast  FungiDB-60_BcinereaB05-10_Genome.fasta > FungiDB-60_BcinereaB05-10_Genome.RM_rmblast.log
