#!/usr/bin/bash -l
#SBATCH -p batch -N 1 -c 24 --mem 64gb --out logs/trim.%A.log

module load fastp
module load parallel
parallel -j 4 fastp -i {} -o trimmed/{/.}.gz -w 6 -p ::: $(ls input/*.gz)

