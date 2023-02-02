#!/usr/bin/bash -l
#SBATCH -p batch -N 1 -c 24 --mem 64gb --out logs/trim.%A.log

module load fastp
module load parallel
parallel -j 4 if \[ ! -f trimmed/{/.}.gz \]\; then fastp -i {} -o trimmed/{/.}.gz -w 6 -p\; fi ::: $(ls input/*.gz)
