#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 2 --mem 32gb --out logs/countreads_byexp_subtractTomato.log

module load biopython

time ./scripts/process_reads_byexp_subtract_Tomato.py
Rscript scripts/plot_sRNA_sizes_uniq.R
Rscript scripts/plot_sRNA_sizes.R
