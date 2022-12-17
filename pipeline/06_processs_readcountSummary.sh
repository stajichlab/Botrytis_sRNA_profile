#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 1 --mem 8gb --out logs/countreads_subtractTomato.log

module load biopython

time ./scripts/process_reads_subtract_Tomato.py -s
