#!/usr/bin/bash -l
#SBATCH -p short -N 1 -c 16 -C xeon --mem 64gb

module load star
GENOME=FungiDB-60_BcinereaB05-10_Genome.fasta
GFF=FungiDB-60_BcinereaB05-10.gff
STAR --runThreadN 16 --runMode genomeGenerate \
	--genomeFastaFiles $GENOME --sjdbGTFfile $GFF \
	--sjdbOverhang 100
