#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 48 --mem 64gb --out logs/align_bylib.%A.log
module load star
module load csvkit
CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

tomato=genomes/tomato/GenomeDir
botrytis=genomes/botrytis/GenomeDir
SAMPLES=samples.csv
INDIR=trimmed
OUTDIR=results_bylib
MANIFEST=manifest.tsv
mkdir -p $OUTDIR

if [[ ! -f $MANIEST || $SAMPLES -nt $MANIEST ]]; then
	csvgrep -c 2 -r miRNA $SAMPLES | csvcut -c 1,5  | perl -p -e 's/,/_1.fastq.gz\t-\t/; s/^/trimmed\//' | tail -n +2 > $MANIFEST
fi
cat $MANIFEST | while read FASTQ BLANK LIBRARY
do
	mkdir -p $OUTDIR/Tomato
	if [ ! -s $OUTDIR/Tomato/${LIBRARY}_Log.final.out ]; then
		STAR --runThreadN $CPU --genomeDir $tomato --readFilesIn $FASTQ --readFilesCommand zcat \
			 --outFileNamePrefix $OUTDIR/Tomato/${LIBRARY}_ --outTmpDir $SCRATCH/Tom
	fi
	mkdir -p $OUTDIR/Botrytis
	if [ ! -s $OUTDIR/Botrytis/${LIBRARY}_Log.final.out ]; then
		STAR --runThreadN $CPU --genomeDir $botrytis --readFilesIn $FASTQ --readFilesCommand zcat \
			 --outFileNamePrefix $OUTDIR/Botrytis/${LIBRARY}_ --outTmpDir $SCRATCH/Bot
	fi
done
