#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 48 --mem 64gb --out logs/align.%A.log
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
OUTDIR=results
MANIFEST=manifest.tsv
mkdir -p $OUTDIR

if [[ ! -f $MANIEST || $SAMPLES -nt $MANIEST ]]; then
	csvgrep -c 2 -r miRNA $SAMPLES | csvcut -c 1,5  | perl -p -e 's/,/_1.fastq.gz\t-\t/; s/^/trimmed\//' | tail -n +2 > $MANIFEST
fi

if [ ! -s $OUTDIR/Tomato_matchLog.final.out ]; then
	STAR --runThreadN $CPU --genomeDir $tomato --readFilesIn $INDIR --readFilesCommand zcat \
		 --outFileNamePrefix $OUTDIR/Tomato_match --readFilesManifest $MANIFEST \
		 --outTmpDir $SCRATCH/Tom
fi

if [ ! -s $OUTDIR/Botrytis_match.final.out ]; then
	STAR --runThreadN $CPU --genomeDir $botrytis --readFilesIn $INDIR --readFilesCommand zcat \
		 --outFileNamePrefix $OUTDIR/Botrytis_match --readFilesManifest $MANIFEST \
		 --outTmpDir $SCRATCH/Bot
fi
