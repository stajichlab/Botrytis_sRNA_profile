#!/usr/bin/bash -l
module load genometools
gt gff3 -addintrons FungiDB-60_BcinereaB05-10.gff | grep intron > introns.gff3
