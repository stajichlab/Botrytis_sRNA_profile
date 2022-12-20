#!/usr/bin/bash -l
module load bcftools
(grep ^"#" BcinereaB05-10.features.gff3; grep -v ^"#" BcinereaB05-10.features.gff3 | grep -v "^$" | grep "\t" | sort -k1,1 -k4,4n) > BcinereaB05-10.features.sort.gff3
bgzip BcinereaB05-10.features.sort.gff3
tabix BcinereaB05-10.features.sort.gff3.gz
