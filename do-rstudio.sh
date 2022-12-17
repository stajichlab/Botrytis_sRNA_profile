#!/usr/bin/bash
sbatch -J rstudio -p short -c 8 --mem=12g --wrap='module unload R; module load R/4.2.0; module load rstudio-server/2022.02.0-443; start-rserver_4.2.0.sh'  --output='rstudio-%J.out'
