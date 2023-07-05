#!/bin/bash
#SBATCH --time=5
#SBATCH --nodes=1
#SBATCH --mem=1GB
# use the current directory
#$ -cwd
#$ -S /bin/bash

#srun /home/users/habnice/R/bin/Rscript test.R
# srun module load R
# /usr/bin/Rscript $*
srun /home/users/habnice/.local/bin/Rscript $*


## submit command
# for n in $(seq 30 50 230); do for b in {1..10}; do sbatch ./submit.sh sim.R $n; done; done
## for b in {1..10}; do for n in $(seq 1e3 1e3 8e3); do for pi0 in .05 .03; do sbatch ./submit.sh sim.R $n 10 8 7 $pi0; done; done; done&
