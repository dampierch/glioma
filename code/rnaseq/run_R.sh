#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1                       #### set to 20 for BiocParallel
#SBATCH --mem=2000                             #### set to 150 for large BiocParallel
#SBATCH --time=1:00:00                          #### hh:mm:ss
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc                   #### alternative: cphg_caseylab
#SBATCH --mail-type=END
#SBATCH --mail-user=chd5n@virginia.edu

## variable
script_name=$1  ## first positional argument after shell script must be name of R script to run
args="$2 $3 $4 $5 $6 $7"  ## second and more positional arguments after shell script must be args to be passed to Rscript

## header
pwd; hostname; date
## environments
module purge
module load gcc/7.1.0 R/3.6.1
## execution
printf "start ${script_name}\n"
touch ${script_name%.R}.flag  ## to serve as output for snakemake
Rscript ${script_name} --args ${args}
## footer
printf "end ${script_name}\n"; date
echo "seff ${SLURM_JOBID}"
echo "sacct -o reqmem,maxrss,elapsed,alloccpus -j ${SLURM_JOBID}"
