#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10000
#SBATCH --time=1:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=END
#SBATCH --mail-user=chd5n@virginia.edu

## requires ~10min, 250MB for bc data

script_name=mqc.sh
in_path=$1
out_path=$2

## header
pwd; hostname; date
printf "start ${script_name}\ninput: ${in_path}\noutput: ${out_path}\n\n"

##########
# set up #
##########

module purge
module load anaconda/5.2.0-py3.6
source activate multiqc_env

###########
# execute #
###########

## use --force to overwright existing summaries, --ignore to skip directory
multiqc ${in_path} --outdir ${out_path}

## footer
printf "end ${script_name}\n"; date
echo "seff ${SLURM_JOBID}"
echo "sacct -o reqmem,maxrss,elapsed,alloccpus -j ${SLURM_JOBID}"
