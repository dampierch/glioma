#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10000
#SBATCH --time=5:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chd5n@virginia.edu

## requires ~2.5h, 230MB for 434 merged bc fq.gz files with 9 arrays of 50

script_name=fqc.sh

## header
pwd; hostname; date

## modifiable parameters
in_path=$1  ## take from positional arguments on command line
out_path=$2
batch_size=$3
fe=.fastq.gz

## stable/auto parameters
printf "start ${script_name}\nbatch size: ${batch_size}\ninput: ${in_path}\noutput: ${out_path}\n"
module load trimgalore
SATID=${SLURM_ARRAY_TASK_ID}
d=$(( ${batch_size} * ( ${SATID} - 1 ) ))
q=$(( ${batch_size} * ${SATID} ))
if [ ${SATID} -eq 1 ]; then
  find ${in_path} -maxdepth 1 -type f -name "*${fe}" | sort | uniq | sed -e "${q}q" > ${in_path}fqcIdx${SATID}.temp
elif [ ${SATID} -gt 1 ]; then
  find ${in_path} -maxdepth 1 -type f -name "*${fe}" | sort | uniq | sed -e "1,${d}d;${q}q" > ${in_path}fqcIdx${SATID}.temp
fi

## execute
for fn in $(<${in_path}fqcIdx${SATID}.temp); do
  fastqc --outdir=${out_path} ${fn}
done
rm ${in_path}fqcIdx${SATID}.temp

## footer
printf "end ${script_name}\n"; date
echo "seff ${SLURM_JOBID}"
echo "sacct -o reqmem,maxrss,elapsed,alloccpus -j ${SLURM_JOBID}"
