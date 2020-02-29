#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=15000
#SBATCH --time=6:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=END
#SBATCH --mail-user=chd5n@virginia.edu

  ## this script assumes setup_salmon.sh script has already been run for
  ## appropriate kmer length and stored in /scratch/chd5n/Reference_genome/salmon_index
  ## recent versions of salmon now perform selective alignment as opposed to the
  ## original quasi-mapping
  ## takes about ? GB, ? min for human

script_name=salgc.sh

## header
pwd; hostname; date

## variables
in_path=$1  ## from python, which in turn is from Snakefile
out_path=$2
batch_size=$3
index_path=/scratch/chd5n/Reference_genome/salmon_index/gentrome_k31/
file_suffix=fastq.gz
cpus_per_task=${SLURM_CPUS_PER_TASK}

## stable/auto parameters
printf "start ${script_name}\nbatch size: ${batch_size}\nindex: ${index_path}\ninput: ${in_path}\noutput: ${out_path}\n"
module load gcc/7.1.0 salmon/1.0.0
SATID=${SLURM_ARRAY_TASK_ID}
d=$(( ${batch_size} * ( ${SATID} - 1 ) ))
q=$(( ${batch_size} * ${SATID} ))
if [ ${SATID} -eq 1 ]; then
  find ${in_path} -maxdepth 1 -type f -name "*${file_suffix}" | awk -v FS="_R" '{print $1}' | sort | uniq | sed -e "${q}q" > ${in_path}salGcIdx${SATID}.tmp
elif [ ${SATID} -gt 1 ]; then
  find ${in_path} -maxdepth 1 -type f -name "*${file_suffix}" | awk -v FS="_R" '{print $1}' | sort | uniq | sed -e "1,${d}d;${q}q" > ${in_path}salGcIdx${SATID}.tmp
fi

## execute
for fn in $(<${in_path}salGcIdx${SATID}.tmp); do  ## fn includes path to input
  sn=`echo ${fn} | awk -v FS='/' '{print $NF}'`
  ## all samples are single-end in glioma dataset
  ## length is 75bp; salmon index is set to k31
  salmon quant \
    -p ${cpus_per_task} \
    -i ${index_path} \
    -l A \
    -r <(gunzip -c ${fn}_R1_001.${file_suffix}) \
    -o ${out_path}${sn} \
    --gcBias \
    --validateMappings
done
rm ${in_path}salGcIdx${SATID}.tmp

## footer
printf "end ${script_name}\n"; date
echo "seff ${SLURM_JOBID}"
echo "sacct -o reqmem,maxrss,elapsed,alloccpus -j ${SLURM_JOBID}"
