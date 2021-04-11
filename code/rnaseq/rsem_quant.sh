#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10000      ## mem-per-node
#SBATCH --time=10:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=END
#SBATCH --mail-user=chd5n@virginia.edu

  ## this script assumes setup_rsem.sh script has already been run from caseylab-home
  ## took 20h, 4GB for 5 samples prior to Thanksgiving with --ci-calc (1GB default)
  ## took max 4h 30m, 3-4GB for 5 samples without --ci-calc
  ## http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html#arguments estimates about 2GB needed for non-ci-calc activities

script_name=rsem_quant.sh

## header
pwd; hostname; date

## modifiable variables
in_path=$1  ## from python, which in turn is from Snakefile
out_path=$2
batch_size=$3
index_path=/scratch/chd5n/reference-genome/rsem_index/
file_suffix=Aligned.toTranscriptome.out.bam

## auto variables
cpus_per_task=${SLURM_CPUS_PER_TASK}
ci_mem=$(( ${SLURM_MEM_PER_NODE} - 4000 ))
SATID=${SLURM_ARRAY_TASK_ID}
d=$(( ${batch_size} * ( ${SATID} - 1 ) ))
q=$(( ${batch_size} * ${SATID} ))
if [ ${SATID} -eq 1 ]; then
  find ${in_path} -mindepth 2 -maxdepth 2 -type f -name "${file_suffix}" | sort | uniq | sed -e "${q}q" > ${in_path}rsemIdx${SATID}.tmp
elif [ ${SATID} -gt 1 ]; then
  find ${in_path} -mindepth 2 -maxdepth 2 -type f -name "${file_suffix}" | sort | uniq | sed -e "1,${d}d;${q}q" > ${in_path}rsemIdx${SATID}.tmp
fi

## execute
printf "start ${script_name}\nbatch size: ${batch_size}\nindex: ${index_path}\ninput: ${in_path}\noutput: ${out_path}\n"
module load gcc/7.1.0 rsem/1.3.0

for fn in $(<${in_path}rsemIdx${SATID}.tmp); do
  sn=`echo ${fn} | awk -v FS='/' '{print $(NF-1)}'`
  echo working on ${sn}
  mkdir ${out_path}${sn}

  rsem-calculate-expression \
    --num-threads ${cpus_per_task} \
    --alignments \
    --fragment-length-mean 75 \
    --fragment-length-sd 0 \
    --estimate-rspd \
    --no-bam-output \
    ${fn} \
    ${index_path}GRCh38 \
    ${out_path}${sn}/${sn}.rsem

    # --calc-ci \
    # --ci-memory ${ci_mem} \

  echo done with ${sn}
done
rm ${in_path}rsemIdx${SATID}.tmp

## footer
printf "end ${script_name}\n"; date
echo "seff ${SLURM_JOBID}"
echo "sacct -o reqmem,maxrss,elapsed,alloccpus,nodelist -j ${SLURM_JOBID}"
