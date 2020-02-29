#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40000
#SBATCH --time=2:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=END
#SBATCH --mail-user=chd5n@virginia.edu

  ## process took 35m, 35GB for 28 compressed single-end samples in
  ## 6 batches of 5 on 10 cores

script_name=star_align_p2.sh

## header
pwd; hostname; date

## modifiable parameters
in_path=$1  ## take from positional arguments on command line
out_path=$2
batch_size=$3
batch=$4
index_dir=/scratch/chd5n/glioma/star/${batch}/star_index/  ## project-specific index for second pass
file_suffix=fastq.gz
read_length=75

## stable/auto parameters
overhang=$(( ${read_length}-1 ))
cpus_per_task=${SLURM_CPUS_PER_TASK}
SATID=${SLURM_ARRAY_TASK_ID}
d=$(( ${batch_size} * ( ${SATID} - 1 ) ))
q=$(( ${batch_size} * ${SATID} ))
if [ ${SATID} -eq 1 ]; then
  find ${in_path} -maxdepth 1 -type f -name "*${file_suffix}" | awk -v FS="_R" '{print $1}' | sort | uniq | sed -e "${q}q" > ${in_path}strAlignIdx${SATID}.tmp
elif [ ${SATID} -gt 1 ]; then
  find ${in_path} -maxdepth 1 -type f -name "*${file_suffix}" | awk -v FS="_R" '{print $1}' | sort | uniq | sed -e "1,${d}d;${q}q" > ${in_path}strAlignIdx${SATID}.tmp
fi

## execute
printf "start ${script_name}\nbatch size: ${batch_size}\nindex: ${index_dir}\ninput: ${in_path}\noutput: ${out_path}\n"
module purge
module load star/2.7.2b
for fn in $(<${in_path}strAlignIdx${SATID}.tmp); do
  sn=`echo ${fn} | awk -v FS='/' '{print $NF}'`
  echo working on ${sn}
  mkdir ${out_path}${sn}
  ## all samples are single-end in glioma dataset
  ## length is 75; index built with 74bp sjdbOverhang
  STAR \
    --genomeDir ${index_dir} \
    --genomeLoad NoSharedMemory \
    --readFilesIn ${fn}_R1_001.${file_suffix} \
    --readFilesCommand zcat \
    --runThreadN ${cpus_per_task} \
    --outFileNamePrefix ${out_path}${sn}/ \
    --sjdbOverhang ${overhang} \
    --limitSjdbInsertNsj 1200000 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJDBoverhangMin 1 \
    --alignSJoverhangMin 8 \
    --chimSegmentMin 15 \
    --chimJunctionOverhangMin 15 \
    --chimMainSegmentMultNmax 1 \
    --chimOutType WithinBAM SoftClip \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.1 \
    --outFilterType BySJout \
    --outFilterScoreMinOverLread 0.33 \
    --outFilterMatchNminOverLread 0.33 \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM Unsorted \
    --outSAMmode Full \
    --outSAMunmapped Within \
    --outSAMattrRGline ID:rg1 SM:sm1 \
    --outSAMattributes NH HI AS nM NM ch \
    --quantMode TranscriptomeSAM GeneCounts
  echo done with ${sn}
done
rm ${in_path}strAlignIdx${SATID}.tmp

## footer
printf "end ${script_name}\n"; date
echo "seff ${SLURM_JOBID}"
echo "sacct -o reqmem,maxrss,elapsed,alloccpus,nodelist -j ${SLURM_JOBID}"
