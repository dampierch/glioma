#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50000
#SBATCH --time=2:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH	--mail-type=END
#SBATCH	--mail-user=chd5n@virginia.edu

  ## this script assumes ENSEMBL reference genome and GTF are already
  ## downloaded to /scratch/chd5n/reference-genome/ and splice junctions are
  ## stored in project directory SJ.out.tab files from first pass; takes
  ## those as input, generates project-specific genome index for star second
  ## pass alignment
  ## takes about 35 GB, 35 min for human

  ## usage: sbatch --output=star_idx_p2_chrom20.out star_idx_p2.sh chrom20

  ## process took ~ 50min, 35GB for GRCh38 filt $7>2

script_name=star_idx_p2.sh

## header
pwd; hostname; date

## variables
work_dir=/scratch/chd5n/glioma/star/
batch=$1
star_idx_dir=star_index/
ref_home=/scratch/chd5n/reference-genome/
ref_fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa
ref_gtf=Homo_sapiens.GRCh38.98.gtf
read_length=75
overhang=$(( ${read_length}-1 ))
cpus_per_task=${SLURM_CPUS_PER_TASK}
align_dir=${work_dir}${batch}/
sjouttabs=`find ${align_dir}p1/ -mindepth 2 -maxdepth 2 -type f -name "SJ.out.tab"`
sjoutfilt=SJ_p1_filtered.tab

## prepare star index directory structure
if [[ -e ${work_dir}${batch}/${star_idx_dir} ]]; then
  :
else
  mkdir ${work_dir}${batch}/${star_idx_dir}
fi

## execute
printf "start ${script_name}\ninput: ${ref_home}${ref_fasta}\noutput: ${work_dir}${batch}/${star_idx_dir}\n"
module purge
module load star/2.7.2b

cat ${sjouttabs} | awk '($5>0 && $7>2 && $6==0)' | cut -f1-6 | sort | uniq > ${align_dir}${sjoutfilt}  ## >2 uniquely mapped reads per sample to pass filter

STAR \
  --runMode genomeGenerate \
  --runThreadN ${cpus_per_task} \
  --genomeDir ${work_dir}${batch}/${star_idx_dir} \
  --genomeFastaFiles ${ref_home}${ref_fasta} \
  --sjdbOverhang ${overhang} \
  --sjdbGTFfile ${ref_home}${ref_gtf} \
  --sjdbFileChrStartEnd ${align_dir}${sjoutfilt} \
  --outTmpDir ${work_dir}${batch}/_STARtmp

## footer
printf "end ${script_name}\n"; date
echo "seff ${SLURM_JOBID}"
echo "sacct -o reqmem,maxrss,elapsed,alloccpus,nodelist -j ${SLURM_JOBID}"
