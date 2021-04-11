#!/bin/bash

## set variables
script_name=setup.sh
project_name=glioma
scratch=/scratch/chd5n/
projects=~/projects/
scratch_subdirs="annotations fastq rdata results fqc mqc quant star rsem"
results_subdirs="plots tables"
proj_subdirs="logs scripts"

## make scratch directory
if [[ -e ${scratch}${project_name} ]]; then
  :
else
  mkdir ${scratch}${project_name}
fi

for i in ${scratch_subdirs}; do
  if [[ -e ${scratch}${project_name}/${i} ]]; then
    :
  else
    mkdir ${scratch}${project_name}/${i}
  fi
done

for i in ${results_subdirs}; do
  if [[ -e ${scratch}${project_name}/results/${i} ]]; then
    :
  else
    mkdir ${scratch}${project_name}/results/${i}
  fi
done

## make project directory
if [[ -e ${projects}${project_name} ]]; then
  :
else
  mkdir ${projects}${project_name}
fi

for i in ${proj_subdirs}; do
  if [[ -e ${projects}${project_name}/${i} ]]; then
    :
  else
    mkdir ${projects}${project_name}/${i}
  fi
done

## download fastq
touch ${script_name%.sh}.flag  ## to serve as output for snakemake
echo starting download
origin=https://bccl.csbc.vcu.edu/external/Graham/
dest=${scratch}${project_name}/fastq
wget -P ${dest} --recursive --level=2 --no-parent --no-host-directories ${origin}
mv ${dest}/MiSeq_191217_1x75/ ${dest}/chrom9/
mv ${dest}/MiSeq_191218_1x75/ ${dest}/chrom20/

# FINISHED --2020-02-24 19:11:52--
# Total wall clock time: 2h 6m 52s
# Downloaded: 47 files, 32G in 2h 6m 50s (4.31 MB/s)
