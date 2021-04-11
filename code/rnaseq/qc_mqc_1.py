#!/usr/bin/env python3


'''
__Tasks__

1. aggregate fastqc output across all samples in each batch with multiqc
'''


import os
import glob
import subprocess
import argparse
import math


def run_multiqc(in_path,out_path,batch,live):
    output = 'mqc_'+batch+'.out'
    os.chdir(scripts_dir)
    open('qc_mqc_1.flag','w').close()  ## to serve as output for snakemake
    shell_args = ' '.join([in_path,out_path])
    cmd = ' '.join(['sbatch','--output='+output,'mqc.sh',shell_args])
    if live == True:
        subprocess.call(cmd, shell=True)
    else:
        print(cmd)


##
## set stable environmental variables
##
home = os.environ['HOME']
scripts_dir = home + '/projects/glioma/code/rnaseq/'
work_dir = '/scratch/chd5n/glioma/'
fq_dir = work_dir + 'fastq/'
fqc_dir = work_dir + 'fqc/'
mqc_dir = work_dir + 'mqc/'
qnt_dir = work_dir + 'quant/'
suff = '.fastq.gz'
##
## parse command line arguments
##
parser = argparse.ArgumentParser(description='multiqc 1')
parser.add_argument('--live',help='execute a live run',action='store_true',dest='run_type')
parser.add_argument('--test',help='execute a test run',action='store_false',dest='run_type')
parser.set_defaults(run_type=False)
args = parser.parse_args()
##
## prepare and execute mqc 1
##
batches = ['chrom9', 'chrom20']
for batch in batches:
    in_path = fqc_dir+batch+'/'
    out_path = mqc_dir+batch+'/'
    if os.path.exists(out_path):
        pass
    else:
        os.mkdir(out_path)
    run_multiqc(in_path,out_path,batch,args.run_type)
