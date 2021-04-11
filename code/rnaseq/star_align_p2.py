#!/usr/bin/env python3

'''
__Task__
1. run star to generate second pass alignments
'''

import os
import glob
import subprocess
import argparse
import math


def run_star(in_path,out_path,batch,batch_size,array_max,live):
    output = 'star_align_p2_'+batch+'_%a.out'
    os.chdir(scripts_dir)
    open('star_align_p2.flag','w').close()  ## to serve as output for snakemake
    shell_args = ' '.join([in_path,out_path,str(batch_size),batch])
    cmd = ' '.join(['sbatch','--output='+output,'--array=1-'+str(array_max),'star_align_p2.sh',shell_args])
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
str_dir = work_dir + 'star/'
suff = '.fastq.gz'
##
## parse command line arguments
##
parser = argparse.ArgumentParser(description='star align 2')
parser.add_argument('--batch_size',help='number of samples to be processed in single thread',action='store',dest='batch_size')
parser.add_argument('--live',help='execute a live run',action='store_true',dest='run_type')
parser.add_argument('--test',help='execute a test run',action='store_false',dest='run_type')
parser.set_defaults(run_type=False,batch_size=10)
args = parser.parse_args()
##
## prepare and execute
##
batches = ['chrom9', 'chrom20']
for batch in batches:
    in_path = fq_dir+batch+'/'
    out_path = str_dir+batch+'/p2/'
    if os.path.exists(str_dir+batch):
        pass
    else:
        os.mkdir(str_dir+batch)
    if os.path.exists(out_path):
        pass
    else:
        os.mkdir(out_path)
    fileset = glob.glob(in_path+'*'+suff)
    array_max = math.ceil(len(fileset)/int(args.batch_size))
    run_star(in_path,out_path,batch,int(args.batch_size),array_max,args.run_type)
