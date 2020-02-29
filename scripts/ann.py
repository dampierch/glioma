#!/usr/bin/env python3

'''
__Task__
1. generate preliminary annotations for analysis
'''

import os
import glob
import argparse
import pandas as pd


def make_ann():
    ## full annotations from mourad
    filename = 'glioma_ann.xlsx'
    ann = pd.read_excel(ann_dir+filename)
    in_fields = ['Sample Name', 'RIN values', 'Sample_pheno', 'chrom', 'day', 'ex_date']
    ann = ann[in_fields]
    fix_field_names = ['sample_id','rin','sample_pheno','chrom','day','ex_date']
    ann.columns = fix_field_names
    ann['cell_line'] = ann.sample_pheno.str.strip().str.split().str[1]
    ann['treatment_group'] = ann.sample_pheno.str.strip().str.split().str[2].str.replace('\d+','')
    ann['pair_id'] = ann.cell_line + '_' + ann.day.astype(str)
    samples_per_exp = 14
    ann['num'] = ann.sample_id.str.slice(-2).astype(int)%samples_per_exp
    ann.loc[ann['num'] == 0, 'num'] = samples_per_exp
    ann['fq_id'] = ann.sample_id.str.replace('l','1') + '_' + 'S' + ann.num.astype(str)
    return(ann)


##
## set stable environmental variables
##
home = os.environ['HOME']
scripts_dir = home + '/projects/glioma/scripts/'
work_dir = '/scratch/chd5n/glioma/'
fq_dir = work_dir + 'fastq/'
fqc_dir = work_dir + 'fqc/'
mqc_dir = work_dir + 'mqc/'
qnt_dir = work_dir + 'quant/'
str_dir = work_dir + 'star/'
rsm_dir = work_dir + 'rsem/'
ann_dir = work_dir + 'annotations/'
suff = '.fastq.gz'
##
## parse command line arguments
##
parser = argparse.ArgumentParser(description='annotations')
parser.add_argument('--live',help='execute a live run',action='store_true',dest='run_type')
parser.add_argument('--test',help='execute a test run',action='store_false',dest='run_type')
parser.set_defaults(run_type=False)
args = parser.parse_args()
open('annotations.flag','w').close()  ## to serve as output for snakemake
##
## prepare and merge on left
##
ann = make_ann()
# qcann = get_qc_metrics_1()
# ann_1 = qcann.merge(ann_629,how='left',on='vm_id')
# ann_1 = ann_1.merge(spann,how='left',on='vm_id')
##
## write annotation
##
if args.run_type == True:
    filepath = ann_dir
    filename = 'ann.tsv'
    ann.to_csv(filepath+filename,sep='\t',header=True,index=False)
else:
    print('test complete')
