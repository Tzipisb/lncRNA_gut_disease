#!/usr/bin/env python
import os
import pandas as pd 
import sys
sys.path.insert(0, '/pita/users/tzipi/code/rnaSeq')
from run_kallisto_funcs import *

kallisto_index = '/pita/pub/bin/kallisto_linux-v0.42.5/index/hg_gencode_v24'

in_path = 'RISK/fastq_files/'
out_path = 'RISK/'

# running kallisto
# kal_out_path = out_path + 'kallisto_s31/'
kal_out_path = out_path + 'kallisto/'
os.system('mkdir ' + kal_out_path)
run_kallisto_paired_kal_42_5_gencode_v24(in_path, kal_out_path, r1_suffix = '_1.fastq', r2_suffix = '_2.fastq', kallisto_index = kallisto_index)


