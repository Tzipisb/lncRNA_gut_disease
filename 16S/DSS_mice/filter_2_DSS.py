#!/usr/bin/env python
import os
import sys
sys.path.insert(0, '/pita/users/tzipi/code/qiime/')
from qiime2_funcs import *

## setting veriables
verbose = True
threads = 20
rar_num = 2000
decompres_flag = True

# map_file = 'data/DB-19_DSS_mice_map.txt'
# out_path = 'res_all/'

# map_file = 'data/DB-19_DSS_mice_DH_map.txt'
# out_path = 'res_DH/'

map_file = 'data/DB-19_DSS_mice_map_noOut2.txt'
out_path = 'res_all_noOut2/'


os.system('mkdir ' + out_path)


## setting input biom files (table, rep-seqs)
in_path = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/'
table_file = in_path + 'table.qza'
rep_seqs_file = in_path + 'rep-seqs.qza'

os.system('cp ' + rep_seqs_file + ' ' + out_path + 'rep-seqs.qza')
os.system('cp ' + table_file + ' ' + out_path + '/table_all_samps.qza')

## filtering general biom to wanted samples
os.system('qiime feature-table filter-samples \
  --i-table ' + out_path + '/table_all_samps.qza \
  --m-metadata-file ' + map_file + ' \
  --o-filtered-table ' + out_path + '/table.qza')

## run downstream analysis
os.system('cp ' + map_file + ' ' + out_path + '/map.txt')
# phylogenetic_analysis(out_path, verbose, threads = threads)
get_final_reads_count(out_path)
update_map_with_reads(out_path)
phylogenetic_analysis_sepp(out_path, verbose, threads = threads)
diversity_analysis(out_path, map_file, rar_num, verbose, decompres_flag, threads = threads)
taxa_analysis(out_path, map_file, verbose, threads = threads)
# clean_big_files(out_path)


open_files_for_R(out_path)

# for picrust
os.system('qiime tools export \
	--input-path ' + out_path + '/rep-seqs.qza \
	--output-path rep-seqs' )


# for beta diversity Genotype
os.system('qiime diversity beta-group-significance \
	--i-distance-matrix ' + out_path + '/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
	--m-metadata-file ' + map_file + ' \
	--m-metadata-column Genotype \
	--p-pairwise \
	--o-visualization ' + out_path + '/core-metrics-results/unweighted-unifrac-Genotype-significance')

os.system('qiime tools export \
	--input-path ' + out_path + '/core-metrics-results/unweighted-unifrac-Genotype-significance.qzv \
	--output-path ' + out_path + '/core-metrics-results/unweighted-unifrac-Genotype-significance')

os.system('biom convert -i ' + out_path + '/biom/feature-table.biom -o ' + out_path + '/biom/feature-table.txt --to-tsv')

