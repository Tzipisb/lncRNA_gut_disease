source('/pita/users/tzipi/code/R_figs/machine_learing_funcs.R')
library(ggplot2)
library(AUC)
library(randomForest)

set.seed(4)

fix_group_size_bias2 = function(data, check_val)
{
  groups = levels( check_val )
  min_group_size = min( summary( check_val ) )
  
  pos = c()
  for (g in groups)
  {
    temp =  which( check_val == g )
    pos = c( pos, temp[1:min_group_size] )
  }
  return(pos)
}

basic_rf_run = function(metadata_df_sf, check_val, tpm_sf)
{  
  output.forest <- randomForest(as.factor(metadata_df_sf[[check_val]]) ~ .,data=tpm_sf)
  print(output.forest)
  
  imp = as.data.frame(importance(output.forest,type = 2))
  imp = data.frame( taxa = rownames(imp), MeanDecreaseGini = imp$MeanDecreaseGini)
  imp = imp[ rev(order(imp$MeanDecreaseGini)) ,]
  
  res = roc(output.forest$votes[,1],factor(1 * (output.forest$y==levels(output.forest$y)[1] )))
  auc_res = auc(res,min = 0, max = 1) 
  
  score_and_lab = data.frame(pred = output.forest$votes[,1], survived = output.forest$y)
  youden = calc_youden(score_and_lab)
  
  roc_data = data.frame(fpr = res$fpr, tpr = res$tpr)
  
  return( list(output.forest, imp, auc_res, score_and_lab, youden, roc_data) )
}

set_data_for_rf = function(tpm, metadata_df, cohort, test_source, to_test, diagnosis_prm, check_val)
{
  tpm = tpm[, make.names( metadata_df$SampleID )]
  
  #             #### choose data to look at (usign DE genes)
  # # filter to wanted cohort
  pos =  metadata_df$Cohort == cohort
  tpm = tpm[,pos]
  metadata_df = metadata_df[pos,]
  
  # extra filteres
  pos = metadata_df$Source %in% test_source
  pos2 = metadata_df[[diagnosis_prm]] %in% to_test
  pos = pos & pos2 
  # pos = pos2
  tpm = tpm[,pos]
  metadata_df = metadata_df[pos,]
  
  
  # tpm2 = tpm[rowSums(is.na(tpm)) != ncol(tpm),] # remove NAs lines
  tpm2 = tpm[rowSums(is.na(tpm)) == 0,] # remove any line with >1 NA
  tpm2 = tpm2[rowSums(tpm2) != 0,] # remove any gene with 0 expression
  # input_data = t(tpm2)
  tpm2 = as.data.frame(t(tpm2))
  names(tpm2) = make.names(names(tpm2))
  
  # pos = fix_group_size_bias2(tpm2, as.factor(metadata_df$Dx))
  pos = !is.na(metadata_df[[check_val]])
  metadata_df_sf = metadata_df[pos,]
  tpm_sf = tpm2[pos,]
  
  return(list(metadata_df_sf, tpm_sf))
}

out_path = 'plots/random_forest_outcome/'
dir.create(out_path)

type = 'lncRNA' 
# type = 'ProteinCoding'

cohort = 'Protect'
source = 'Rectal'
test_source = source
diagnosis_prm = 'Dx'
to_test = c('UC','Control')

name = 'PROTECT'

# check_val = 'protect_W4R'
# check_val = 'protect_3years_colectomy'
check_val = 'protect_CSFREE_REMISSION_WK52_KM'

# genes_file = sprintf( 'plots/venns/lists/%s_%s_TPM_filter.txt', cohort,  type )
# g_df = read.table(file = genes_file, header = F,sep = '\t')
# genes = g_df$V1
genes_file = sprintf('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/Protect/res_PUCAI/Protect_%s_protect_BASELINE_Pucai_Category_1_vs_3_deseq2_res_filtered.tsv', type)
g_df = read.table(file = genes_file, header = T,sep = '\t')
genes = g_df$Gene
name = sprintf('%s_severityGenes', name)

# genes_file = '/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/res/WGCNA/Protect_gene_module_table_lncRNA.txt'
# g_df = read.table(file = genes_file, header = T,sep = '\t')
# genes = g_df$Gene[g_df$module_color == 'red']

tpm_file =  sprintf('data/All_cohorts_%s_TPM_v7.txt',type)
tpm = read.table(file = tpm_file, header = T, row.names = 1)
tpm = tpm[as.character(genes),]

metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v7.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na','_'))
metadata_df = metadata_df[order(metadata_df$SampleID),]

metadata_df = metadata_df[metadata_df$protect_BASELINE_Pucai_Category %in% c(2,3) & !is.na(metadata_df$protect_BASELINE_Pucai_Category),]
name = sprintf('%s_noMild',name)

res = set_data_for_rf(tpm, metadata_df, cohort, test_source, to_test, diagnosis_prm, check_val)
metadata_df_sf = res[[1]]
tpm_sf = res[[2]]

res = basic_rf_run(metadata_df_sf, check_val, tpm_sf)
output.forest_lnc = res[[1]]
imp_lnc = res[[2]]
auc_res_lnc = res[[3]]
score_and_lab_lnc = res[[4]]
youden_lnc = res[[5]]
roc_data_lnc = res[[6]]
  


type = 'ProteinCoding'

genes_file = sprintf('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/Protect/res_PUCAI/Protect_%s_protect_BASELINE_Pucai_Category_1_vs_3_deseq2_res_filtered.tsv', type)
g_df = read.table(file = genes_file, header = T,sep = '\t')
genes = g_df$Gene
name = sprintf('%s_severityGenes', name)

# genes_file = '/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/res/WGCNA/Protect_gene_module_table_lncRNA.txt'
# g_df = read.table(file = genes_file, header = T,sep = '\t')
# genes = g_df$Gene[g_df$module_color == 'red']

tpm_file =  sprintf('data/All_cohorts_%s_TPM_v7.txt',type)
tpm = read.table(file = tpm_file, header = T, row.names = 1)
tpm = tpm[as.character(genes),]

metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v7.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na','_'))
metadata_df = metadata_df[order(metadata_df$SampleID),]

metadata_df = metadata_df[metadata_df$protect_BASELINE_Pucai_Category %in% c(2,3) & !is.na(metadata_df$protect_BASELINE_Pucai_Category),]
name = sprintf('%s_noMild',name)

res = set_data_for_rf(tpm, metadata_df, cohort, test_source, to_test, diagnosis_prm, check_val)
metadata_df_sf = res[[1]]
tpm_sf = res[[2]]

res = basic_rf_run(metadata_df_sf, check_val, tpm_sf)
output.forest_pc = res[[1]]
imp_pc = res[[2]]
auc_res_pc = res[[3]]
score_and_lab_pc = res[[4]]
youden_pc = res[[5]]
roc_data_pc = res[[6]]


roc_data_pc$Type = 'ProteinCoding'
roc_data_lnc$Type = 'lncRNA'
roc_data = rbind( roc_data_lnc, roc_data_pc )

p_roc = ggplot(roc_data) + 
  geom_line(aes(fpr,tpr, group = Type, color = Type), size=1.5) + xlab("FPR") + ylab("TPR") + theme_bw() + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="darkgrey", linetype="dashed") + ggtitle(sprintf('%s ROC', check_val)) +
  scale_colour_brewer(palette = 'Set2', name = '') + 
  theme(panel.grid = element_blank(), legend.position =c(0.73,0.22)) 
# ggsave(sprintf('%s/RF_%s_%s_lnc_auc%f_pc_auc%f.tiff', out_path, cohort, check_val, auc_res_lnc, auc_res_pc),plot = p_roc, device = 'tiff', width=3,height =3, compression  = 'lzw')


type = 'lncPrtCd'

genes_file = sprintf('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/Protect/res_PUCAI/Protect_%s_protect_BASELINE_Pucai_Category_1_vs_3_deseq2_res_filtered.tsv', 'lncRNA')
g_df = read.table(file = genes_file, header = T,sep = '\t')
genes_lnc = g_df$Gene

genes_file = sprintf('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/Protect/res_PUCAI/Protect_%s_protect_BASELINE_Pucai_Category_1_vs_3_deseq2_res_filtered.tsv', 'ProteinCoding')
g_df = read.table(file = genes_file, header = T,sep = '\t')
genes_pc = g_df$Gene

genes = c(genes_lnc, genes_pc)

tpm_file =  sprintf('data/All_cohorts_%s_TPM_v7.txt',type)
tpm = read.table(file = tpm_file, header = T, row.names = 1)
tpm = tpm[as.character(genes),]

metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v7.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na','_'))
metadata_df = metadata_df[order(metadata_df$SampleID),]

metadata_df = metadata_df[metadata_df$protect_BASELINE_Pucai_Category %in% c(2,3) & !is.na(metadata_df$protect_BASELINE_Pucai_Category),]
name = sprintf('%s_noMild',name)

res = set_data_for_rf(tpm, metadata_df, cohort, test_source, to_test, diagnosis_prm, check_val)
metadata_df_sf = res[[1]]
tpm_sf = res[[2]]

res = basic_rf_run(metadata_df_sf, check_val, tpm_sf)
output.forest_pc = res[[1]]
imp_lpc = res[[2]]
auc_res_lpc = res[[3]]
score_and_lab_lpc = res[[4]]
youden_lpc = res[[5]]
roc_data_lpc = res[[6]]

roc_data_lpc$Type = 'lncRNA +\nProteinCoding'
roc_data_all = rbind( roc_data, roc_data_lpc )

p_roc2 = ggplot(roc_data_all) + 
  geom_line(aes(fpr,tpr, group = Type, color = Type), size=1.5) + xlab("FPR") + ylab("TPR") + theme_bw() + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="darkgrey", linetype="dashed") + ggtitle(sprintf('%s ROC', check_val)) +
  scale_colour_brewer(palette = 'Set2', name = '') + 
  theme(panel.grid = element_blank(), legend.position =c(0.73,0.22), legend.background = element_blank()) 
ggsave(sprintf('%s/RF_%s_%s_lnc_auc%f_pc_auc%f_both_auc%f.tiff', out_path, cohort, check_val, auc_res_lnc, auc_res_pc,auc_res_lpc),plot = p_roc2, device = 'tiff', width=3,height =3, compression  = 'lzw')

p_roc_lpc = ggplot(roc_data_lpc) + 
  geom_line(aes(fpr,tpr, group = Type, color = Type), size=1.5) + xlab("FPR") + ylab("TPR") + theme_bw() + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="darkgrey", linetype="dashed") + ggtitle(sprintf('%s ROC', check_val)) +
  scale_colour_brewer(palette = 'Set2', name = '') + 
  theme(panel.grid = element_blank(), legend.position =c(0.70,0.22)) 
ggsave(sprintf('%s/RF_%s_%s_both_auc%f.tiff', out_path, cohort, check_val,auc_res_lpc),plot = p_roc_lpc, device = 'tiff', width=3,height =3, compression  = 'lzw')

