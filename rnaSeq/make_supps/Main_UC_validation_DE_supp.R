source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/lncRNA_funcs.R')

# type = 'lncRNA'
type = 'ProteinCoding'

# set needed tables
metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v8.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)

TPM_file = sprintf('data/All_cohorts_%s_TPM_v8.txt',type)
TPM = read.table(file = TPM_file, header = T,sep = '\t', stringsAsFactors = F)

# ## Rectal UC
Name = 'Rectal_UC'
names = c('Protect_FC','RISK_Rectal_UC_FC')

cohort = 'Protect'; cohort_name = 'Protect'
source = 'Rectal'
diagnosis_prm = 'Dx'
to_test = c('UC','Control')
DE_name = sprintf('../%s/res_main3TPMf_TPMf_v8/%s_%s_%s_%s_vs_%s_deseq2_res_filtered.tsv', cohort, cohort, type, source, to_test[1],to_test[2])
DE_protect = read.table(file = DE_name, header = T,sep = '\t', stringsAsFactors = F)
genes_protect = DE_protect$Gene
df_protect = cohort_FC_pqval_avgTPM_table_by_geneList(genes_protect, DE_protect, TPM, metadata_df, cohort, cohort_name, source, diagnosis_prm, to_test )
names(DE_protect)[1:6] = sprintf('Protect__%s',names(DE_protect))[1:6]

cohort = 'RISK'; cohort_name = 'RISK_rectalUC'
source = 'Rectal'
diagnosis_prm = 'Dx_specific'
to_test = c('UC','Control')
DE_name = sprintf('../%s/res_main3TPMf_TPMf_v8/%s_%s_%s_%s_vs_%s_deseq2_res.tsv', cohort, cohort, type, source, to_test[1],to_test[2])
DE_RISK_rectalUC = read.table(file = DE_name, header = T,sep = '\t', stringsAsFactors = F)
df_RISK_UC = cohort_FC_pqval_avgTPM_table_by_geneList(genes_protect, DE_RISK_rectalUC, TPM, metadata_df, cohort, cohort_name, source, diagnosis_prm, to_test )
# names(DE_RISK_rectalUC) = sprintf('RISK_rectalUC__%s',names(DE_RISK_rectalUC))

df_protect = df_protect[,c('Gene','Protect__UC_averageTPM','Protect__Control_averageTPM')]
uc_df = merge(x = DE_protect, y=df_protect, by = 'Gene')
uc_df = merge(x = uc_df, y=df_RISK_UC, by = 'Gene')
# uc_df$Score = calc_score_v2(uc_df, c('Protect__log2FoldChange','RISK_rectalUC__logFC') )
uc_df$Score = calc_score_v3(uc_df, c('Protect','RISK_rectalUC') )

print( summary( as.factor( uc_df$Score ) ) )

out_file = sprintf('res/supps_v8/UC_main_validation_DE_%s_supp.txt',type)
dir.create('res/supps_v8/')
write.table(x = uc_df, file = out_file, quote = F, sep = '\t',row.names = F)

# Name = 'Ileum_CD'
# names = c('SOURCE_FC','RISK_Ileum_iCD_FC')
# DE_files = c('../../SOURCE/res/SOURCE_lncRNA_terminal_ileum_CD_vs_Control_deseq2_res_filtered.tsv',
# '../../RISK/res/RISK_lncRNA_ileal_iCD_vs_Control_deseq2_res_filtered.tsv')

# # celiac
# Name = 'Duodenum_Celiac'
# names = c('SEEM_FC','Celiac_Leonard_FC')
# DE_files = c('../../SEEM/res/SEEM_lncRNA_Doudenum_Celiac_vs_Control_deseq2_res_filtered.tsv',
#              '../../Celiac_Leonard/res/Celiac_Leonard_lncRNA_Duodenum_Celiac_active_vs_Control_deseq2_res_filtered.tsv')
