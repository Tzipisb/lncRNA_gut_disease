source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/lncRNA_funcs.R')

type = 'lncRNA'
type = 'ProteinCoding'

# set needed tables
metadata_file =  '../../metadata/lncRNA_meta_analysis_megamap_v4.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)

TPM_file = sprintf('../data/All_cohorts_%s_TPM.txt',type)
TPM = read.table(file = TPM_file, header = T,sep = '\t', stringsAsFactors = F)

# ## Rectal UC
Name = 'Ileum_CD'
names = c('SOURCE_FC','RISK_Ileum_iCD_FC')

cohort = 'SOURCE'; cohort_name = 'SOURCE'
source = 'terminal_ileum'
diagnosis_prm = 'Dx'
to_test = c('CD','Control')
DE_name = sprintf('../../%s/res_main3TPMf_TPMf/%s_%s_%s_%s_vs_%s_deseq2_res_filtered.tsv', cohort, cohort, type, source, to_test[1],to_test[2])
DE_SOURCE = read.table(file = DE_name, header = T,sep = '\t', stringsAsFactors = F)
genes_SOURCE = DE_SOURCE$Gene
df_SOURCE = cohort_FC_pqval_avgTPM_table_by_geneList(genes_SOURCE, DE_SOURCE, TPM, metadata_df, Cohort, cohort_name, source, diagnosis_prm, to_test )
names(DE_SOURCE)[1:6] = sprintf('SOURCE__%s',names(DE_SOURCE))[1:6]

cohort = 'RISK'; cohort_name = 'RISK_ilealCD'
source = 'ileal'
diagnosis_prm = 'Dx_specific'
to_test = c('iCD','Control')
DE_name = sprintf('../../%s/res_main3TPMf_TPMf/%s_%s_%s_%s_vs_%s_deseq2_res.tsv', cohort, cohort, type, source, to_test[1],to_test[2])
DE_RISK_ilealCD = read.table(file = DE_name, header = T,sep = '\t', stringsAsFactors = F)
df_RISK_ilealCD = cohort_FC_pqval_avgTPM_table_by_geneList(genes_SOURCE, DE_RISK_ilealCD, TPM, metadata_df, Cohort, cohort_name, source, diagnosis_prm, to_test )
# names(DE_RISK_rectalUC) = sprintf('RISK_rectalUC__%s',names(DE_RISK_rectalUC))

df_SOURCE = df_SOURCE[,c('Gene','SOURCE__CD_averageTPM','SOURCE__Control_averageTPM')]
cd_df = merge(x = DE_SOURCE, y=df_SOURCE, by = 'Gene')
cd_df = merge(x = cd_df, y=df_RISK_ilealCD, by = 'Gene')
# cd_df$Score = calc_score_v2(cd_df, c('SOURCE__log2FoldChange','RISK_ilealCD__logFC') )
cd_df$Score = calc_score_v3(cd_df, c('SOURCE','RISK_ilealCD') )

print( summary( as.factor( cd_df$Score ) ) )

out_file = sprintf('../res/supps/CD_main_validation_DE_%s_supp.txt',type)
write.table(x = cd_df, file = out_file, quote = F, sep = '\t',row.names = F)

# # celiac
# Name = 'Duodenum_Celiac'
# names = c('SEEM_FC','Celiac_Leonard_FC')
# DE_files = c('../../SEEM/res/SEEM_lncRNA_Doudenum_Celiac_vs_Control_deseq2_res_filtered.tsv',
#              '../../Celiac_Leonard/res/Celiac_Leonard_lncRNA_Duodenum_Celiac_active_vs_Control_deseq2_res_filtered.tsv')
