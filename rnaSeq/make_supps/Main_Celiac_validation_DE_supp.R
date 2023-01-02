source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/lncRNA_funcs.R')

type = 'lncRNA'
type = 'ProteinCoding'

# set needed tables
metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v8.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)

TPM_file = sprintf('data/All_cohorts_%s_TPM_v8.txt',type)
TPM = read.table(file = TPM_file, header = T,sep = '\t', stringsAsFactors = F)
names(TPM) = make.names(names(TPM))

# ## Rectal UC
Name = 'Duodenum_Celiac'
names = c('SEEM_FC','Celiac_Leonard_FC')

cohort = 'SEEM'; cohort_name = 'SEEM'
diagnosis_prm = 'Dx_specific'
source = 'Duodenum'
to_test = c('Celiac','Control')
DE_name = sprintf('../%s/res_main3TPMf_TPMf_v8/%s_%s_%s_%s_vs_%s_deseq2_res_filtered.tsv', cohort, cohort, type, source, to_test[1],to_test[2])
DE_SEEM = read.table(file = DE_name, header = T,sep = '\t', stringsAsFactors = F)
genes_SEEM = DE_SEEM$Gene
df_SEEM = cohort_FC_pqval_avgTPM_table_by_geneList(genes_SEEM, DE_SEEM, TPM, metadata_df, cohort, cohort_name, source, diagnosis_prm, to_test )
names(DE_SEEM)[1:6] = sprintf('SEEM__%s',names(DE_SEEM))[1:6]

cohort = 'Celiac_Leonard'; cohort_name = 'Leonard'
source = 'Duodenum'
diagnosis_prm = 'Dx_specific'
to_test = c('Celiac_active','Control')
DE_name = sprintf('../%s/res_main3TPMf_TPMf_v8/%s_%s_%s_%s_vs_%s_deseq2_res.tsv', cohort, cohort, type, source, to_test[1],to_test[2])
DE_Leonard = read.table(file = DE_name, header = T,sep = '\t', stringsAsFactors = F)
df_Leonard = cohort_FC_pqval_avgTPM_table_by_geneList(genes_SEEM, DE_Leonard, TPM, metadata_df, cohort, cohort_name, source, diagnosis_prm, to_test )
# names(DE_Leonard) = sprintf('Leonard__%s',names(DE_Leonard))

df_SEEM = df_SEEM[,c('Gene','SEEM__Celiac_averageTPM','SEEM__Control_averageTPM')]
Celiac_df = merge(x = DE_SEEM, y=df_SEEM, by = 'Gene')
Celiac_df = merge(x = Celiac_df, y=df_Leonard, by = 'Gene')
# Celiac_df$Score = calc_score_v2(Celiac_df, c('SEEM__log2FoldChange','Leonard__logFC') )
Celiac_df$Score = calc_score_v3(Celiac_df, c('SEEM','Leonard') )

print( summary( as.factor( Celiac_df$Score ) ) )

out_file = sprintf('res/supps_v8/Celiac_main_validation_DE_%s_supp.txt',type)
write.table(x = Celiac_df, file = out_file, quote = F, sep = '\t',row.names = F)

# Name = 'Ileum_CD'
# names = c('SOURCE_FC','RISK_Ileum_iCD_FC')
# DE_files = c('../../SOURCE/res/SOURCE_lncRNA_terminal_ileum_CD_vs_Control_deseq2_res_filtered.tsv',
# '../../RISK/res/RISK_lncRNA_ileal_iCD_vs_Control_deseq2_res_filtered.tsv')

# # celiac
# Name = 'Duodenum_Celiac'
# names = c('SEEM_FC','Celiac_Leonard_FC')
# DE_files = c('../../SEEM/res/SEEM_lncRNA_Doudenum_Celiac_vs_Control_deseq2_res_filtered.tsv',
#              '../../Celiac_Leonard/res/Celiac_Leonard_lncRNA_Duodenum_Celiac_active_vs_Control_deseq2_res_filtered.tsv')
