source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/lncRNA_funcs.R')

type = 'lncRNA'
# type = 'ProteinCoding'

# set needed tables
metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v4.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)

TPM_file = sprintf('data/All_cohorts_%s_TPM.txt',type)
TPM = read.table(file = TPM_file, header = T,sep = '\t', stringsAsFactors = F)

## get gene list to look at - validated from 3 main
UC_file = sprintf('res/supps/UC_main_validation_DE_%s_supp.txt',type)
CD_file = sprintf('res/supps/CD_main_validation_DE_%s_supp.txt',type)
Celiac_file = sprintf('res/supps/Celiac_main_validation_DE_%s_supp.txt',type)
UC = read.table(UC_file, header = T, stringsAsFactors = F)
CD = read.table(CD_file, header = T, stringsAsFactors = F)
Celiac = read.table(Celiac_file, header = T, stringsAsFactors = F)
genes = unique( c( UC$Gene[UC$Score == 'fit'], 
                   CD$Gene[CD$Score == 'fit'],
                   Celiac$Gene[Celiac$Score == 'fit'] ) )


# ## Rectal UC
Name = 'Duodenum_Celiac'
names = c('SEEM_FC','Celiac_Leonard_FC')

cohort = 'SEEM'; cohort_name = 'SEEM'
diagnosis_prm = 'Dx_specific'
source = 'Doudenum'
to_test = c('Celiac','Control')
DE_name = sprintf('../../%s/res_main3TPMf_TPMf/%s_%s_%s_%s_vs_%s_deseq2_res.tsv', cohort, cohort, type, source, to_test[1],to_test[2])
DE_SEEM = read.table(file = DE_name, header = T,sep = '\t', stringsAsFactors = F)
genes_SEEM = DE_SEEM$Gene
df_SEEM = cohort_FC_pqval_avgTPM_table_by_geneList(genes, DE_SEEM, TPM, metadata_df, Cohort, cohort_name, source, diagnosis_prm, to_test )

cohort = 'Celiac_Leonard'; cohort_name = 'Leonard'
source = 'Duodenum'
diagnosis_prm = 'Dx_specific'
to_test = c('Celiac_active','Control')
DE_name = sprintf('../../%s/res_main3TPMf_TPMf/%s_%s_%s_%s_vs_%s_deseq2_res.tsv', cohort, cohort, type, source, to_test[1],to_test[2])
DE_Leonard = read.table(file = DE_name, header = T,sep = '\t', stringsAsFactors = F)
df_Leonard = cohort_FC_pqval_avgTPM_table_by_geneList(genes, DE_Leonard, TPM, metadata_df, Cohort, cohort_name, source, diagnosis_prm, to_test )


# ## Rectal UC
Name = 'Rectal_UC'
names = c('Protect_FC','RISK_Rectal_UC_FC')

cohort = 'Protect'; cohort_name = 'Protect'
source = 'Rectal'
diagnosis_prm = 'Dx'
to_test = c('UC','Control')
DE_name = sprintf('../../%s/res_main3TPMf_TPMf/%s_%s_%s_%s_vs_%s_deseq2_res.tsv', cohort, cohort, type, source, to_test[1],to_test[2])
DE_protect = read.table(file = DE_name, header = T,sep = '\t', stringsAsFactors = F)
genes_protect = DE_protect$Gene
df_protect = cohort_FC_pqval_avgTPM_table_by_geneList(genes, DE_protect, TPM, metadata_df, Cohort, cohort_name, source, diagnosis_prm, to_test )

cohort = 'RISK'; cohort_name = 'RISK_rectalUC'
source = 'Rectal'
diagnosis_prm = 'Dx_specific'
to_test = c('UC','Control')
DE_name = sprintf('../../%s/res_main3TPMf_TPMf/%s_%s_%s_%s_vs_%s_deseq2_res.tsv', cohort, cohort, type, source, to_test[1],to_test[2])
DE_RISK_rectalUC = read.table(file = DE_name, header = T,sep = '\t', stringsAsFactors = F)
df_RISK_UC = cohort_FC_pqval_avgTPM_table_by_geneList(genes, DE_RISK_rectalUC, TPM, metadata_df, Cohort, cohort_name, source, diagnosis_prm, to_test )

# ## Rectal UC
Name = 'Ileum_CD'
names = c('SOURCE_FC','RISK_Ileum_iCD_FC')

cohort = 'SOURCE'; cohort_name = 'SOURCE'
source = 'terminal_ileum'
diagnosis_prm = 'Dx'
to_test = c('CD','Control')
DE_name = sprintf('../../%s/res_main3TPMf_TPMf/%s_%s_%s_%s_vs_%s_deseq2_res.tsv', cohort, cohort, type, source, to_test[1],to_test[2])
DE_SOURCE = read.table(file = DE_name, header = T,sep = '\t', stringsAsFactors = F)
genes_SOURCE = DE_SOURCE$Gene
df_SOURCE = cohort_FC_pqval_avgTPM_table_by_geneList(genes, DE_SOURCE, TPM, metadata_df, Cohort, cohort_name, source, diagnosis_prm, to_test )

cohort = 'RISK'; cohort_name = 'RISK_ilealCD'
source = 'ileal'
diagnosis_prm = 'Dx_specific'
to_test = c('iCD','Control')
DE_name = sprintf('../../%s/res_main3TPMf_TPMf/%s_%s_%s_%s_vs_%s_deseq2_res.tsv', cohort, cohort, type, source, to_test[1],to_test[2])
DE_RISK_ilealCD = read.table(file = DE_name, header = T,sep = '\t', stringsAsFactors = F)
df_RISK_ilealCD = cohort_FC_pqval_avgTPM_table_by_geneList(genes, DE_RISK_ilealCD, TPM, metadata_df, Cohort, cohort_name, source, diagnosis_prm, to_test )


## setting the final merged df
df = merge(x = df_SEEM, y=df_Leonard, by = 'Gene')
df = merge(x = df, y=df_protect, by = 'Gene')
df = merge(x = df, y=df_RISK_UC, by = 'Gene')
df = merge(x = df, y=df_SOURCE, by = 'Gene')
df = merge(x = df, y=df_RISK_ilealCD, by = 'Gene')

# calculate scores for different combinations
cohorts = c('Protect','RISK',
            'SOURCE','RISK_ilealCD',
            'SEEM','Leonard')

df$Score_UC_CD_Celiac = calc_score_v3(df, c('Protect', 'SOURCE','SEEM') )

df$Score_UC_CD = calc_score_v3(df, c('Protect', 'SOURCE' ) )

df2 = df
df2$SEEM__logFC = -1*df2$SEEM__logFC
df$Score_UC_CD_noCeliac = calc_score_specific(df2, c('Protect', 'SOURCE'),'SEEM' )

df2 = df
df2$SEEM__logFC = -1*df2$SEEM__logFC
df2$SOURCE__logFC = -1*df2$SOURCE__logFC
df$Score_UC_noCD_noCeliac = calc_score_specific(df2, c('Protect'),c('SEEM', 'SOURCE') )

df2 = df
df2$SEEM__logFC = -1*df2$SEEM__logFC
df2$Protect__logFC = -1*df2$Protect__logFC
df$Score_noUC_CD_noCeliac = calc_score_specific(df2, c('SOURCE'), c('Protect','SEEM')  )

df2 = df
df2$SOURCE__logFC = -1*df2$SOURCE__logFC
df2$Protect__logFC = -1*df2$Protect__logFC
df$Score_noUC_noCD_Celiac = calc_score_specific(df2, c('SEEM'),c('Protect', 'SOURCE') )

df2 = df
df2$Protect__logFC = -1*df2$Protect__logFC
df$Score_noUC_CD_Celiac = calc_score_specific(df2, c('SEEM','SOURCE'),c('Protect') )

df2 = df
df2$SOURCE__logFC = -1*df2$SOURCE__logFC
df$Score_UC_noCD_Celiac = calc_score_specific(df2, c('SEEM','Protect'),c('SOURCE') )


out_file = sprintf('../res/supps/main3_validated_genes_cohorts_%s_supp.txt',type)
write.table(x = df, file = out_file, quote = F, sep = '\t',row.names = F)
