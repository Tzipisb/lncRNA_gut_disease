
cohorts = c('Protect','RISK','Celiac_Nurit','Celiac_Leonard','IBD_Howell','Crohn_Yim','SOURCE','SEEM','caco','blood')

type = 'lncRNA'
# type = 'ProteinCoding'
type = 'lncPrtCd'

cohort = cohorts[1]
tpm_file_format = '../%s/%s_txi_%s_TPM.txt'
# tpm_file_format = '../%s/%s_txi_%s_TPM_filtered_TPM.txt'
tpm_file =  sprintf(tpm_file_format,cohort,cohort, type)
# tpm_file =  sprintf('../%s/%s_txi_%s_TPM_filtered_TPM.txt',cohort,cohort, type)
merged_tpm = read.table(file = tpm_file, header = T)
for ( i in 2:length(cohorts) )
{
  cohort = cohorts[i]
  tpm_file =  sprintf(tpm_file_format,cohort,cohort, type)
  # if (cohort == 'RISK')
  # {
  #   tpm_file =  sprintf('../%s/%s_s31_txi_%s_TPM_filtered_TPM.txt',cohort,cohort, type)
  # }
  tpm2 = read.table(file = tpm_file, header = T)

  merged_tpm = merge(x = merged_tpm, y = tpm2, by = 'Gene', all=T)
  # merged_tpm = merge(x = merged_tpm, y = tpm2, by = 'Gene')
}

# names(merged_tpm)

# out_file = sprintf('data/All_cohorts_%s_1cTPM_filtered_TPM.txt', type)
out_file = sprintf('data/All_cohorts_%s_TPM_v7.txt', type)
# out_file = sprintf('data/All_cohorts_%s_1cTPM_filtered_TPM_s31Try.txt', type)
# out_file = sprintf('data/All_cohorts_%s_1cTPM_filtered_TPM_joinedGenes.txt', type)


write.table(x = merged_tpm, file = out_file,quote = F,sep = '\t',row.names = F)




