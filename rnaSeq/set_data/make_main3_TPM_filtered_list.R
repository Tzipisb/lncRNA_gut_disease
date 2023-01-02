
type = 'lncRNA'
# type = 'ProteinCoding'
# type = 'lncPrtCd'

# TPM_file = sprintf('../data/All_cohorts_%s_TPM.txt',type)
TPM_file = sprintf('data/All_cohorts_%s_TPM_v8.txt',type)
TPM = read.table(file = TPM_file, header = T,sep = '\t', stringsAsFactors = F)

# metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v4.txt' 
# metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)

cohort = 'Protect'
TPM_filter_file = sprintf('../%s/%s_txi_%s_TPM_filtered_TPM.txt',cohort,cohort, type)
TPM_filter = read.table(file = TPM_filter_file, header = T,sep = '\t', stringsAsFactors = F)
protect_genes = TPM_filter$Gene

cohort = 'SOURCE'
TPM_filter_file = sprintf('../%s_v2/%s_txi_%s_TPM_filtered_TPM.txt',cohort,cohort, type)
TPM_filter = read.table(file = TPM_filter_file, header = T,sep = '\t', stringsAsFactors = F)
SOURCE_genes = TPM_filter$Gene

cohort = 'SEEM'
TPM_filter_file = sprintf('../%s/%s_txi_%s_TPM_filtered_TPM.txt',cohort,cohort, type)
TPM_filter = read.table(file = TPM_filter_file, header = T,sep = '\t', stringsAsFactors = F)
SEEM_genes = TPM_filter$Gene


genes = unique( c(protect_genes, SOURCE_genes, SEEM_genes) )

TPM = TPM[TPM$Gene %in% genes,]
# write.table(x = TPM, file = sprintf('res/for_yael/%s_TPM_filtered_by_main_20perTPM_v7.txt',type),quote = F, sep = '\t',row.names = F)

temp = TPM[,1]
genes_list_file = sprintf('res/gene_lists/main3_TPMfilterd_%s_list_v8.txt',type)
write.table(x = temp, file = genes_list_file, quote = F,row.names = F, col.names=F)

