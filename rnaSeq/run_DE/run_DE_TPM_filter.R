library('DESeq2')
source('/pita/users/tzipi/code/rnaSeq/DESeq2_DE_funcs.r')

filter_by_TPM_flag = TRUE
# filter_by_TPM_flag = FALSE

filter_2_main_genes = TRUE
# filter_2_main_genes = FALSE

type = 'lncRNA'
# type = 'lncPrtCd'
type = 'ProteinCoding'

# cohort = 'Protect'
# source = 'Rectal'
# diagnosis_prm = 'Dx'
# to_test = c('UC','Control')

# cohort = 'Celiac_Leonard'
# source = 'Duodenum'
# diagnosis_prm = 'Dx_specific'
# to_test = c('Celiac_active','Control')

cohort = 'RISK'
diagnosis_prm = 'Dx_specific'
source = 'ileal'
to_test = c('iCD','Control')
# source = 'Rectal'
# to_test = c('UC','Control')

# cohort = 'SEEM'
# diagnosis_prm = 'Dx_specific'
# source = 'Doudenum'
# to_test = c('Celiac','Control')

# cohort = 'Crohn_Yim'
# diagnosis_prm = 'Dx_specific'
# source = 'ileal_fibroblast'
# to_test = c('CD_inflamed','Control')

# cohort = 'IBD_Howell'
# diagnosis_prm = 'Dx'
# source = 'terminal_ileum'
# to_test = c('CD','Control')
# # source = 'sigmoid_colon'
# # to_test = c('UC','Control')

# cohort = 'SOURCE'
# source = 'terminal_ileum'
# diagnosis_prm = 'Dx'
# to_test = c('CD','Control')

metadata_file =  '../../metadata/lncRNA_meta_analysis_megamap_v4.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)
metadata_df = metadata_df[order(metadata_df$SampleID),]
metadata_df = metadata_df[metadata_df$Cohort == cohort,]

pos_source = metadata_df$Source == source 
pos_dx = metadata_df[[diagnosis_prm]] == to_test[1] | metadata_df[[diagnosis_prm]] == to_test[2]
pos = pos_source & pos_dx
metadata_df = metadata_df[pos,]

rds_file = sprintf('../../%s/%s_txi_%s.rds',cohort, cohort, type)

txi = readRDS(rds_file)
txi$abundance = txi$abundance[,pos]
txi$counts = txi$counts[,pos]
txi$length = txi$length[,pos]

# only uses the genes that passes TPM filtering in one of the main 2 cohorts.
if ( filter_2_main_genes )
{
  gene_list = read.table(sprintf('../res/gene_lists/main3_TPMfilterd_%s_list.txt', type))  
  txi$abundance = txi$abundance[row.names(txi$abundance) %in% gene_list$V1, ]
  txi$counts = txi$counts[row.names(txi$counts) %in% gene_list$V1, ]
  txi$length = txi$length[row.names(txi$length) %in% gene_list$V1, ]
}

if ( filter_by_TPM_flag )
{
  txi_genes = row.names(txi$abundance)
  wanted_genes_pos = c()
  for ( i in 1:length(txi_genes) )
  {
    per = sum(txi$abundance[i,] > 1) / dim(txi$counts)[2]
    if ( per > 0.2 ) {wanted_genes_pos = c(wanted_genes_pos, i) }
  }
  txi$abundance = txi$abundance[wanted_genes_pos, ]
  txi$counts = txi$counts[wanted_genes_pos, ]
  txi$length = txi$length[wanted_genes_pos, ]
}

diagnosis = factor(as.factor(metadata_df[[diagnosis_prm]]), levels=rev(to_test))
coldata = data.frame(row.names=metadata_df$Samples, diagnosis=diagnosis)
# coldata$row.names=metadata_df$Samples
# colData(ddsHTSeq)$diagnosis<-factor(colData(ddsHTSeq)$diagnosis, levels=rev(to_test) )

ddsHTSeq <- DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~ diagnosis)

dds<-DESeq(ddsHTSeq)

res<-results(dds)
res<-res[order(res$padj),]
head(res)

res_name = sprintf('%s_%s_%s_%s_vs_%s', cohort, type, source, to_test[1],to_test[2])
out_path = sprintf('../../%s/res_v2/', cohort)
if (filter_2_main_genes & filter_by_TPM_flag)
{
  out_path = sprintf('../../%s/res_main3TPMf_TPMf/', cohort)
}
dir.create(out_path) # v2 = TPM filtered for only the genes in this DE
write_DE_result(dds, name = res_name, out_path = out_path, type_flag = F)
res_df = write_DE_result(dds, name = res_name, filter_flag = T, FDR_cutoff = 0.05, FC_cutoff = log2(1.5), out_path = out_path, type_flag = F)

