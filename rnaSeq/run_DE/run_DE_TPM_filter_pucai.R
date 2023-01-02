library('DESeq2')
source('/pita/users/tzipi/code/rnaSeq/DESeq2_DE_funcs.r')

filter_by_TPM_flag = TRUE
# filter_by_TPM_flag = FALSE

cohort = 'Protect'

# source = 'Rectal'
# diagnosis_prm = 'Dx'
# to_test = c('UC','Control')

source = 'Rectal'

diagnosis_prm = 'protect_BASELINE_Pucai_Category'
to_test = c('1','3')

# diagnosis_prm = 'protect_TOTAL_MAYO_C3'
# to_test = c('1','3')

# diagnosis_prm = 'protect_ENDO_MAYO'
# to_test = c('1','3')

# diagnosis_prm = 'protect_W4R'
# to_test = c('0','1')

# diagnosis_prm = 'protect_3years_colectomy'
# to_test = c('0','1')

# diagnosis_prm = 'protect_REMISSION_WK52'
3 to_test = c('0','1')

type = 'lncRNA'
# type = 'lncPrtCd'
# type = 'ProteinCoding'

metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v7.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)
metadata_df = metadata_df[order(metadata_df$SampleID),]
metadata_df = metadata_df[metadata_df$Cohort == cohort,]

pos_source = metadata_df$Source == source 
pos_dx = metadata_df[[diagnosis_prm]] == to_test[1] | metadata_df[[diagnosis_prm]] == to_test[2]
pos = pos_source & pos_dx
metadata_df = metadata_df[pos,]

rds_file = sprintf('%s_txi_%s.rds',cohort, type)

txi = readRDS(rds_file)
txi$abundance = txi$abundance[,pos]
txi$counts = txi$counts[,pos]
txi$length = txi$length[,pos]

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

res_name = sprintf('%s_%s_%s_%s_vs_%s', cohort, type, diagnosis_prm, to_test[1],to_test[2])
out_path = 'res_PUCAI/'
# dir.create(out_path) # v2 = TPM filtered for only the genes in this DE
write_DE_result(dds, name = res_name, out_path = out_path, type_flag = F)
res_df = write_DE_result(dds, name = res_name, filter_flag = T, FDR_cutoff = 0.05, FC_cutoff = log2(1.5), out_path = out_path, type_flag = F)

