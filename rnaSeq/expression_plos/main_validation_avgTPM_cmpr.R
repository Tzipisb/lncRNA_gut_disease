library(ggplot2)
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

type = 'lncRNA'
# type = 'ProteinCoding'

cohorts = c('PROTECT','SEEM','SOURCE','RISK_UC_Rectum','RISK_iCD_Ileum','Celiac_Leonard')

if (!exists('disease'))
  disease = 'UC'

if ( disease == 'UC' )
{
  x_var ='PROTECT__Control'
  y_var = 'RISK_UC_Rectum__Control'
  genes_file = sprintf('plots/venns/lists_v2/PROTECT_%s_TPM_control.txt',type)
}
if ( disease == 'CD' )
{
  x_var ='SOURCE__Control'
  y_var = 'RISK_iCD_Ileum__Control'
  genes_file = sprintf('plots/venns/lists_v2/SOURCE_%s_TPM_control.txt',type)
}
if ( disease == 'Celiac' )
{
  x_var ='SEEM__Control'
  y_var = 'Celiac_Leonard__Control'
  genes_file = sprintf('plots/venns/lists_v2/SEEM_%s_TPM_control.txt',type)
}

tpm_file =  sprintf('data/All_cohorts_%s_TPM_v7.txt',type)
tpm = read.table(file = tpm_file, header = T, row.names = 1)
# tpm[is.na(tpm)] = 0

metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v7.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))

metadata_df = metadata_df[order(metadata_df$SampleID),]
tpm = tpm[, make.names( metadata_df$SampleID )]

# ## filter to genes that passed main3 TPM filtering
# genes_file = sprintf('res/gene_lists/main3_TPMfilterd_%s_list.txt',type)
# g_df = read.table(file = genes_file, header = F,sep = '\t')
# genes = g_df$V1
# tpm = tpm[as.character(genes),]

g_df = read.table(file = genes_file, header = F,sep = '\t')
genes = g_df$V1
tpm = tpm[as.character(genes),]


new_chr = metadata_df$Cohort
new_chr[metadata_df$Cohort == 'RISK' &  
          metadata_df$Source == 'Rectal' & 
          metadata_df$Dx %in% c('UC','Control')] = 'RISK_UC_Rectum'
new_chr[metadata_df$Cohort == 'RISK' &  
          metadata_df$Source == 'ileal' & 
          metadata_df$Dx_specific %in% c('iCD','Control')] = 'RISK_iCD_Ileum'
new_chr[new_chr == 'Protect'] = 'PROTECT'
metadata_df$new_chr = new_chr

healthy = metadata_df$Dx
healthy[healthy != 'Control'] = 'Disease'
metadata_df$healthy = healthy

df = data.frame(row.names = row.names(tpm))
for (chr in cohorts)
{
  for ( dx in c('Control','Disease') )
  {
    temp = tpm[,chr == metadata_df$new_chr & metadata_df$healthy == dx]
    avg_tpm = rowMeans(temp)
    df[[ sprintf('%s__%s',chr,dx) ]] = avg_tpm
  }
}

df_ctrl = df[,grepl(pattern = 'Control',x = names(df))]


c = gsub(pattern = '__.*','',names(df_ctrl))

g = plot_variable_and_stats(taxa =df_ctrl, check_val = x_var, boxplot_test_type = 'spearman',
                            y_val =y_var, y_lb = y_var, sig_ast_flag = T, print_pvals = F)
g = g +
  scale_x_log10() + scale_y_log10() +
  xlab(sprintf('%s Avg TPM',x_var)) + ylab(sprintf('%s Avg TPM',y_var))
out_path = 'plots/gene_expression/control_tpm_noise/'
# ggsave(sprintf('%s/control_tpm_%s_%s.tiff', out_path, x_var, y_var),plot = g, device = 'tiff', 
#        width = 6,height = 4, compression = 'lzw')
res = cor.test(x = df_ctrl[[x_var]], y=df_ctrl[[y_var]], method = 'spearman')

