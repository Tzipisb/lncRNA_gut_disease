## PCA
# library(car)
# library(devtools)
# library(ggbiplot)
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

type = 'lncRNA'
# type = 'ProteinCoding'
# type = 'lncPrtCd'

cohorts = c('Protect','RISK','Celiac_Nurit','Celiac_Leonard','IBD_Howell','Crohn_Yim','SEEM','SOURCE')
cohorts = c('Protect','SEEM','SOURCE')

c = 'Protect'

# type = 'ProteinCoding'
# c = 'SEEM'

# for ( c in cohorts)
{
  
  # tpm_file =  sprintf('data/All_cohorts_%s_1cTPM_filtered_TPM_joinedGenes.txt',type)
  tpm_file = sprintf('../%s/%s_txi_%s_TPM_filtered_TPM.txt',c,c,type)
  tpm = read.table(file = tpm_file, header = T, row.names = 1)
  # tpm[is.na(tpm)] = 0
  
  metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v7.txt' 
  metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)
  
  # cht = gsub(pattern = '__.*','',names(tpm))=='RISK' ## ,make sure its the same in both versions = no mistakes
  cht =  metadata_df$Cohort
  
  ids = make.names(metadata_df$SampleID[metadata_df$Cohort == c])
  # filter to wanted cohort
  tpm = tpm[,ids]
  metadata_df = metadata_df[make.names(metadata_df$SampleID) %in% ids,]
  
  # # extra filteres
  # pos = metadata_df$Source == 'ileal'
  # pos2 = metadata_df$Dx_specific == 'Control' | metadata_df$Dx_specific == 'iCD'
  # pos = pos & pos2
  # # pos = pos2
  # tpm = tpm[,pos]
  # metadata_df = metadata_df[pos,]
  
  input_data = t(tpm)
  
  df.pca = prcomp(input_data, center = T, scale = T)
  imp = summary(df.pca)$importance[2,]
  
  data = as.data.frame(df.pca$x)
  
  lvl = unique(metadata_df$Dx_specific)
  pos = which(lvl == 'Control')
  pos2 = which(lvl != 'Control')
  Dx = factor(metadata_df$Dx_specific, levels = c(lvl[pos], lvl[pos2]))
  
  # metadata_df$Cohort[metadata_df$Cohort == 'Protect'] = 'PROTECT'
  
  if (type == 'ProteinCoding')
    data$PC2 = -1*data$PC2
  
  # col_prm_name = 'Baseline Mayo'
  # col_prm = 'protect_TOTAL_MAYO_C3'
  # metadata_df$protect_TOTAL_MAYO_C3 = sprintf('UC Mayo %s',metadata_df$protect_TOTAL_MAYO_C3)
  # metadata_df$protect_TOTAL_MAYO_C3[metadata_df$Dx == 'Control'] = 'Control'
  
  # col_prm_name = 'Endoscopic Mayo'
  # col_prm = 'protect_ENDO_MAYO'
  # metadata_df$protect_ENDO_MAYO = sprintf('UC Mayo %s',metadata_df$protect_ENDO_MAYO)
  # metadata_df$protect_ENDO_MAYO[metadata_df$Dx == 'Control'] = 'Control'
  # 
  # col_prm = 'protect_W4R'
  # col_prm_name = 'Week 4 remission'
  # metadata_df$protect_W4R = sprintf('UC %s',metadata_df$protect_W4R)
  # metadata_df$protect_W4R[metadata_df$Dx == 'Control'] = 'Control'

  # col_prm = 'protect_3years_colectomy'
  # col_prm_name = '3 years colectomy'
  # metadata_df$protect_3years_colectomy = sprintf('UC %s',metadata_df$protect_3years_colectomy)
  # metadata_df$protect_3years_colectomy[metadata_df$Dx == 'Control'] = 'Control'

  col_prm_name = 'Baseline PUCAI'
  col_prm = 'protect_BASELINE_Pucai_Category'
  metadata_df$protect_BASELINE_Pucai_Category = sprintf('UC PUCAI %s',metadata_df$protect_BASELINE_Pucai_Category)
  metadata_df$protect_BASELINE_Pucai_Category[metadata_df$Dx == 'Control'] = 'Control'
  
  g = ggplot(data, aes(x=PC1, y=PC2)) +
    # geom_point(shape = 21, aes( fill = metadata_df$protect_BASELINE_Pucai_Category ), size=2) +
    geom_point(shape = 21, aes( fill = metadata_df$protect_BASELINE_Pucai_Category ), size=3) +
    xlab(sprintf('PC1 (%.2f%%)', imp[1]*100)) + ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) +
    # scale_fill_manual(values = c('#0066CC','#FF0033')) + 
    scale_fill_manual(values = c('#0066CC','gold','darkorange1','firebrick')) +
    # scale_fill_manual(values = c('#0066CC','gold','firebrick')) +
    guides(fill=guide_legend(title=col_prm_name, override.aes = list(shape = 21))) + 
    theme_bw()
  # + scale_colour_brewer(palette = 'Paired')
  g
  # g = ggplot(data, aes(x=data$PC1, y=data$PC2)) +
  #   geom_point(aes( colour = metadata_df$Source, shape = metadata_df$Source)) # + scale_colour_brewer(palette = 'Paired')
  
  # out_path = sprintf('plots/%s/',type)
  out_path = sprintf('plots/PCA/')
  dir.create(out_path)
  # ggsave(filename = sprintf('%s/PCA_%s_%s_dx_%s.pdf', out_path, c, type, col_prm_name), plot = g, device = 'pdf', width = 6,height = 4)
  
  b = ggplot(data) + 
    # geom_boxplot(aes( x = as.factor(metadata_df$protect_BASELINE_Pucai_Category), 
    #                   ill = as.factor(metadata_df$protect_BASELINE_Pucai_Category), y=PC2)) +
    geom_boxplot(aes( x = as.factor(metadata_df[[col_prm]]),
                      fill = as.factor(metadata_df[[col_prm]]), y=PC2)) +
    scale_fill_manual(values = c('#0066CC','gold','darkorange1','firebrick')) +
    theme_bw() +
    # scale_fill_manual(values = c('#0066CC','gold','firebrick')) +
    theme(legend.position = 'none',axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    xlab(col_prm_name)
  res = add_significant_asterix_to_plot_BH_v2(p = b, var = as.factor(metadata_df[[col_prm]]), val = data$PC2, test_type = 'wilcox', print_pvals = T)
  b = res[[1]]
  # ggsave(filename = sprintf('%s/PCA_%s_%s_dx_%s_PC2_boxplot.pdf', out_path, c, type, col_prm_name), plot = b, device = 'pdf', width = 4,height = 4)
  
}
