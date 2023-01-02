## PCA
# library(car)
# library(devtools)
# library(ggbiplot)
library(ggplot2)

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
  
  # g = ggplot(data, aes(x=PC1, y=PC2)) +
  #   geom_point(aes( colour = Dx )) +
  #   scale_color_manual(values = c('#0066CC','#FF0033')) + theme_bw()
  #   # + scale_colour_brewer(palette = 'Paired')
  
  g = ggplot(data, aes(x=PC1, y=PC2)) +
    geom_point(shape = 21, aes( fill = Dx ), size=3) +
    xlab(sprintf('PC1 (%.2f%%)', imp[1]*100)) + ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) +
    scale_fill_manual(values = c('#0066CC','#FF0033')) + theme_bw()
  # + scale_colour_brewer(palette = 'Paired')
  g
  # g = ggplot(data, aes(x=data$PC1, y=data$PC2)) +
  #   geom_point(aes( colour = metadata_df$Source, shape = metadata_df$Source)) # + scale_colour_brewer(palette = 'Paired')
  
  # out_path = sprintf('plots/%s/',type)
  out_path = sprintf('plots/PCA/')
  dir.create(out_path)
  # ggsave(filename = sprintf('%s/PCA_%s_%s_dx.pdf', out_path, c, type), plot = g, device = 'pdf', width = 6,height = 4)
  
  # cs = c(1,2)
  # g <- ggbiplot(df.pca, obs.scale = 1, var.scale = 1,
  #               groups = metadata_df$Dx_specific, ellipse = F,
  #               circle = F, choices = cs,  labels = metadata_df$SampleID,
  #               varname.size = 2.5, varname.adjust = 1.1, var.axes = F)
  
}