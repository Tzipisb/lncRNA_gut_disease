library(ggplot2)

type = 'lncRNA'
type = 'ProteinCoding'

df = data.frame()
df2 = data.frame()

types = c('lncRNA','ProteinCoding')
for (type in types)
{
  tpm_file =  sprintf('data/All_cohorts_%s_TPM_v7.txt',type)
  tpm = read.table(file = tpm_file, header = T, row.names = 1)
  # tpm = tpm[as.character(genes),]
  
  metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v7.txt' 
  metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)
  metadata_df = metadata_df[order(metadata_df$SampleID),]
  tpm = tpm[, make.names( metadata_df$SampleID )]
  
  cohorts = unique(metadata_df$Cohort)
  
  cohorts = c('Protect','SEEM','SOURCE','RISK Rectal','RISK Ileal','Celiac_Leonard')
  new_chr = metadata_df$Cohort
  new_chr[metadata_df$Cohort == 'RISK' &  
            metadata_df$Source == 'Rectal' & 
            metadata_df$Dx %in% c('UC','Control')] = 'RISK Rectal'
  new_chr[metadata_df$Cohort == 'RISK' &  
            metadata_df$Source == 'ileal' & 
            metadata_df$Dx_specific %in% c('iCD','Control')] = 'RISK Ileal'
  metadata_df$new_chr = new_chr
  for (chr in cohorts)
  {
    tpm_f = tpm[,chr == metadata_df$new_chr]
    metadata_df_f = metadata_df[chr == metadata_df$new_chr,]
    metadata_df_f = metadata_df_f[metadata_df_f$Dx == 'Control',]
    tpm_f = tpm_f[, make.names(metadata_df_f$SampleID)]
    ## filter TPM to 1 TPM>0.2%
    tpm_genes = row.names(tpm_f)
    wanted_genes_pos = c()
    for ( i in 1:length(tpm_genes) )
    {
      per = sum(tpm_f[i,] > 1) / dim(tpm_f)[2]
      if ( per > 0.2 ) {wanted_genes_pos = c(wanted_genes_pos, i) }
    }
    tpm_f = tpm_f[wanted_genes_pos, ]
    
    smr=summary(apply(tpm_f,1,mean))
    # smr=summary(apply(tpm_f,1,var))
    smr = data.frame(t(unclass(smr)), check.names = FALSE, stringsAsFactors = FALSE)
    smr$n = dim(tpm_f)[1]
    smr$Cohort = chr
    smr$Type = type
    # row.names(smr)[1] = sprintf('%s_%s',chr,type)
    df = rbind(df, smr)
    
    tpm_med = apply(tpm_f,1,mean)
    t=data.frame(tpm = tpm_med, var = 'Mean',Cohort = chr, Type = type)
    df2 = rbind(df2, t)
    
  }
}



df2$Cohort[df2$Cohort == 'Celiac_Leonard'] = 'PRJNA52875'
df2$Cohort[df2$Cohort == 'Protect'] = 'PROTECT'
df2$Cohort = factor(df2$Cohort, levels = rev(c('PROTECT','SEEM','SOURCE','RISK Rectal','RISK Ileal','PRJNA52875')))

# g=ggplot(df, aes(y=Cohort, x=df$Median)) + geom_point() # + facet_wrap(.~Type, scales = 'free')
df2$Location = as.character(df2$Cohort)
df2$Location[df2$Cohort %in% c('PROTECT','RISK Rectal')] = 'Rectal'
df2$Location[df2$Cohort %in% c('SOURCE','RISK Ileal')] = 'Ileal'
df2$Location[df2$Cohort %in% c('SEEM','PRJNA52875')] = 'Duodenum'


df2$Location = factor(df2$Location, levels = c('Rectal','Ileal','Duodenum'))
df2$Cohort = factor(df2$Cohort, levels = c('PROTECT','SEEM','SOURCE','RISK Rectal','RISK Ileal','PRJNA52875'))
# lim = boxplot.stats(df2$tpm)$stats[c(1, 5)] + .Machine$double.eps
# lim[2] = 50000
lim = c(0.3,1200)
g3 = ggplot(df2) + geom_boxplot(outlier.shape = NA,aes(y=tpm, x=Cohort, fill = Type)) +
  coord_cartesian(ylim = ( lim )* c( 1.05) ) + 
  ylab(sprintf('TPM %s',df2$var[1])) + theme_bw() + scale_y_log10() +
  scale_fill_brewer(palette = 'Set2') + 
  facet_grid(.~Location, scales = 'free_x') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


g3
# g2=ggplot(df, aes(x=sprintf('%s+%s',df$Cohort, df$Type) )) +
#   geom_boxplot(stat="identity", aes(lower=df$`1st Qu.`,
#                                     middle=df$Median,
#                                     upper=df$`3rd Qu.`,ymin=df$`1st Qu.`,ymax=df$`3rd Qu.`)) + 
#   coord_flip() + # scale_y_log10()+ 
#   facet_wrap(.~Type, scales = 'free')

out_path = 'plots/gene_expression/'
out_file = sprintf('%s/gene_%s_%s_per_cohort_controls.tiff',out_path,'Mean','n')
# ggsave(filename = out_file,plot = g , device = 'tiff',width = 8,height = 6, compression ='lzw')

out_file3 = sprintf('%s/gene_%s_per_cohort_controls_fixed_boxplot3.pdf',out_path,df2$var[1])
# ggsave(filename = out_file3,plot = g3 , device = 'pdf',width = 6,height = 3)


for (cohort in unique(df2$Cohort))
{
  posLnc = df2$Cohort == cohort & df2$Type == 'lncRNA'
  posPc = df2$Cohort == cohort & df2$Type == 'ProteinCoding'

  res = wilcox.test(x = df2$tpm[posLnc], y=df2$tpm[posPc])

  print(cohort)
  print(res)
}

