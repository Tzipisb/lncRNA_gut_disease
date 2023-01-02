library(ggplot2)

type = 'lncRNA'
type = 'ProteinCoding'


cutoff = 0.1

df = data.frame()
df2 = data.frame()
pers_df = data.frame()

types = c('lncRNA','ProteinCoding')
for (type in types)
{
  tpm_file =  sprintf('data/All_cohorts_%s_TPM_v8.txt',type)
  tpm = read.table(file = tpm_file, header = T, row.names = 1)
  # tpm = tpm[as.character(genes),]
  names(tpm) = make.names(names(tpm))
  
  metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v8.txt' 
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
    pers0 = c()
    pers1 = c()
    tpm_f = tpm[,chr == metadata_df$new_chr]
    metadata_df_f = metadata_df[chr == metadata_df$new_chr,]
    metadata_df_f = metadata_df_f[metadata_df_f$Dx == 'Control',]
    tpm_f = tpm_f[, make.names(metadata_df_f$SampleID)]
    ## filter TPM to 1 TPM>0.2%
    tpm_genes = row.names(tpm_f)
    wanted_genes_pos = c()
    for ( i in 1:length(tpm_genes) )
    {
      per1 = sum(tpm_f[i,] > 1) / dim(tpm_f)[2]
      per0 = sum(tpm_f[i,] > cutoff ) / dim(tpm_f)[2]
      pers0 = c(pers0, per0)
      pers1 = c(pers1, per1)
      if ( per1 > 0.2 ) {wanted_genes_pos = c(wanted_genes_pos, i) }
    }
    # passed_filtering = i:dim(tpm_f)[1] %in% wanted_genes_pos
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
    
    t2 = data.frame(pers0 = pers0, pers1 = pers1, cohort = chr, Type = type, gene = row.names(tpm))
    pers_df = rbind(pers_df, t2)
  }
}

pers_df$cohort[pers_df$cohort == 'Protect'] = 'PROTECT'
pers_df$cohort[pers_df$cohort == 'Celiac_Leonard'] = 'PRJNA52875'
pers_df$cohort = factor(pers_df$cohort, levels = rev(c('PROTECT','SEEM','SOURCE','RISK Rectal','RISK Ileal','PRJNA52875')))

# g=ggplot(df, aes(y=Cohort, x=df$Median)) + geom_point() # + facet_wrap(.~Type, scales = 'free')
pers_df$Location = as.character(pers_df$cohort)
pers_df$Location[pers_df$cohort %in% c('PROTECT','RISK Rectal')] = 'Rectal'
pers_df$Location[pers_df$cohort %in% c('SOURCE','RISK Ileal')] = 'Ileal'
pers_df$Location[pers_df$cohort %in% c('SEEM','PRJNA52875')] = 'Duodenum'


pers_df$Location = factor(pers_df$Location, levels = c('Rectal','Ileal','Duodenum'))
pers_df$cohort = factor(pers_df$cohort, levels = c('PROTECT','SEEM','SOURCE','RISK Rectal','RISK Ileal','PRJNA52875'))
# lim = boxplot.stats(df2$tpm)$stats[c(1, 5)] + .Machine$double.eps
# lim[2] = 50000
# lim = c(0.1,1000)
# g3 = ggplot(df2) + geom_boxplot(outlier.shape = NA,aes(y=tpm, x=Cohort, fill = Type)) +
#   coord_cartesian(ylim = ( lim )* c( 1.05) ) + 
#   ylab(sprintf('TPM %s',df2$var[1])) + theme_bw() + scale_y_log10() +
#   scale_fill_brewer(palette = 'Set2') + 
#   facet_grid(.~Location, scales = 'free_x') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# pers_df_f = pers_df[pers_df$passed_filter==T,]
pers_df_f = pers_df[pers_df$pers1<=0.2,]
g4 = ggplot(pers_df_f) + geom_boxplot(aes(y=pers0*100, x=cohort, fill = Type)) +
  ylab('Percantage') + xlab('Cohort') + 
  theme_bw() + # scale_y_log10() +
  scale_fill_brewer(palette = 'Set2')  + 
  facet_grid(.~Location, scales = 'free_x') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g4

# pers_df$cut = ifelse(pers_df$pers1>0.2, 'TPM 1>20%', 'TPM 1<=20%')
pers_df$cut = ifelse(pers_df$pers1>0.2, 'Expressed', 'Discarded')
pers_df$cut = factor(pers_df$cut, levels = rev(unique(pers_df$cut)) )
g5 = ggplot(pers_df) + geom_boxplot(aes(y=pers0*100, x=cohort, fill = Type), outlier.alpha = 0) +
  ylab(sprintf('Percentage of samples\nwith TPM > %s', cutoff)) + xlab('Cohort') + 
  theme_bw() + # scale_y_log10() +
  scale_fill_brewer(palette = 'Set2')  + 
  facet_grid(cut~Location, scales = 'free_x') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
g5

out_path = 'plots/gene_expression/'
# ggsave(filename = sprintf('%s/gene_exists_by_filterPass_cohort_controlOnly.pdf', out_path), plot = g5, device = 'pdf', width = 6,height = 3.5)
# ggsave(filename = sprintf('%s/gene_exists_by_filterPass_cohort_controlOnly_withoutOutlier_exist%s.pdf', out_path, cutoff), plot = g5, device = 'pdf', width = 6,height = 3.5)

# g2=ggplot(df, aes(x=sprintf('%s+%s',df$Cohort, df$Type) )) +
#   geom_boxplot(stat="identity", aes(lower=df$`1st Qu.`,
#                                     middle=df$Median,
#                                     upper=df$`3rd Qu.`,ymin=df$`1st Qu.`,ymax=df$`3rd Qu.`)) + 
#   coord_flip() + # scale_y_log10()+ 
#   facet_wrap(.~Type, scales = 'free')


out_path = 'plots/gene_expression/'
out_file = sprintf('%s/gene_%s_%s_per_cohort_controls.tiff',out_path,'Mean','n')
# ggsave(filename = out_file,plot = g , device = 'tiff',width = 8,height = 6, compression ='lzw')

out_file3 = sprintf('%s/gene_%s_per_cohort_controls_boxplot3.pdf',out_path,df2$var[1])
# ggsave(filename = out_file3,plot = g3 , device = 'pdf',width = 6,height = 3)


for (cohort in unique(df2$Cohort))
{
  posLnc = df2$Cohort == cohort & df2$Type == 'lncRNA'
  posPc = df2$Cohort == cohort & df2$Type == 'ProteinCoding'
  
  res = wilcox.test(x = df2$tpm[posLnc], y=df2$tpm[posPc])
  
  print(cohort)
  print(res)
}



