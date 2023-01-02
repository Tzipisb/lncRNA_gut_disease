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
# {
  
  # tpm_file =  sprintf('data/All_cohorts_%s_1cTPM_filtered_TPM_joinedGenes.txt',type)
  tpm_file = sprintf('../%s/%s_txi_%s_TPM_filtered_TPM.txt',c,c,type)
  tpm = read.table(file = tpm_file, header = T, row.names = 1, na.strings = c('NA','na'))
  # tpm[is.na(tpm)] = 0
  
  metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v7.txt' 
  metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))
  
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
  
  data = as.data.frame(df.pca$x)
  
  lvl = unique(metadata_df$Dx_specific)
  pos = which(lvl == 'Control')
  pos2 = which(lvl != 'Control')
  Dx = factor(metadata_df$Dx_specific, levels = c(lvl[pos], lvl[pos2]))
  
  out_path = sprintf('plots/PCA/')

# }
  
protect_metadata = c('Age_years','Gender','protect_BL_BMIZ',
                       'Dx','protect_BASELINE_Pucai',
                       'protect_BASELINE_Pucai_Category',
                       'protect_ENDO_MAYO',
                       'protect_TOTAL_MAYO_EEF',
                       'protect_TOTAL_MAYO_C3',
                       'protect_PARIS_CLASS',
                       # 'protect_Nancy_total_histology',
                       # 'protect_CBH_VILLIFORM',
                       # 'protect_CBH_INFLAM',
                       # 'protect_CBH_EOSCOUNT',
                       # 'protect_EOSGRADE_C2',
                       'protect_ALB_C2',
                       'protect_W4R',
                       # 'protect_W12_CSfreeR',
                       'protect_CSFREE_REMISSION_WK52_KM',
                       'protect_3years_colectomy')
mets = metadata_df[,protect_metadata]

mets$Dx = ifelse(mets$Dx=="Control", 0, 1)
mets$Gender = ifelse(mets$Gender=="Female", 2, 1)
mets$Age_years = as.numeric(mets$Age_years)  
mets$protect_PARIS_CLASS = gsub(pattern = 'E',replacement = '',x = mets$protect_PARIS_CLASS)
mets$protect_PARIS_CLASS = as.numeric(mets$protect_PARIS_CLASS)
mets$protect_3years_colectomy = as.numeric(mets$protect_3years_colectomy)
mets$protect_CSFREE_REMISSION_WK52_KM = 1 - as.numeric(mets$protect_CSFREE_REMISSION_WK52_KM)
mets$protect_W4R = 1 - as.numeric(mets$protect_W4R)
YN_prms = c('protect_ALB_C2')
for ( prm in YN_prms )
{
  mets[[prm]] = ifelse(mets[[prm]]=="Y", 1, 0)
}
colnames(mets)[colnames(mets) == 'Age_years'] = 'Age'
colnames(mets)[colnames(mets) == 'protect_BL_BMIZ'] = 'BMI'
colnames(mets)[colnames(mets) == 'Dx'] = 'protect_Ulcerative colitis'
colnames(mets)[colnames(mets) == 'protect_BASELINE_Pucai'] = 'protect_PUCAI'
colnames(mets)[colnames(mets) == 'protect_BASELINE_Pucai_Category'] = 'protect_PUCAI category'
colnames(mets)[colnames(mets) == 'protect_ENDO_MAYO'] = 'protect_Endoscopic Mayo'
colnames(mets)[colnames(mets) == 'protect_TOTAL_MAYO_EEF'] = 'protect_Total Mayo'
colnames(mets)[colnames(mets) == 'protect_TOTAL_MAYO_C3'] = 'protect_Total Mayo category'
colnames(mets)[colnames(mets) == 'protect_PARIS_CLASS'] = 'protect_Paris classification'
# colnames(mets)[colnames(mets) == 'protect_Nancy_total_histology'] = 'Baseline Nancy histology score'
# colnames(mets)[colnames(mets) == 'protect_CBH_VILLIFORM'] = 'Baseline histology villiform changes'
# colnames(mets)[colnames(mets) == 'protect_CBH_INFLAM'] = 'Baseline histology inflammation'
# colnames(mets)[colnames(mets) == 'protect_CBH_EOSCOUNT'] = 'Baseline histology eosinophils count'
# colnames(mets)[colnames(mets) == 'protect_EOSGRADE_C2'] = 'Baseline histology eosinophils grade'
colnames(mets)[colnames(mets) == 'protect_ALB_C2'] = 'protect_Albumin'
colnames(mets)[colnames(mets) == 'protect_W4R'] = 'protect_W4 no response'
# colnames(mets)[colnames(mets) == 'protect_W12_CSfreeR'] = 'protect_W12 CS free remsion '
colnames(mets)[colnames(mets) == 'protect_CSFREE_REMISSION_WK52_KM'] = 'protect_W52 no CS free remssion '
colnames(mets)[colnames(mets) == 'protect_3years_colectomy'] = 'protect_Colectomy within 3Y'

names(mets) = gsub('protect_','',names(mets))

cor_r_p1 = vector(mode = 'numeric',length = dim(mets)[2] )
cor_p_p1 = vector(mode = 'numeric',length = dim(mets)[2] )
cor_r_p2 = vector(mode = 'numeric',length = dim(mets)[2] )
cor_p_p2 = vector(mode = 'numeric',length = dim(mets)[2] )
for ( i in 1:dim(mets)[2] )
{
  res_p1 = cor.test(x = data$PC1, y=mets[,i], method = 'spearman')
  res_p2 = cor.test(x = data$PC2, y=mets[,i], method = 'spearman')
  cor_r_p1[i] = res_p1$estimate
  cor_r_p2[i] = res_p2$estimate
  cor_p_p1[i] = res_p1$p.value
  cor_p_p2[i] = res_p2$p.value
}
df = data.frame(cor_r_p1 = cor_r_p1, cor_r_p2 = cor_r_p2, cor_p_p1 = cor_p_p1, cor_p_p2 = cor_p_p2, parameter = names(mets))
# df = df[ ( df$cor_p_p1 <= 0.05 & abs(df$cor_r_p1) > 0.2)  |  ( df$cor_p_p2 <= 0.05 & abs(df$cor_r_p2) > 0.2), ]
df = df[ df$cor_p_p1 <= 0.05 | df$cor_p_p2 <= 0.05, ]

temp = abs(df$cor_p_p1)+abs(df$cor_p_p2)
df$parameter = factor(df$parameter, levels = df$parameter[order(temp)] )
ggplot(df, aes(x=cor_r_p1, y=cor_r_p2)) + geom_point()


# df$Angle = ((180/pi) * atan(df$cor_r_p2/df$cor_r_p1))
# # df$Angle[df$Angle<0] = df$Angle[df$Angle<0]+180
# # df$Offset <- ((-2 * sign(df$cor_r_p1))/2) 
# g = ggplot(df) + geom_segment(data = df, aes(x = 0, y = 0, xend = (cor_r_p1), yend = (cor_r_p2)),
#              arrow = arrow(length = unit(1/2, "picas")), color = "darkgray") + theme_void() + 
#   # geom_point(aes(x=cor_r_p1, y=cor_r_p2), colour='red') + 
#   # geom_text(data = df, aes(label = parameter, x = (cor_r_p1*(4/5)), y = (cor_r_p2)*(4/5)), 
#   geom_text(data = df, aes(label = parameter, x = (cor_r_p1)*(2/3), y = (cor_r_p2)*(2/3)), 
#             color = "black", size = 4,  angle = df$Angle)#, vjust = df$Offset )

unique_cols = c('#e6194b', '#3cb44b', '#ffe119',
                '#4363d8', '#f58231', '#911eb4', '#46f0f0',
                '#f032e6', '#bcf60c', '#9a6324', '#008080',
                '#e6beff', '#aaffc3', '#800000', '#fabebe',
                '#fffac8', '#808000', '#ffd8b1', '#000075',
                '#808080', '#000000') # '#ffffff' = white
# unique_cols = c('#e6194b', '#3cb44b', '#ffe119', 
#                 '#4363d8', '#f58231', '#911eb4', '#46f0f0', 
#                 '#f032e6', '#bcf60c', '#fabebe', '#008080', 
#                 '#e6beff', '#9a6324', '#808000', '#800000', 
#                 '#fffac8', '#aaffc3', '#ffd8b1', '#000075', 
#                 '#808080', '#000000') # '#ffffff' = white
gc = ggplot(df) + geom_segment(data = df, aes(x = 0, y = 0, xend = (cor_r_p1), yend = (cor_r_p2), colour=parameter),size = 1,
                              arrow = arrow(length = unit(1, "picas"))) + theme_bw() + 
  xlab('PC1 corelation R') + ylab('PC2 corelation R') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_colour_manual(name = '',values = unique_cols) # + 
  # geom_point(aes(x=cor_r_p1, y=cor_r_p2,colour = parameter)) # + 
  # geom_text(data = df, aes(label = parameter, x = (cor_r_p1*(4/5)), y = (cor_r_p2)*(4/5)), 
  # geom_text(data = df, aes(label = parameter, x = (cor_r_p1), y = (cor_r_p2)), http://172.35.15.24http://172.35.15.242:8798/graphics/plot_zoom_png?width=541&height=3172:8798/graphics/plot_zoom_png?width=541&height=317
  #           color = "black", size = 4,  angle = df$Angle)#, vjust = df$Offset )


out_path = 'plots/PCA/'
# ggsave(sprintf('%s/PCA_corrs_%s_%s.pdf', out_path,c,type),plot = g, device = 'pdf', width = 5,height = 3)

