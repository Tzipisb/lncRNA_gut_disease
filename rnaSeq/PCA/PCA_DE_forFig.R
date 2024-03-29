source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
## PCA
library(ggplot2)
library(stats)

type = 'lncRNA'
# type = 'ProteinCoding'
# type = 'lncPrtCd'

if (!exists('diseas'))
  diseas = 'UC'

if ( diseas == 'UC')
{
  c = c('Protect','RISK')
  cohort = 'Protect'
  source = 'Rectal'
  test_source = source
  diagnosis_prm = 'Dx'
  to_test = c('UC','Control')
}

if ( diseas == 'CD')
{
  c= c('SOURCE', 'RISK')
  cohort = 'SOURCE'
  diagnosis_prm = 'Dx_specific'
  source = 'terminal_ileum'
  test_source = c('terminal_ileum','ileal')
  to_test = c('CD','Control','iCD')
}

if ( diseas == 'Celiac')
{
  c = c('SEEM','Celiac_Leonard')
  cohort = 'SEEM'
  diagnosis_prm = 'Dx_specific'
  source = 'Duodenum'
  test_source = source
  to_test = c('Celiac','Control','Celiac_active')
}

# cohort = 'RISK'
# diagnosis_prm = 'Dx_specific'
# source = 'ileal'
# to_test = c('iCD','Control')
# source = 'Rectal'
# to_test = c('UC','Control')
# 
# cohort = 'Crohn_Yim'
# diagnosis_prm = 'Dx_specific'
# source = 'ileal_fibroblast'
# to_test = c('CD_inflamed','Control')
# to_test = c('CD_non_inflamed','Control')
# to_test = c('CD_stenotic','Control')
# to_test = c('CD_inflamed','CD_non_inflamed')
# to_test = c('CD_inflamed','CD_stenotic')
# to_test = c('CD_non_inflamed','CD_stenotic')
# 
# cohort = 'IBD_Howell'
# diagnosis_prm = 'Dx'
# source = 'terminal_ileum'
# to_test = c('CD','Control')
# 
# cohort = 'Celiac_Nurit'
# diagnosis_prm = 'Dx_specific'
# source = 'Duodenum_FFPE'
# to_test = c('Celiac','Control')
# 
# cohort = 'Celiac_Leonard'
# diagnosis_prm = 'Dx_specific'
# source = 'Duodenum'
# # to_test = c('Celiac_remission','Control')
# to_test = c('Celiac_active','Control')
# # to_test = c('Celiac_active','Celiac_remission')

DE_name = sprintf('%s_%s_%s_%s_vs_%s', gsub('SOURCE','SOURCE_v2',cohort), type, source, to_test[1],to_test[2])
DE_file = sprintf('../%s/res_main3TPMf_TPMf_v8/%s_deseq2_res_filtered.tsv',cohort, DE_name)
DE_df = read.table(file = DE_file, header = T,sep = '\t')
genes = DE_df$Gene   

tpm_file =  sprintf('data/All_cohorts_%s_TPM_v8.txt',type)
tpm = read.table(file = tpm_file, header = T, row.names = 1)
tpm = tpm[as.character(genes),]
names(tpm) = make.names(names(tpm))

metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v8.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)
metadata_df = metadata_df[order(metadata_df$SampleID),]
tpm = tpm[, make.names( metadata_df$SampleID )]

# cht = gsub(pattern = '__.*','',names(tpm))=='RISK' ## ,make sure its the same in both versions = no mistakes
cht =  metadata_df$Cohort

#             #### choose data to look at (usign DE genes)
# # filter to wanted cohort
tpm = tpm[,cht %in% c]
metadata_df = metadata_df[cht %in% c,]


# extra filteres
pos = metadata_df$Source %in% test_source
pos2 = metadata_df[[diagnosis_prm]] %in% to_test
pos = pos & pos2 
# pos = pos2
tpm = tpm[,pos]
metadata_df = metadata_df[pos,]


# tpm2 = tpm[rowSums(is.na(tpm)) != ncol(tpm),] # remove NAs lines
tpm2 = tpm[rowSums(is.na(tpm)) == 0,] # remove any line with >1 NA
tpm2 = tpm2[rowSums(tpm2) != 0,] # remove any gene with 0 expression
input_data = t(tpm2)

df.pca = prcomp(input_data, center = T, scale = T)
imp = summary(df.pca)$importance[2,]

data = as.data.frame(df.pca$x)

lvl = unique(metadata_df$Dx)
pos = which(lvl == 'Control')
pos2 = which(lvl != 'Control')
metadata_df$Dx = factor(metadata_df$Dx, levels = c(lvl[pos], lvl[pos2]))

# g = ggplot(data, aes(x=data$PC1, y=data$PC2)) + 
# geom_point(aes( colour = metadata_df$RISK_Type_R_Bx, shape = metadata_df$Dx_specific, size =4, alpha = 0.5)) + scale_colour_brewer(palette = 'Paired')

# g = ggplot(data, aes(x=data$PC1, y=data$PC2)) +
#   geom_point(aes( colour = metadata_df$Dx_specific, shape = metadata_df$Source)) # + scale_colour_brewer(palette = 'Paired')
if ( cohort == 'Protect')
    data$PC2 = -1*data$PC2

metadata_df$Cohort[metadata_df$Cohort == 'Protect'] = 'PROTECT'
c[c=='Protect'] = 'PROTECT'

metadata_df$Cohort[metadata_df$Cohort == 'Celiac_Leonard'] = 'PRJNA52875'
c[c=='Celiac_Leonard'] = 'PRJNA52875'

metadata_df$Cohort = factor(metadata_df$Cohort, levels = c)

g = ggplot(data, aes(x=PC1, y=PC2)) +
  geom_point(aes( fill = metadata_df$Dx, shape = metadata_df$Cohort ), size=2) +
  scale_fill_manual(values = c('#0066CC','#FF0033')) + 
  scale_shape_manual(values = c(21,23), 'Cohort') + 
  xlab(sprintf('PC1 (%.2gf%%)', imp[1]*100)) + ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) +
  guides(fill=guide_legend(title='Dx', override.aes = list(shape = 21))) + 
  theme_bw() # + scale_colour_brewer(palette = 'Paired')

metadata_df$Dx_cohort = sprintf('%s %s', metadata_df$Dx, metadata_df$Cohort)
metadata_df$Dx_cohort = factor(metadata_df$Dx_cohort, levels = c( sprintf( '%s %s', levels(metadata_df$Dx)[1],c[1] ),
                                                                  sprintf( '%s %s', levels(metadata_df$Dx)[1],c[2] ),
                                                                  sprintf( '%s %s', levels(metadata_df$Dx)[2],c[1] ),
                                                                  sprintf( '%s %s', levels(metadata_df$Dx)[2],c[2] )))
cols = c('#0E3EDA','#42C2FF','red','#FFC300')
metadata_df$Dx_cohort = factor(metadata_df$Dx_cohort, levels = c( sprintf( '%s %s', levels(metadata_df$Dx)[1],c[1] ),
                                                                  sprintf( '%s %s', levels(metadata_df$Dx)[2],c[1] ),
                                                                  sprintf( '%s %s', levels(metadata_df$Dx)[1],c[2] ),
                                                                  sprintf( '%s %s', levels(metadata_df$Dx)[2],c[2] )))
cols = c('#0E3EDA','red','#42C2FF','#FFC300')

g2 = ggplot(data, aes(x=PC1, y=PC2)) +
  geom_point(aes(fill = metadata_df$Dx_cohort ), size=2, shape = 21) +
  scale_fill_manual(values = cols, name ='') + 
  # geom_point(aes(colour = metadata_df$Dx_cohort ), size=2) +
  # scale_colour_manual(values = cols, name ='') + 
  # scale_shape_manual(values = c(21), '') + 
  xlab(sprintf('PC1 (%.2f%%)', imp[1]*100)) + ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) +
  # guides(fill=guide_legend(title='Dx', override.aes = list(shape = 21))) + 
  theme_bw() # + scale_colour_brewer(palette = 'Paired')

library(ggside)
g3 = ggplot(data, aes(x=PC1, y=PC2)) +
  geom_point(aes(fill = metadata_df$Dx_cohort ), size=2, shape = 21, colour = 'gray25') +
  scale_fill_manual(values = cols, name ='') + 
  # geom_point(aes(colour = metadata_df$Dx_cohort ), size=2) +
  scale_colour_manual(values = cols, name ='') + 
  scale_shape_manual(values = c(21), '') + 
  geom_xsideboxplot(aes(y = metadata_df$Dx_cohort, fill = metadata_df$Dx_cohort), orientation = "y", show.legend = F, outlier.alpha = 0) +
  scale_xsidey_discrete(labels = NULL) +
  geom_ysideboxplot(aes(x = metadata_df$Dx_cohort, fill = metadata_df$Dx_cohort), orientation = "x", show.legend = F, outlier.alpha = 0) +
  scale_ysidex_discrete(labels = NULL) +
  xlab(sprintf('PC1 (%.2f%%)', imp[1]*100)) + ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) +
  # geom_xsidedensity(aes(col = metadata_df$Dx_cohort)) +
  # geom_ysidedensity(aes(col = metadata_df$Dx_cohort)) +
  # guides(fill=guide_legend(title='Dx', override.aes = list(shape = 21))) + 
  theme_bw() + theme(ggside.panel.scale = 0.25, ggside.axis.ticks = element_blank(), panel.grid = element_blank()) 

out_path = 'plots/PCA_by_DE/'
out_file = sprintf('%s/PCA_%s_and_%s_by_%s_DE_%s.pdf', out_path, c[1], c[2], cohort, type)

pos_pc1_c1 = metadata_df$Cohort == c[1]
g_pc1_c1 = ggplot(data[pos_pc1_c1,]) + 
  geom_boxplot(aes(x=metadata_df$Dx_cohort[pos_pc1_c1], y=PC1, fill =metadata_df$Dx_cohort[pos_pc1_c1] )) + 
  scale_fill_manual(values = cols, name ='', breaks = levels(metadata_df$Dx_cohort)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  ylab('PC1') + xlab('') 
res = add_significant_asterix_to_plot_BH_v2(p = g_pc1_c1, var = droplevels(as.factor(metadata_df$Dx_cohort[pos_pc1_c1])), 
                                            val = data$PC1[pos_pc1_c1], test_type = 'wilcox')
g_pc1_c1 = res[[1]]

pos_pc1_c2 = metadata_df$Cohort == c[2]
g_pc1_c2 = ggplot(data[pos_pc1_c2,]) + 
  geom_boxplot(aes(x=metadata_df$Dx_cohort[pos_pc1_c2], y=PC1, fill =metadata_df$Dx_cohort[pos_pc1_c2] )) + 
  scale_fill_manual(values = cols, name ='', breaks = levels(metadata_df$Dx_cohort)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  ylab('PC1') + xlab('') 
res = add_significant_asterix_to_plot_BH_v2(p = g_pc1_c2, var = droplevels(as.factor(metadata_df$Dx_cohort[pos_pc1_c2])), 
                                            val = data$PC1[pos_pc1_c2], test_type = 'wilcox')
g_pc1_c2 = res[[1]]

pos_pc2_c1 = metadata_df$Cohort == c[1]
g_pc2_c1 = ggplot(data[pos_pc2_c1,]) + 
  geom_boxplot(aes(x=metadata_df$Dx_cohort[pos_pc2_c1], y=PC2, fill =metadata_df$Dx_cohort[pos_pc2_c1] )) + 
  scale_fill_manual(values = cols, name ='', breaks = levels(metadata_df$Dx_cohort)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  ylab('PC2') + xlab('') 
res = add_significant_asterix_to_plot_BH_v2(p = g_pc2_c1, var = droplevels(as.factor(metadata_df$Dx_cohort[pos_pc2_c1])), 
                                            val = data$PC2[pos_pc2_c1], test_type = 'wilcox')
g_pc2_c1 = res[[1]]

pos_pc2_c2 = metadata_df$Cohort == c[2]
g_pc2_c2 = ggplot(data[pos_pc2_c2,]) + 
  geom_boxplot(aes(x=metadata_df$Dx_cohort[pos_pc2_c2], y=PC2, fill =metadata_df$Dx_cohort[pos_pc2_c2] )) + 
  scale_fill_manual(values = cols, name ='', breaks = levels(metadata_df$Dx_cohort)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  ylab('PC2') + xlab('') 
res = add_significant_asterix_to_plot_BH_v2(p = g_pc2_c2, var = droplevels(as.factor(metadata_df$Dx_cohort[pos_pc2_c2])), 
                                            val = data$PC2[pos_pc2_c2], test_type = 'wilcox')
g_pc2_c2 = res[[1]]

# g_pc2 = ggplot(data) + geom_boxplot(aes(x=metadata_df$Dx_cohort, y=PC2, fill =metadata_df$Dx_cohort )) + 
#   scale_fill_manual(values = c('#0E3EDA','#42C2FF','red','#FFC300'), name ='') + theme_bw() + 
#   ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) + xlab('') +theme(axis.text.x = element_text(angle = 45, hjust = 1))
# res = add_significant_asterix_to_plot_BH_v2(p = g_pc2, var = as.factor(metadata_df$Dx_cohort), val = data$PC2, test_type = 'wilcox')
# g_pc2 = res[[1]]

# dir.create(out_path)
# ggsave(filename = out_file, plot = g, device = 'jpg', width = 6,height = 4)
# ggsave(out_file, plot = g, device = 'pdf', width = 5,height = 3)


