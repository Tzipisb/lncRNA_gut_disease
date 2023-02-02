# library(ade4)
# library(rgl)
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
# library("plot3D")
# library("scatterplot3d") 

met_cols = c('SampleID','dysbiosis_index','faith_pd',
                  'BarcodeSequence', 'LinkerPrimerSequence', 'Barcode',
                  'Database', 'Cohort','number','Cage','Genotype',
                  'Group','Crossbreeding_origin_.Mother.s_earing_.',
                  'Date_of_parturition','Age_weeks','Sample','Time',
                  'sample_ID','DSS.status','Day',
                  'remark','pn_ID','Age')

met_cols = c('SampleID','DSS.status','Genotype','Geno2','Day', 'Cage')
# col_parm = 'Geno2'
col_parm = 'DSS_status'
col_parm = 'DSS_status'
# cols = 1:21

weight = 'unweighted'

name = 'DSS_Mice'
if ( name == 'DSS_Mice')
{
  path = '../res_d1_noOut2/'
  path = '../res_all_noOut2'
  # path = '../res_WT/'
  # path = '../res_HET_noOut2/'
  # path = '../res_DH/'
  taxa_file = sprintf('%s/biom/%s_deblur_map_L7.1.txt',path, name)
  SV_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, name)
  a_file = sprintf('%s/core-metrics-results/faith_pd_vector/alpha-diversity.tsv',path)
  b_file = sprintf('%s/core-metrics-results/unweighted_unifrac_distance_matrix/distance-matrix.tsv',path)
  pcoa_file = sprintf('%s/core-metrics-results/unweighted_unifrac_pcoa_results/ordination.txt',path)
  out_path =  sprintf('%s/R_res',path)
}
dir.create(out_path)
na_str = c('no_data','_','NA','unknown', 'other','na','No_followup')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str, comment.char = '', quote="\"")



out_path = sprintf('%s/PCoA',out_path)
dir.create(out_path)

res = pcoa_set_data(pcoa_file, taxa_file, met_cols)
pcoa = res[[1]]
dims = res[[2]]

pcoa$DSS_status = ifelse(pcoa$DSS.status == 'Yes','DSS','Water')
pcoa$Day = sprintf('Day %s', pcoa$Day)
pcoa$Day = factor(pcoa$Day, levels = c('Day 1','Day 6','Day 10'))

dim_labs = dims
for ( i in 1:3 )
  dim_labs[[i]] = sprintf("PC%d (%.2f%%)", i, dims[i]*100)

# levels(pcoa[[col_parm]]) = c(levels(pcoa[[col_parm]]), 'America','Israeli_healthy','Israeli_hospitalized','Malawi','Venezuela')
# temp = set_color(pcoa[[col_parm]], cols, data_wanted_levels = c('America','Israeli_healthy','Israeli_hospitalized','Malawi','Venezuela')  )

# tiff(sprintf('plots/97_figs/dEG_PCoA_%s_%s_try.tiff',weight, age), width = 1800, height = 1500, units = "px", res = 300, type='cairo')

# # will not cut edges when saving.
# # par(xpd=TRUE)
# ## ploting
# scatter3D(pcoa$V4, pcoa$V3, pcoa$V2, bg=temp, surface=FALSE, colvar = NA, colkey = F, 
#           pch = 21, alpha=0.6, xlab=dim_labs[3], ylab=dim_labs[2],  zlab=dim_labs[1], cex=0.8, 
# #          theta = tht, phi = pi, main = age)
#           theta = tht, phi = pi)
# # scatter3d(x = pcoa$V2, y = pcoa$V4, z = pcoa$V3, point.col = temp, surface=FALSE)
# # ord = c(3,2,4,5,1)
# # legend("topright",p1, legend = legs[ord], fill = cols[ord], cex=0.8, inset=c(-0.05,0), pt.cex = 1)

# dev.off()

# unqiue colours
library(RColorBrewer)
n <- 40
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=sample(col_vector, n))
col_vector = col_vector[1:21]

# pcoa = pcoa[pcoa$after_feb18 == 'no',]
col_parm = 'Geno2'
all_fig =# ggplot(pcoa, aes_string(x='V2', y='V3' ,colour = col_parm, fill = col_parm)) + 
  ggplot(pcoa, aes_string(x='V2', y='V3' ,colour = col_parm, fill = col_parm, shape = 'DSS_status')) + 
  # ggplot(pcoa, aes_string(x='V2', y='V3' ,colour = col_parm, fill = col_parm, shape = 'Day')) + 
  geom_point(size = 2, alpha= 1) + theme_classic() + 
  theme_classic() + xlab(dim_labs[[1]]) + ylab(dim_labs[[2]]) +
  scale_colour_manual(values = c('#b07fc0','#9f044d','#808080'))
  # scale_colour_manual(values = c('red','blue'))
# scale_color_manual(values=col_vector[c(5:7,9)]) # + scale_fill_manual(values=col_vector)
p = format_fig(fig = all_fig, set_ttl_flag = F, bw_flag = F) 
# # ggsave(sprintf( '%s/%s.pdf', out_path, sprintf('DSS d1 no bad 2 Alpha Diversity %s', check_val) ),plot = p, device = 'pdf', width = 5,height = 4)
# # ggsave(sprintf( '%s/%s.pdf', out_path, sprintf('DSS d1 noOut2 PCoA %s', check_val) ),plot = p, device = 'pdf', width = 6,height = 4)
ggsave(sprintf( '%s/%s.pdf', out_path, sprintf('DSS noOut2 PCoA %s', col_parm) ),plot = p, device = 'pdf', width = 5,height = 3)

# ggsave(sprintf('%s/%s_PCoA_PC12.jpg', out_path, col_parm),plot = all_fig, device = 'jpg', width = 6,height = 4)

cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

all_fig_cage =# ggplot(pcoa, aes_string(x='V2', y='V3' ,colour = col_parm, fill = col_parm)) + 
  ggplot(pcoa, aes_string(x='V2', y='V3' ,colour = 'Cage', shape = 'DSS_status')) + 
  # ggplot(pcoa, aes_string(x='V2', y='V3' ,colour = col_parm, fill = col_parm, shape = 'Day')) + 
  geom_point(size = 2, alpha= 1) +theme_classic() + 
  scale_colour_manual(values = cols) + # guides(colour="none") + 
  xlab(dim_labs[[1]]) + ylab(dim_labs[[2]]) 
p_cage = format_fig(fig = all_fig_cage, set_ttl_flag = F, bw_flag = F) 
ggsave(sprintf( '%s/%s.pdf', out_path, sprintf('DSS noOut2 PCoA %s', 'Cage') ),plot = p_cage, device = 'pdf', width = 5,height = 3)

col_parm = 'DSS.status'
g1 = ggplot(pcoa) + geom_boxplot(aes_string(x=col_parm,y = 'V2'), outlier.alpha = 1, fill = 'gray75') +
  xlab(col_parm) +  ylab(dim_labs[[1]]) +
  theme_bw() + theme(legend.position = "none")

res = add_significant_asterix_to_plot_BH_v2(p = g1, var = as.factor(pcoa[[col_parm]]), val = pcoa$V2, test_type = 'wilcox',print_pvals = T)
g1 = format_fig(fig = res[[1]], set_ttl_flag = F) 
# ggsave(sprintf( '%s/%s.pdf', out_path, sprintf('%s PC1 boxplot', col_parm) ),plot = g1, device = 'pdf', width = 4,height = 3)

g2 = ggplot(pcoa) + geom_boxplot(aes_string(x=col_parm,y = 'V3'), outlier.alpha = 1, fill = 'gray75') +
  xlab(col_parm) +  ylab(dim_labs[[2]]) +
  theme_bw() + theme(legend.position = "none")

res = add_significant_asterix_to_plot_BH_v2(p = g2, var = as.factor(pcoa[[col_parm]]), val = pcoa$V3, test_type = 'wilcox',print_pvals = T)
g2 = format_fig(fig = res[[1]], set_ttl_flag = F) 
# ggsave(sprintf( '%s/%s.pdf', out_path, sprintf('%s PC2 boxplot', col_parm) ),plot = g2, device = 'pdf', width = 4,height = 3)


for ( gn in unique(pcoa$Geno2))
{
  pcoa_f = pcoa[pcoa$Geno2 == gn & !is.na(pcoa$DSS.status),]
  g3 = ggplot(pcoa_f) + geom_boxplot(aes_string(x=col_parm,y = 'V3'), outlier.alpha = 1, fill = 'gray75') +
    xlab(col_parm) +  ylab(dim_labs[[2]]) +
    theme_bw() + theme(legend.position = "none") + ggtitle(gn)
  
  res = add_significant_asterix_to_plot_BH_v2(p = g3, var = as.factor(pcoa_f[[col_parm]]), val = pcoa_f$V3, test_type = 'wilcox',print_pvals = T, show_pval_as_asterix = F, label_size = 4, asterix_scale_var = 2.2)
  g3 = format_fig(fig = res[[1]], set_ttl_flag = F) + scale_y_continuous(limits = c(-0.25, 0.3))
  # ggsave(sprintf( '%s/%s.pdf', out_path, sprintf('%s %s PC2 boxplot', col_parm, gn) ),plot = g3, device = 'pdf', width = 3.5,height = 3)
}



# all_fig = ggplot(pcoa, aes_string(x='V2', y='V3' ,colour = 'age_months', size = 'Group')) + 
#   geom_point(size = 4, alpha= 0.7) +
#   theme_bw() + xlab(dim_labs[[1]]) + ylab(dim_labs[[2]])  # +
# 
