# library(ade4)
# library(rgl)
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
# library("plot3D")
# library("scatterplot3d") 

time2 = 6

met_cols = c('SampleID','dysbiosis_index','faith_pd',
             'BarcodeSequence', 'LinkerPrimerSequence', 'Barcode',
             'Database', 'Cohort','number','Cage','Genotype',
             'Group','Crossbreeding_origin_.Mother.s_earing_.',
             'Date_of_parturition','Age_weeks','Sample','Time',
             'sample_ID','DSS.status','Day',
             'remark','pn_ID','Age')

met_cols = c('SampleID','DSS.status','Genotype','Geno2','Day')
# col_parm = 'Geno2'
col_parm = 'DSS_status'
# cols = 1:21

weight = 'unweighted'

name = 'DSS_Mice'
if ( name == 'DSS_Mice')
{
  path = '../res_all_noOut2'
  taxa_file = sprintf('%s/biom/%s_deblur_map_L7.1.txt',path, name)
  SV_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, name)
  a_file = sprintf('%s/core-metrics-results/faith_pd_vector/alpha-diversity.tsv',path)
  b_file = sprintf('%s/core-metrics-results/unweighted_unifrac_distance_matrix/distance-matrix.tsv',path)
  pcoa_file = sprintf('%s/core-metrics-results/unweighted_unifrac_pcoa_results/ordination.txt',path)
  out_path =  sprintf('%s/R_res',path)
}
# dir.create(out_path)
na_str = c('no_data','_','NA','unknown', 'other','na','No_followup')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str, comment.char = '', quote="\"")

b_df = read.table(b_file, header = T, row.names = 1)

pn_ID = gsub('D[0-9]+_','',taxa$sample_ID)
pn_ID = gsub('D[0-9]+$','',pn_ID)
pn_ID = gsub('_[XYZ][0-9]+$','',pn_ID)

taxa$pn_ID = pn_ID

taxa = taxa[taxa$Day %in% c(1,time2),]
taxa = taxa[taxa$Day == 1 | taxa$DSS.status == 'Yes',]
taxa = taxa[taxa$pn_ID %in% row.names(table(taxa$pn_ID))[table(taxa$pn_ID) == 2], ]

# taxa = taxa[taxa$Geno2!='Het',]

pns = unique(taxa$pn_ID)
df = data.frame(pn_ID = pns, Genotype = 'BAD', distance = -1)

for (i in 1:length(pns))
{
  df$Genotype[i] = taxa$Geno2[taxa$pn_ID == pns[i] & taxa$Day == 1]
  id1 = taxa$SampleID[ taxa$pn_ID == pns[i] & taxa$Day == 1 ]
  id2 = taxa$SampleID[ taxa$pn_ID == pns[i] & taxa$Day != 1 ]
  df$distance[i] = b_df[id1, id2]
}
check_val = 'Genotype'
alpha_var = 'distance'; y_lb = 'Distance'
g = ggplot(df) + geom_boxplot(aes_string(x=check_val,y = alpha_var), outlier.alpha = 1, fill = 'gray75') +
  xlab(check_val) +  ylab(y_lb) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")

r = plot_variable_and_stats(df, check_val, y_val =alpha_var, y_lb = y_lb, sig_ast_flag = F, print_pvals = F, jitter_flag = F)

res = add_significant_asterix_to_plot_BH_v2(p = r, var = as.factor(df[[check_val]]), val = df[[alpha_var]], test_type = 'wilcox',print_pvals = T, show_pval_as_asterix = F, label_size = 4, asterix_scale_var = 2.2,p_val_cutoff =0.4)

p = format_fig(fig = res[[1]], set_ttl_flag = F, x_angle = 30) + xlab('') + scale_fill_manual(values = rep('gray75',10)) + guides(fill = 'none')

# ggsave(sprintf( '%s/%s.pdf', out_path, sprintf('%s_distance_paired_by_mice_between_time_1 and_time_%s_with_DSS', weight, time2) ),plot = p, device = 'pdf', width = 3.5,height = 3.5)
# ggsave(sprintf( '%s/%s.tiff', out_path, sprintf('%s_distance_paired_by_mice_between_time_1 and_time_%s_with_DSS', weight, time2) ),plot = p, device = 'tiff', width = 3.5,height = 3.5, compression='lzw', dpi = 300)
# # ggsave(sprintf( '%s/%s.pdf', out_path, sprintf('%s_distance_paired_by_mice_between_time_1 and_time_%s_with_DSS_noHet.pfd', weight, time2) ),plot = p, device = 'pdf', width = 4,height = 3)


