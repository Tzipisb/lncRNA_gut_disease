
library(ggplot2)
library(stringr)

source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

source('/pita/users/tzipi/code/R_figs/maaslin_funcs.R')
source('/pita/users/tzipi/code/R_figs/lay_out.R')

test_type = 'wilcox'


check_val = 'Genotype'
# check_val = 'smoking'

weight = 'unweighted'

name = 'DSS_Mice'
if ( name == 'DSS_Mice')
{
  path = '../res_all/'
  path = '../res_d1_noOut2/'
  # path = '../res_all_noOut2/'
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
# SV = read.table(SV_file,sep="\t", header=TRUE, na.strings = na_str, comment.char = '', quote="\"")

alpha_var = 'faith_pd'; y_lb = 'Faith\'s Phylogenetic Diversity'

taxa = taxa[taxa$Day == 1,]
# taxa = taxa[!taxa$SampleID %in% c('DB18.040','DB18.041'),]

check_val = 'Geno2'
# taxa$Genotype = factor(taxa$Genotype, levels = c('W.T','Het','D.H.'))
taxa$Geno2 = factor(taxa$Geno2, levels = c('WT','Het','KO'))
# alpha_res = plot_variable_and_stats(taxa, check_val, y_val =alpha_var, y_lb = y_lb, sig_ast_flag = T, print_pvals = T, jitter_flag = F, boxplot_p_text_flag = F, show_pval_as_asterix=F)

g = ggplot(taxa) + geom_boxplot(aes_string(x=check_val,y = alpha_var), outlier.alpha = 1, fill = 'gray75') +
  xlab(check_val) +  ylab(y_lb) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")

res = add_significant_asterix_to_plot_BH_v2(p = g, var = as.factor(taxa[[check_val]]), val = taxa[[alpha_var]], test_type = 'wilcox',print_pvals = T, show_pval_as_asterix = F, label_size = 4, asterix_scale_var = 2.2)

p = format_fig(fig = res[[1]], set_ttl_flag = F, bw_flag = F) + xlab('') 
ggsave(sprintf( '%s/%s.pdf', out_path, sprintf('DSS d1 no bad 2 Alpha Diversity %s', check_val) ),plot = p, device = 'pdf', width = 4,height = 4)
# ggsave(sprintf( '%s/%s.pdf', out_path, sprintf('DSS d1 Alpha Diversity %s', check_val) ),plot = p, device = 'pdf', width = 5,height = 4)


# kruskal.test(Genotype ~ faith_pd, data = taxa, )
# 
# kruskal.test(Genotype ~ faith_pd, data = taxa[taxa$Genotype %in% c('W.T','Het'), ])
# kruskal.test(Genotype ~ faith_pd, data = taxa[taxa$Genotype %in% c('D.H.','Het'), ])
