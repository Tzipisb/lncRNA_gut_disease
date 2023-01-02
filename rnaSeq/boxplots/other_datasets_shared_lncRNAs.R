library(ggplot2)
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

type = 'lncRNA'


tpm_file =  sprintf('data/All_cohorts_%s_TPM_v8.txt',type)
tpm = read.table(file = tpm_file, header = T, row.names = 1)
# tpm[is.na(tpm)] = 0
names(tpm) = make.names(names(tpm))

metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v8.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))


shared_list_file = 'res/supps_v8/main3_validated_genes_cohorts_lncRNA_supp.txt'
shared_list = read.table(file = shared_list_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))
shared_list = shared_list$Gene[shared_list$Score_UC_CD_Celiac == 'fit']


for (gene in shared_list)
{
  # ## blood
  # mt = metadata_df[metadata_df$Cohort == 'blood',]
  # tpm_f = as.data.frame(t(tpm[, mt$SampleID]))
  # 
  # g_blood = ggplot(tpm_f, aes(x=as.factor(mt$Source), y=tpm_f[[gene]])) + 
  #   geom_boxplot(fill = 'gray75', outlier.alpha = 0) + geom_jitter() + 
  #   theme_bw() + ylab(gene) + xlab('') +
  #   theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  # # ggsave(filename = sprintf('code_for_paper/rnaSeq_source2/figs/adde_datasets_boxplots/blood_boxplot_%s.pdf', gene), plot =g_blood, device = 'pdf', width = 4.5,height = 2.5)
  # 
  # mt = metadata_df[metadata_df$Cohort == 'Crohn_Yim',]
  # tpm_f = as.data.frame(t(tpm[, mt$SampleID]))
  # 
  # g_yim = ggplot(tpm_f) + 
  #   geom_boxplot(fill = 'gray75', outlier.alpha = 0, aes(x=as.factor(mt$Dx), y=tpm_f[[gene]])) + 
  #   geom_jitter(aes(x=as.factor(mt$Dx), y=tpm_f[[gene]])) + 
  #   theme_bw() + ylab(gene) + xlab('') + ggtitle('Ileal fibroblast') + 
  #   theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  # res = add_significant_asterix_to_plot_BH_v2(p = g_yim, var = as.factor(mt$Dx), 
  #                                             val = tpm_f[[gene]], test_type = 'wilcox',
  #                                             show_pval_as_asterix = F, print_pvals = T,
  #                                             p_val_cutoff = 1, label_size = 4, asterix_scale_var = 3)
  # g_yim = res[[1]]
  # # ggsave(filename = sprintf('code_for_paper/rnaSeq_source2/figs/adde_datasets_boxplots/yim_boxplot_%s.pdf', gene), plot =g_yim, device = 'pdf', width = 4,height = 3.5)
  # 
  # 
  # mt = metadata_df[metadata_df$Cohort == 'IBD_Howell' & metadata_df$Source == 'terminal_ileum',]
  # tpm_f = as.data.frame(t(tpm[, mt$SampleID]))
  # 
  # g_Howell = ggplot(tpm_f) + 
  #   geom_boxplot(fill = 'gray75', outlier.alpha = 0, aes(x=as.factor(mt$Dx_specific), y=tpm_f[[gene]])) + 
  #   geom_jitter(aes(x=as.factor(mt$Dx_specific), y=tpm_f[[gene]])) + 
  #   theme_bw() + ylab(gene) + xlab('') + ggtitle(unique(mt$Source)) + 
  #   theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  # res = add_significant_asterix_to_plot_BH_v2(p = g_Howell, var = as.factor(mt$Dx_specific), 
  #                                             val = tpm_f[[gene]], test_type = 'wilcox',
  #                                             show_pval_as_asterix = F, print_pvals = T,
  #                                             p_val_cutoff = 0.05, label_size = 4, asterix_scale_var = 3)
  # g_Howell = res[[1]]
  # # ggsave(filename = sprintf('code_for_paper/rnaSeq_source2/figs/adde_datasets_boxplots/howell_TI_boxplot_%s.pdf', gene), plot =g_Howell, device = 'pdf', width = 4,height = 3)
  # 
  # 
  mt = metadata_df[ (metadata_df$Cohort %in% 'blood' | metadata_df$Cohort %in% 'IBD_Howell' ) & metadata_df$Dx %in% 'Control',]
  tpm_f = as.data.frame(t(tpm[, mt$SampleID]))
  
  ## howell, blood order. 
  g_source = ggplot(tpm_f) + 
    geom_boxplot(fill = 'gray75', outlier.alpha = 0, aes(x=as.factor(mt$Source), y=tpm_f[[gene]])) + 
    geom_jitter(aes(x=as.factor(mt$Source), y=tpm_f[[gene]])) + 
    theme_bw() + ylab(gene) + xlab('') +#  ggtitle(unique(mt$Source)) + 
    theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
    facet_wrap(mt$Cohort~., scales = 'free_x' )
  ggsave(filename = sprintf('code_for_paper/rnaSeq_source2/figs/adde_datasets_boxplots/boxplot_source_all3_%s.pdf', gene), plot =g_source, device = 'pdf', width = 4,height = 3)
  
}

