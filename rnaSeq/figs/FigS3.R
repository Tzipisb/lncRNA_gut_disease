
type = 'lncRNA'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/PCA/PCA_main3genes_flexible_forFig_TPMfixed.R')
pca_location_lnc = format_fig(fig = g_location, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = sprintf('code_for_paper/rnaSeq_source2/figs/figS3/PCA_allControls_%s_location.pdf', type), plot =pca_location_lnc, device = 'pdf', width = 5,height = 3)
pca_platform_lnc = format_fig(fig = g_Platform, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = sprintf('code_for_paper/rnaSeq_source2/figs/figS3/PCA_allControls_%s_platform.pdf', type), plot =pca_platform_lnc, device = 'pdf', width = 5.5,height = 3)

type = 'ProteinCoding'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/PCA/PCA_main3genes_flexible_forFig_TPMfixed.R')
pca_location_pc = format_fig(fig = g_location, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = sprintf('code_for_paper/rnaSeq_source2/figs/figS3/PCA_allControls_%s_location.pdf', type), plot =pca_location_pc, device = 'pdf', width = 5,height = 3)
pca_platform_pc = format_fig(fig = g_Platform, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = sprintf('code_for_paper/rnaSeq_source2/figs/figS3/PCA_allControls_%s_platform.pdf', type), plot =pca_platform_pc, device = 'pdf', width = 5.5,height = 3)

## long
# source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/expression_plos/check_expression_main_controls.R')
expressed_genes = format_fig(fig = g5, set_ttl_flag = F, x_angle = 45) + theme(panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS3/expressed_genes_control_boxplots.pdf', plot =expressed_genes, device = 'pdf', width = 7,height = 4)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/gene_number/lncRNA_pc_bar_v2.R')
cohort_counts = format_fig(fig = g, set_ttl_flag = F) # + theme(panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS3/Cohorts_counts.pdf', plot =cohort_counts, device = 'pdf', width = 7,height = 3)


disease = 'UC'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/expression_plos/main_validation_avgTPM_cmpr.R')
UC_noise = format_fig(fig = g, lab = sprintf('Spearman\'s rho=%.2f, p<E-10', res$estimate), title_size = 15) + 
  theme(panel.grid = element_blank()) + xlab('PROTECT mean TPM') + ylab('RISK mean TPM')
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS3/UC_noise.pdf', plot =UC_noise, device = 'pdf', width = 5,height = 4)

disease = 'CD'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/expression_plos/main_validation_avgTPM_cmpr.R')
CD_noise = format_fig(fig = g, lab = sprintf('Spearman\'s rho=%.2f, p<E-10', res$estimate), title_size = 15) + 
  theme(panel.grid = element_blank()) + xlab('SOURCE mean TPM') + ylab('RISK mean TPM')
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS3/CD_noise.pdf', plot =CD_noise, device = 'pdf', width = 5,height = 4)

disease = 'Celiac'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/expression_plos/main_validation_avgTPM_cmpr.R')
Celiac_noise = format_fig(fig = g, lab = sprintf('Spearman\'s rho=%.2f, p<E-10', res$estimate), title_size = 15) + 
  theme(panel.grid = element_blank()) + xlab('SEEM mean TPM') + ylab('PRJNA52875 mean TPM')
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS3/Celiac_noise.pdf', plot =Celiac_noise, device = 'pdf', width = 5,height = 4)
