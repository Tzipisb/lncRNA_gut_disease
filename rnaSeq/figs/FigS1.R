
diseas = 'UC'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/expression_plos/main_validation_compare.R')
UC_FC_FC = format_fig(fig = g, set_ttl_flag = F) + theme(panel.grid = element_blank()) +
  ylab('RISK log2(FC)') + 
  # annotate("text", x=-1, y=4, size = 5, label= sprintf('Spearman rho=%.2f, p=%.2e', cor$estimate, cor$p.value)) 
  annotate("text", x=-1, y=4, size = 5, label= sprintf('Spearman rho=%.2f, p<E-10', cor$estimate)) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/UC_FC_FC.pdf', plot =UC_FC_FC, device = 'pdf', width = 6,height = 4)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/PCA/PCA_DE_forFig.R')
# pc1 = format_fig(fig = g_pc1, set_ttl_flag = F, x_angle = 90) + theme(panel.grid = element_blank(), legend.position = 'none') 
# ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_UC_pc1_boxplot.tiff', plot =pc1, device = 'tiff', compression = 'lzw', width = 4,height =4)
pc1_c1 = format_fig(fig = g_pc1_c1, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x =element_blank(), axis.ticks.x =  element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_UC_pc1_c1_boxplot.tiff', plot =pc1_c1, device = 'tiff', compression = 'lzw', width = 2,height =1.7)
pc1_c2 = format_fig(fig = g_pc1_c2, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x =element_blank(), axis.ticks.x =  element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_UC_pc1_c2_boxplot.tiff', plot =pc1_c2, device = 'tiff', compression = 'lzw', width = 2,height =1.7)
pc2_c1 = format_fig(fig = g_pc2_c1, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x =element_blank(), axis.ticks.x =  element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_UC_pc2_c1_boxplot.tiff', plot =pc2_c1, device = 'tiff', compression = 'lzw', width = 2,height =1.7)
pc2_c2 = format_fig(fig = g_pc2_c2, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x =element_blank(), axis.ticks.x =  element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_UC_pc2_c2_boxplot.tiff', plot =pc2_c2, device = 'tiff', compression = 'lzw', width = 2,height =1.7)


diseas = 'CD'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/expression_plos/main_validation_compare.R')
CD_FC_FC = format_fig(fig = g, set_ttl_flag = F) + theme(panel.grid = element_blank()) +
  ylab('RISK log2(FC)') + 
  # annotate("text", x=-3, y=4, size = 5, label= sprintf('Spearman rho=%.2f, p=%.2e', cor$estimate, cor$p.value)) 
  annotate("text", x=-0.5, y=4.5, size = 5, label= sprintf('Spearman rho=%.2f, p<E-10', cor$estimate)) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/CD_FC_FC.pdf', plot =CD_FC_FC, device = 'pdf', width = 6,height = 4)
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/PCA/PCA_DE_forFig.R')
# pc1 = format_fig(fig = g_pc1, set_ttl_flag = F, x_angle = 90) + theme(panel.grid = element_blank(), legend.position = 'none') 
# ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_UC_pc1_boxplot.tiff', plot =pc1, device = 'tiff', compression = 'lzw', width = 4,height =4)
pc1_c1 = format_fig(fig = g_pc1_c1, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x =element_blank(), axis.ticks.x =  element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_CD_pc1_c1_boxplot.tiff', plot =pc1_c1, device = 'tiff', compression = 'lzw', width = 2,height =1.7)
pc1_c2 = format_fig(fig = g_pc1_c2, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x =element_blank(), axis.ticks.x =  element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_CD_pc1_c2_boxplot.tiff', plot =pc1_c2, device = 'tiff', compression = 'lzw', width = 2,height =1.7)
pc2_c1 = format_fig(fig = g_pc2_c1, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x =element_blank(), axis.ticks.x =  element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_CD_pc2_c1_boxplot.tiff', plot =pc2_c1, device = 'tiff', compression = 'lzw', width = 2,height =1.7)
pc2_c2 = format_fig(fig = g_pc2_c2, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x =element_blank(), axis.ticks.x =  element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_CD_pc2_c2_boxplot.tiff', plot =pc2_c2, device = 'tiff', compression = 'lzw', width = 2,height =1.7)


diseas = 'Celiac'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/expression_plos/main_validation_compare.R')
celiac_FC_FC = format_fig(fig = g, set_ttl_flag = F) + theme(panel.grid = element_blank()) +
  ylab('PRJNA52875 log2(FC)') + 
  # annotate("text", x=-1, y=4, size = 5, label= sprintf('Spearman rho=%.2f, p=%.2e', cor$estimate, cor$p.value)) 
  annotate("text", x=-1, y=4, size = 5, label= sprintf('Spearman rho=%.2f, p<E-10', cor$estimate))
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/celiac_FC_FC.pdf', plot =celiac_FC_FC, device = 'pdf', width = 6,height = 4)
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/PCA/PCA_DE_forFig.R')
# pc1 = format_fig(fig = g_pc1, set_ttl_flag = F, x_angle = 90) + theme(panel.grid = element_blank(), legend.position = 'none') 
# ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_UC_pc1_boxplot.tiff', plot =pc1, device = 'tiff', compression = 'lzw', width = 4,height =4)
pc1_c1 = format_fig(fig = g_pc1_c1, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x =element_blank(), axis.ticks.x =  element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_Celiac_pc1_c1_boxplot.tiff', plot =pc1_c1, device = 'tiff', compression = 'lzw', width = 2,height =1.7)
pc1_c2 = format_fig(fig = g_pc1_c2, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x =element_blank(), axis.ticks.x =  element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_Celiac_pc1_c2_boxplot.tiff', plot =pc1_c2, device = 'tiff', compression = 'lzw', width = 2,height =1.7)
pc2_c1 = format_fig(fig = g_pc2_c1, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x =element_blank(), axis.ticks.x =  element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_Celiac_pc2_c1_boxplot.tiff', plot =pc2_c1, device = 'tiff', compression = 'lzw', width = 2,height =1.7)
pc2_c2 = format_fig(fig = g_pc2_c2, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x =element_blank(), axis.ticks.x =  element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/figS1/PCA_Celiac_pc2_c2_boxplot.tiff', plot =pc2_c2, device = 'tiff', compression = 'lzw', width = 2,height =1.7)


