source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/funcs/figs_funcs.R')

## long
# source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/expression_plos/expression_main_cohorts_controls_v2.R')
# control_lncPc_expression = format_fig(fig = g3, set_ttl_flag = F, x_angle = 45) + theme(legend.title=element_blank(), panel.grid = element_blank())
# ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/controls_gene_meanTPM.pdf', plot =control_lncPc_expression, device = 'pdf', width = 7,height = 3)

diseas = 'UC'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/expression_plos/main_validation_compare.R')
UC_vol_fc = format_fig(fig = g2, set_ttl_flag = F) + theme(legend.position =c(0.29,0.78),panel.grid = element_blank(),legend.box.background = element_rect(colour = "gray")) + 
  xlab('RISK log2(FC)') + ylab('RISK -log10(Q)') + labs(color='PROTECT DE') 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/volcanoDE_UC.pdf', plot =UC_vol_fc, device = 'pdf', width = 3.5,height = 3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/run_DE/Dx_volcano_plot.R')
UC_vol = format_fig(fig = g3, set_ttl_flag = F) + theme(legend.position ='none',panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/volcano_UC.pdf', plot =UC_vol, device = 'pdf', width = 3.5,height = 3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/PCA/PCA_DE_forFig.R')
UC_PCA = format_fig(fig = g, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/PCA_UC.pdf', plot =UC_PCA, device = 'pdf', width = 5,height = 3)
UC_PCA2 = format_fig(fig = g2, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/PCA_UC_v2.pdf', plot =UC_PCA2, device = 'pdf', width = 5.5,height = 3)

cht = 'Protect'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/random_forest/random_forest_DE_all_v3.R')
UC_rf =  format_fig(fig = p_roc2, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position =c(0.86,0.23)) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/RF_UC.pdf', plot =UC_rf, device = 'pdf', width = 6,height = 3)



diseas = 'CD'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/expression_plos/main_validation_compare.R')
CD_vol_fc = format_fig(fig = g2, set_ttl_flag = F) + theme(legend.position =c(0.29,0.78),panel.grid = element_blank(),legend.box.background = element_rect(colour = "gray")) + 
  xlab('RISK log2(FC)') + ylab('RISK -log10(Q)') + labs(color='SOURCE DE') 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/volcanoDE_CD.pdf', plot =CD_vol_fc, device = 'pdf', width = 3.5,height = 3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/run_DE/Dx_volcano_plot.R')
CD_vol = format_fig(fig = g3, set_ttl_flag = F) + theme(legend.position ='none',panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/volcano_CD.pdf', plot =CD_vol, device = 'pdf', width = 3.5,height = 3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/PCA/PCA_DE_forFig.R')
CD_PCA = format_fig(fig = g, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/PCA_CD.pdf', plot =CD_PCA, device = 'pdf', width = 5,height = 3)
CD_PCA2 = format_fig(fig = g2, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/PCA_CD_v2.pdf', plot =CD_PCA2, device = 'pdf', width = 5.5,height = 3)

cht = 'SOURCE'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/random_forest/random_forest_DE_all_v3.R')
CD_rf =  format_fig(fig = p_roc2, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position =c(0.86,0.23)) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/RF_CD.pdf', plot =CD_rf, device = 'pdf', width = 6,height = 3)


diseas = 'Celiac'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/expression_plos/main_validation_compare.R')
celiac_vol_fc = format_fig(fig = g2, set_ttl_flag = F) + theme(legend.position =c(0.29,0.78),panel.grid = element_blank(),legend.box.background = element_rect(colour = "gray")) + 
  xlab('PRJNA52875 log2(FC)') + ylab('PRJNA52875 -log10(Q)') + labs(color='SEEM DE') 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/volcanoDE_Celiac.pdf', plot =celiac_vol_fc, device = 'pdf', width = 3.5,height = 3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/run_DE/Dx_volcano_plot.R')
Celiac_vol = format_fig(fig = g3, set_ttl_flag = F) + theme(legend.position ='none',panel.grid = element_blank())
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/volcano_celiac.pdf', plot =Celiac_vol, device = 'pdf', width = 3.5,height = 3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/PCA/PCA_DE_forFig.R')
celiac_PCA = format_fig(fig = g, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/PCA_celiac.pdf', plot =celiac_PCA, device = 'pdf', width = 5,height = 3)
celiac_PCA2 = format_fig(fig = g2, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/PCA_celiac_v2.pdf', plot =celiac_PCA2, device = 'pdf', width = 5.5,height = 3)

cht = 'SEEM'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/random_forest/random_forest_DE_all_v3.R')
Celiac_rf =  format_fig(fig = p_roc2, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position =c(0.83,0.22)) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig3/RF_Celiac.pdf', plot =Celiac_rf, device = 'pdf', width = 6,height = 3)

