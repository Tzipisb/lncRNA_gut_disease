source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/funcs/figs_funcs.R')
library(patchwork)

diseas = 'UC'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/expression_plos/main_validation_compare.R')
UC_vol_fc = format_fig(fig = g2, set_ttl_flag = F) + theme(legend.position =c(0.29,0.78),panel.grid = element_blank(),legend.box.background = element_rect(colour = "gray")) + 
  xlab('RISK log2(FC)') + ylab('RISK -log10(Q)') + labs(color='PROTECT DE') 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/volcanoDE_UC.tiff', plot =UC_vol_fc, device = 'tiff', compression = 'lzw', width = 3.5,height = 3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/run_DE/Dx_volcano_plot.R')
UC_vol = format_fig(fig = g3, set_ttl_flag = F) + theme(legend.position ='none',panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/volcano_UC.tiff', plot =UC_vol, device = 'tiff', compression = 'lzw', width = 3.5,height = 3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/PCA/PCA_DE_forFig.R')
# UC_PCA = format_fig(fig = g, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
# ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/PCA_UC.tiff', plot =UC_PCA, device = 'tiff', compression = 'lzw', width = 5,height = 3)
UC_PCA2 = format_fig(fig = g2, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/PCA_UC_v2.tiff', plot =UC_PCA2, device = 'tiff', compression = 'lzw', width = 5.5,height = 3)
UC_PCA3 = format_fig(fig = g3, set_ttl_flag = F) + theme(ggside.panel.scale = 0.25, ggside.axis.ticks = element_blank(), panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/PCA_UC_v3.tiff', plot =UC_PCA3, device = 'tiff', compression = 'lzw', width = 6,height = 3.5)

cht = 'Protect'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/random_forest/random_forest_DE_all_v3.R')
UC_rf =  format_fig(fig = p_roc2, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position =c(0.86,0.23)) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/RF_UC.tiff', plot =UC_rf, device = 'tiff', compression = 'lzw', width = 6,height = 3)
UC_rf_f =  UC_rf + facet_wrap(Type~., nrow = 2) + theme(panel.grid = element_blank(), legend.position =c(0.7,0.1)) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/RF_UC_flip.tiff', plot =UC_rf_f, device = 'tiff', compression = 'lzw', width = 3.5,height = 5.5)

# # layout =
# #   '#
# #    A
# #    B
# #    C
# #    D'
# # fig1 =  wrap_plots(A = UC_vol_fc, B = UC_vol, C = UC_PCA, D =UC_rf, design = layout)
# layout ='#A#D
#          #BCD'
# fig1 =  wrap_plots(A = UC_vol_fc + ggtitle('B'), B = UC_vol + ggtitle('C'), 
#                    C = UC_PCA + theme(legend.position = 'none')+ ggtitle('D'), 
#                    D =(UC_rf + facet_wrap(Type~., nrow = 2))+ ggtitle('E'), design = layout)


diseas = 'CD'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/expression_plos/main_validation_compare.R')
CD_vol_fc = format_fig(fig = g2, set_ttl_flag = F) + theme(legend.position =c(0.29,0.78),panel.grid = element_blank(),legend.box.background = element_rect(colour = "gray")) + 
  xlab('RISK log2(FC)') + ylab('RISK -log10(Q)') + labs(color='SOURCE DE') 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/volcanoDE_CD.tiff', plot =CD_vol_fc, device = 'tiff', compression = 'lzw', width = 3.5,height = 3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/run_DE/Dx_volcano_plot.R')
CD_vol = format_fig(fig = g3, set_ttl_flag = F) + theme(legend.position ='none',panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/volcano_CD.tiff', plot =CD_vol, device = 'tiff', compression = 'lzw', width = 3.5,height = 3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/PCA/PCA_DE_forFig.R')
# CD_PCA = format_fig(fig = g, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
# ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/PCA_CD.tiff', plot =CD_PCA, device = 'tiff', compression = 'lzw', width = 5,height = 3)
CD_PCA2 = format_fig(fig = g2, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/PCA_CD_v2.tiff', plot =CD_PCA2, device = 'tiff', compression = 'lzw', width = 5.5,height = 3)
CD_PCA3 = format_fig(fig = g3, set_ttl_flag = F) + theme(ggside.panel.scale = 0.25, ggside.axis.ticks = element_blank(), panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/PCA_CD_v3.tiff', plot =CD_PCA3, device = 'tiff', compression = 'lzw', width = 6,height = 3.5)

cht = 'SOURCE'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/random_forest/random_forest_DE_all_v3.R')
CD_rf =  format_fig(fig = p_roc2, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position =c(0.86,0.23)) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/RF_CD.tiff', plot =CD_rf, device = 'tiff', compression = 'lzw', width = 6,height = 3)
CD_rf_f =  CD_rf + facet_wrap(Type~., nrow = 2) + theme(panel.grid = element_blank(), legend.position =c(0.7,0.1)) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/RF_CD_flip.tiff', plot =CD_rf_f, device = 'tiff', compression = 'lzw', width = 3.5,height = 5.5)


diseas = 'Celiac'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/expression_plos/main_validation_compare.R')
celiac_vol_fc = format_fig(fig = g2, set_ttl_flag = F) + theme(legend.position =c(0.29,0.78),panel.grid = element_blank(),legend.box.background = element_rect(colour = "gray")) + 
  xlab('PRJNA52875 log2(FC)') + ylab('PRJNA52875 -log10(Q)') + labs(color='SEEM DE') 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/volcanoDE_Celiac.tiff', plot =celiac_vol_fc, device = 'tiff', compression = 'lzw', width = 3.5,height = 3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/run_DE/Dx_volcano_plot.R')
Celiac_vol = format_fig(fig = g3, set_ttl_flag = F) + theme(legend.position ='none',panel.grid = element_blank())
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/volcano_celiac.tiff', plot =Celiac_vol, device = 'tiff', compression = 'lzw', width = 3.5,height = 3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/PCA/PCA_DE_forFig.R')
# celiac_PCA = format_fig(fig = g, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
# ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/PCA_celiac.tiff', plot =celiac_PCA, device = 'tiff', compression = 'lzw', width = 5,height = 3)
celiac_PCA2 = format_fig(fig = g2, set_ttl_flag = F) + theme(panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/PCA_celiac_v2.tiff', plot =celiac_PCA2, device = 'tiff', compression = 'lzw', width = 5.5,height = 3)
celiac_PCA3 = format_fig(fig = g3, set_ttl_flag = F) + theme(ggside.panel.scale = 0.25, ggside.axis.ticks = element_blank(), panel.grid = element_blank()) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/PCA_celiac_v3.tiff', plot =celiac_PCA3, device = 'tiff', compression = 'lzw', width = 6,height = 3.5)

cht = 'SEEM'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/random_forest/random_forest_DE_all_v3.R')
Celiac_rf =  format_fig(fig = p_roc2, set_ttl_flag = F) + theme(panel.grid = element_blank(), legend.position =c(0.83,0.22)) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/RF_Celiac.tiff', plot =Celiac_rf, device = 'tiff', compression = 'lzw', width = 6,height = 3)
Celiac_rf_f =  Celiac_rf + facet_wrap(Type~., nrow = 2) + theme(panel.grid = element_blank(), legend.position =c(0.7,0.1)) 
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig1/RF_Celiac_flip.tiff', plot =Celiac_rf_f, device = 'tiff', compression = 'lzw', width = 3.5,height = 5.5)

