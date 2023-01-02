# library(patchwork)
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/PCA/PCA_forFig_Dx.R')
protect_dx_pca = format_fig(fig = g, set_ttl_flag = F) + theme(legend.position =c(0.84,0.86), legend.background = element_blank(),
                           legend.box.background = element_rect(colour = "gray"), legend.title=element_blank(), panel.grid = element_blank())
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig1/PROTECT_pca_dx.pdf', plot =protect_dx_pca, device = 'pdf', width = 4.5,height = 4)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/PCA/PCA_cor.R')
protect_pca_corrs = format_fig(fig = gc, set_ttl_flag = F) + theme(legend.text = element_text(size=16), panel.grid = element_blank())
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig1/PROTECT_pca_corr.pdf', plot =protect_pca_corrs, device = 'pdf', width = 7,height = 4)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/PCA/PCA_forFig_PUCAI.R')
protect_pucai_pca = format_fig(fig = g, set_ttl_flag = F) + theme(legend.position =c(0.84,0.81), legend.background = element_blank(),
                                                                  legend.box.background = element_rect(colour = "gray",), 
         # ggsave(filename = sprintf('%s/PCA_%s_%s_dx.pdf', out_path, c, type), plot = g, device = 'pdf', width = 6,height = 4)
                                                           legend.text=element_text(size=8), legend.title=element_blank(), panel.grid = element_blank())
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig1/PROTECT_pca_pucai.pdf', plot =protect_pucai_pca, device = 'pdf', width = 4.5,height = 4)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/PCA/PCA_forFig_PC2_boxplot.R')
col_prm = 'protect_BASELINE_Pucai_Category'
col_prm_name = 'Baseline PUCAI'
protect_pucai_pc2 = make_pc2_prm_boxplot(c, col_prm, col_prm_name, type)
protect_pucai_pc2 = format_fig(protect_pucai_pc2, x_angle = 45, set_ttl_flag = F) + theme(legend.position = 'none', panel.grid = element_blank())
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig1/PROTECT_pc2_pucai.pdf', plot =protect_pucai_pc2, device = 'pdf', width = 4.5,height = 4)

col_prm_name = 'Endoscopic Mayo'
col_prm = 'protect_ENDO_MAYO'
protect_mayo_pc2 = make_pc2_prm_boxplot(c, col_prm, col_prm_name, type)
protect_mayo_pc2 = format_fig(protect_mayo_pc2, x_angle = 45, set_ttl_flag = F) + theme(legend.position = 'none', panel.grid = element_blank())
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig1/PROTECT_pc2_mayo.pdf', plot =protect_mayo_pc2, device = 'pdf', width = 4.5,height = 4)

col_prm = 'protect_W4R'
col_prm_name = 'W4 response'
protect_w4r_pc2 = make_pc2_prm_boxplot(c, col_prm, col_prm_name, type)
protect_w4r_pc2 = format_fig(protect_w4r_pc2, x_angle = 45, set_ttl_flag = F) + theme(legend.position = 'none', panel.grid = element_blank())
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig1/PROTECT_pc2_w4r_v2.pdf', plot =protect_w4r_pc2, device = 'pdf', width = 4.5,height = 4)

col_prm = 'protect_3years_colectomy'
col_prm_name = '3 years colectomy'
protect_3yc_pc2 = make_pc2_prm_boxplot(c, col_prm, col_prm_name, type)
protect_3yc_pc2 = format_fig(protect_3yc_pc2, x_angle = 45, set_ttl_flag = F) + theme(legend.position = 'none', panel.grid = element_blank())
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig1/PROTECT_pc2_3yc.pdf', plot =protect_3yc_pc2, device = 'pdf', width = 4.5,height = 4)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/run_DE/pritect_pucaiDE_volcano.R')
protect_pucai_DE = format_fig(g, set_ttl_flag = F) + theme(panel.grid = element_blank())
ggsave(filename = 'code_for_paper/rnaSeq/figs/fig1/PROTECT_pucai_volcano.pdf', plot =protect_pucai_DE, device = 'pdf', width = 6.5,height = 4)

## WGCNA, a long run. save in the script
# source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/WCGNA/WGCNA_protect_v3.R')

# fig1 = plot_spacer() + protect_dx_pca + protect_pucai_pca + protect_pca_corrs + 
#   protect_pucai_pc2 + protect_mayo_pc2 + plot_spacer() + plot_spacer() + 
#   protect_w4r_pc2 + protect_3yc_pc2 + plot_spacer() + plot_spacer() + 
#   # protect_pucai_DE + plot_spacer() + plot_spacer() + plot_spacer() + 
#   plot_layout(ncol = 4)
# 
# out_path = 'code_for_paper/rnaSeq/figs/'
# ggsave(sprintf('%s/fig1.pdf', out_path),plot = fig1, device = 'pdf', width = 18,height = 12)

