
## long
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq_source2/expression_plos/expression_main_cohorts_controls_v2.R')
control_lncPc_expression = format_fig(fig = g3, set_ttl_flag = F, x_angle = 45) + theme(legend.title=element_blank(), panel.grid = element_blank())
# ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig2/controls_gene_meanTPM.pdf', plot =control_lncPc_expression, device = 'pdf', width = 7,height = 3)
ggsave(filename = 'code_for_paper/rnaSeq_source2/figs/fig2/controls_gene_meanTPM.tiff', plot =control_lncPc_expression, device = 'tiff',compression = 'lzw', width = 7,height = 3)
