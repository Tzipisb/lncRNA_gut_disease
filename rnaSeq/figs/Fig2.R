

## WGCNA, a long run. save in the script
# source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/WCGNA/WGCNA_protect_v3.R')

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/WCGNA/module_function_heatmap_filtered.R')
module_heatmap = g
ggsave('code_for_paper/rnaSeq/figs/fig2/module_function_yaelFiltered_v2.pdf',plot = module_heatmap, device = 'pdf', width = 7,height = 5)


source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/WCGNA/protect_lncRNA_pc_module_cmp.R')
module_per_heatmap = g
ggsave('code_for_paper/rnaSeq/figs/fig2/lncPrtCd_lnc_modules_shared_lncRNAs.pdf',plot = module_per_heatmap, device = 'pdf', width = 5,height =3)
module_per_heatmap2 = g2
ggsave('code_for_paper/rnaSeq/figs/fig2/lncPrtCd_lnc_modules_shared_lncRNAs_v2.pdf',plot = module_per_heatmap2, device = 'pdf', width = 5,height =3)

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/WCGNA/module_lncPC_nums.R')
lncPC_counts = g +ylab('Module')
ggsave('code_for_paper/rnaSeq/figs/fig2/lncPrtCd_module_coutns.pdf',plot = lncPC_counts, device = 'pdf', width = 5.5,height =2.5)


library(patchwork)

layout =
'#A
 #A
 BC'
fig2 =  wrap_plots(A = module_heatmap, B = module_per_heatmap, C = lncPC_counts, design = layout)
