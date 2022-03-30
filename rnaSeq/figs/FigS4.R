
diseas = 'UC'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/expression_plos/main_validation_compare.R')
UC_FC_FC = format_fig(fig = g, set_ttl_flag = F) + theme(panel.grid = element_blank()) +
  ylab('RISK log2(FC)') + 
  # annotate("text", x=-1, y=4, size = 5, label= sprintf('Spearman rho=%.2f, p=%.2e', cor$estimate, cor$p.value)) 
  annotate("text", x=-1, y=4, size = 5, label= sprintf('Spearman rho=%.2f, p<E-10', cor$estimate)) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/figS4/UC_FC_FC.pdf', plot =UC_FC_FC, device = 'pdf', width = 6,height = 4)


diseas = 'CD'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/expression_plos/main_validation_compare.R')
CD_FC_FC = format_fig(fig = g, set_ttl_flag = F) + theme(panel.grid = element_blank()) +
  ylab('RISK log2(FC)') + 
  # annotate("text", x=-3, y=4, size = 5, label= sprintf('Spearman rho=%.2f, p=%.2e', cor$estimate, cor$p.value)) 
  annotate("text", x=-3.5, y=4, size = 5, label= sprintf('Spearman rho=%.2f, p<E-10', cor$estimate)) 
ggsave(filename = 'code_for_paper/rnaSeq/figs/figS4/CD_FC_FC.pdf', plot =CD_FC_FC, device = 'pdf', width = 6,height = 4)


diseas = 'Celiac'
source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/code_for_paper/rnaSeq/expression_plos/main_validation_compare.R')
celiac_FC_FC = format_fig(fig = g, set_ttl_flag = F) + theme(panel.grid = element_blank()) +
  ylab('PRJNA52875 log2(FC)') + 
  # annotate("text", x=-1, y=4, size = 5, label= sprintf('Spearman rho=%.2f, p=%.2e', cor$estimate, cor$p.value)) 
  annotate("text", x=-1, y=4, size = 5, label= sprintf('Spearman rho=%.2f, p<E-10', cor$estimate))
ggsave(filename = 'code_for_paper/rnaSeq/figs/figS4/celiac_FC_FC.pdf', plot =celiac_FC_FC, device = 'pdf', width = 6,height = 4)


