

# type = 'ProteinCoding'
type = 'lncRNA'

# diseas = 'UC'
# diseas = 'CD'
# diseas = 'Celiac'

if (!exists('diseas'))
  diseas = 'UC'

if ( diseas == 'UC')
{
  main = 'Protect'
  val = 'RISK_rectalUC'
}

if ( diseas == 'CD')
{
  main = 'SOURCE'
  val = 'RISK_ilealCD'
}

if ( diseas == 'Celiac')
{
  main = 'SEEM'
  val = 'Leonard'
}

in_file = sprintf('res/supps/%s_main_validation_DE_%s_supp.txt',diseas,type)
df = read.table(file = in_file, header = T, sep = '\t')

main_FC = df[[sprintf('%s__log2FoldChange',main)]]
val_FC = df[[sprintf('%s__logFC',val)]]

# main_FC = df[[sprintf('%s__Control_averageTPM',main)]]
# val_FC = df[[sprintf('%s__Control_averageTPM',val)]]

g = ggplot(df, aes(x=main_FC, y=val_FC)) + geom_point() + theme_bw() +
  xlab(sprintf('%s log2(FC)',main)) +  ylab(sprintf('%s log2(FC)',val))
cor = cor.test(main_FC, y=val_FC, method = 'spearman')
out_path = 'plots/main_val_cmp/'
# ggsave(sprintf('%s/tpm_FC_%s_%s_%s_cor%.2f_%e.tiff', out_path, main, val, type, cor$estimate, cor$p.value),plot = g, device = 'tiff', width = 6,height = 4, compression = 'lzw')

# ggplot(df, aes(x=main_FC, y=-1*log10(df$Protect__padj))) + geom_point() + theme_bw() + 
#   xlab('log2(FC)') + ylab('-log10(Q)')

main_DE = ifelse(main_FC>0, 'Up','Down' )
main_DE = factor(main_DE, levels = c('Up','Down'))
val_Q = df[[sprintf('%s__Qvalue',val)]]
g2 = ggplot(df, aes(x=val_FC, y=-1*log10(val_Q))) + 
  geom_point( aes(colour = main_DE)) + 
  theme_bw() + xlab('log2(FC)') + ylab('-log10(Q)')  + 
  scale_colour_manual(values = c('red','blue'))
# ggsave(sprintf('%s/direction_volcano_%s_%s_%s.tiff', out_path, val, main, type ),plot = g2, device = 'tiff', width = 6,height = 4, compression = 'lzw')

# protect_de_file = '../Protect/res_v2/Protect_ProteinCoding_Rectal_UC_vs_Control_deseq2_res.tsv'
# protect_de_file = '../Protect/res_v2/Protect_lncRNA_Rectal_UC_vs_Control_deseq2_res.tsv'

# if (main == 'SOURCE')
#   main_de_file = sprintf('../%s/res_v2/%s_%s_terminal_ileum_CD_vs_Control_deseq2_res.tsv',main, main, type)
# if (main == 'SEEM')
#   main_de_file = sprintf('../%s/res_v2/%s_%s_Doudenum_Celiac_vs_Control_deseq2_res.tsv',main, main, type)
# if (main == 'Protect')
#   main_de_file = sprintf('../%s/res_v2/%s_%s_Rectal_UC_vs_Control_deseq2_res.tsv',main, main, type)
# 
# df = read.table(file = main_de_file, header = T, sep = '\t')
# 
# 
# DE = ifelse(abs(df$log2FoldChange) > log2(1.5) & df$padj < 0.05, 'Up','No')
# DE[DE == 'Up' & df$log2FoldChange<0] = 'Down'
# DE = factor(DE, levels = c('Up','Down','No'))
# g3 = ggplot(df, aes(x=df$log2FoldChange, y=-1*log10(df$padj), colour= DE)) + 
#   geom_point() +
#   scale_colour_manual(values = c('red','blue','gray')) + 
#   theme_bw() + xlab('log2(FC)') + ylab('-log10(Q)')
# # ggsave(sprintf('%s/volcano_%s_%s.tiff', out_path, main, type ),plot = g3, device = 'tiff', width = 6,height = 4, compression = 'lzw')
