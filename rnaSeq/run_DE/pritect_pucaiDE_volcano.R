
tpm_de_file = '../Protect/res_PUCAI/Protect_lncRNA_protect_BASELINE_Pucai_Category_1_vs_3_deseq2_res.tsv'
tpm_de = read.table(tpm_de_file, header = T, sep = '\t')

tpm_de$log2FoldChange = -1*tpm_de$log2FoldChange

status = ifelse(tpm_de$padj <= 0.05 & abs(tpm_de$log2FoldChange) >= log2(1.5), yes = 'val','Not significant')
status[ status == 'val' & tpm_de$log2FoldChange > 0] = gsub('val','Up', status[ status == 'val' & tpm_de$log2FoldChange > 0])
status[ status == 'val' & tpm_de$log2FoldChange < 0] = gsub('val','Down', status[ status == 'val' & tpm_de$log2FoldChange < 0])
tpm_de$status = status

tpm_de$status = factor(tpm_de$status, levels = c('Up','Not significant','Down'))

show_genes = c('GATA6-AS1','CDKN2B-AS1','RP11-396O20.2','CTB-61M7.2','LINC01272','RP11-290L1.3','MIR4435-2HG','DPP10-AS1')
tpm_de$Gene2 = tpm_de$Gene
tpm_de$Gene2[!tpm_de$Gene2 %in% show_genes] = NA

library(ggrepel)
g = ggplot(tpm_de, aes(x=log2FoldChange, y=-1*log10(padj), label=Gene2)) + 
  geom_point( aes(color = status ) ) + 
  geom_text_repel(seed = 4,nudge_y =0.2, min.segment.length = unit(0, 'lines')) +
  theme_bw() + xlab('log2(FC)') + ylab('-log10(Q)')  + 
  scale_colour_manual(values = c('red','gray','blue'), name = '')
g
out_path = '../Protect/res_PUCAI/'
# ggsave(sprintf('%s/Protect_lncRNA_protect_BASELINE_Pucai_Category_3_vs_1_DE_volcano.tiff', out_path),plot = g, device = 'tiff', width = 6,height = 4, compression = 'lzw')

# val_Q = df[[sprintf('%s__Qvalue',val)]]
# g2 = ggplot(df, aes(x=val_FC, y=-1*log10(val_Q))) +
#   geom_point( aes(colour =main_FC>0 )) +
#   theme_bw() + xlab('log2(FC)') + ylab('-log10(Q)')  +
#   scale_colour_manual(values = c('blue','red'))
# ggsave(sprintf('%s/direction_volcano_%s_%s_%s.tiff', out_path, val, main, type ),plot = g2, device = 'tiff', width = 6,height = 4, compression = 'lzw')
