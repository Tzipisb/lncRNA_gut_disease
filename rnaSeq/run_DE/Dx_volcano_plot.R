
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

type = 'lncRNA'

if (main == 'SOURCE')
  main_de_file = sprintf('../%s/res_v2/%s_%s_terminal_ileum_CD_vs_Control_deseq2_res.tsv',main, main, type)
if (main == 'SEEM')
  main_de_file = sprintf('../%s/res_v2/%s_%s_Doudenum_Celiac_vs_Control_deseq2_res.tsv',main, main, type)
if (main == 'Protect')
  main_de_file = sprintf('../%s/res_v2/%s_%s_Rectal_UC_vs_Control_deseq2_res.tsv',main, main, type)

df = read.table(file = main_de_file, header = T, sep = '\t')

if (main == 'Protect')
{
  main = 'PROTECT'
  # show_genes = c('GATA6-AS1','CDKN2B-AS1','RP11-396O20.2','CTB-61M7.2','LINC01272','RP11-290L1.3','MIR4435-2HG','DPP10-AS1')
  show_genes = c('GATA6-AS1','RP11-396O20.2','MEG3','LINC01272','RP11-363E7.4','LUCAT1')
  df$Gene2 = df$Gene
  df$Gene2[!df$Gene2 %in% show_genes] = NA
  nudge_y = 3; nudge_x = 2
}

if (main == 'SOURCE')
{
  show_genes = c('GATA6-AS1','LINC01272','RP11-16E12.2','CTC-490G23.2','RP11-249C24.11')
  df$Gene2 = df$Gene
  df$Gene2[!df$Gene2 %in% show_genes] = NA
  nudge_y = 0.5; nudge_x =3
}

if (main == 'SEEM')
{
  show_genes = c('LINC01272','RP11-291B21.2','HCP5','LINC01260','RP11-451G4.3','CTD-3080P12.3')
  df$Gene2 = df$Gene
  df$Gene2[!df$Gene2 %in% show_genes] = NA
  nudge_y = 0; nudge_x =0
}

DE = ifelse(abs(df$log2FoldChange) > log2(1.5) & df$padj < 0.05, 'Up','No')
DE[DE == 'Up' & df$log2FoldChange<0] = 'Down'
DE = factor(DE, levels = c('Up','Down','No'))

library(ggrepel)
g3 = ggplot(df, aes(x=log2FoldChange, y=-1*log10(padj), label = Gene2)) + 
  geom_point(aes(colour= DE)) +
  scale_colour_manual(values = c('red','blue','gray')) + 
  theme_bw() +  xlab(sprintf('%s log2(FC)', main)) + ylab(sprintf('%s -log10(Q)', main))
# ggsave(sprintf('%s/volcano_%s_%s.tiff', out_path, main, type ),plot = g3, device = 'tiff', width = 6,height = 4, compression = 'lzw')

# if (main == 'SEEM')
# {
#   g3 = g3+geom_text_repel(seed = 4,size = 4,nudge_y = ifelse(df$Gene =='GATA6-AS1',2,0),
#                           max.overlaps = Inf) 
# } else
#   g3 = g3+geom_text_repel(seed = 4,size = 4,
#                         nudge_x = ifelse(df$log2FoldChange>0, 1, -1)*nudge_x,nudge_y =nudge_y,
#                         max.overlaps = Inf) 
