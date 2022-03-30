
lnc_pc_file = 'res/WGCNA/Protect_gene_module_table_lncPrtCd.txt'
lnc_pc = read.table(file = lnc_pc_file, header = T, sep = '\t')

lnc_file = 'res/WGCNA/Protect_gene_module_table_lncRNA.txt'
lnc = read.table(file = lnc_file, header = T, sep = '\t')

modules_lnc_pc = unique(lnc_pc$module_color)
modules_lnc = unique(lnc$module_color)

lnc_pc = lnc_pc[lnc_pc$Gene_type == 'lncRNA', ]

md_pc = c()
md_lnc = c()
n = c()
per = c()
for ( module_lnc in modules_lnc )
{
  for ( module_lnc_pc in modules_lnc_pc )
  {
    md_pc = c(md_pc, module_lnc_pc)
    md_lnc = c(md_lnc, module_lnc)
    
    lnc_list = lnc$Gene[lnc$module_color == module_lnc]
    pos = sum( lnc_list %in% lnc_pc$Gene[ lnc_pc$module_color == module_lnc_pc ]  )
    
    n = c(n, sum(pos))
    per = c(per, sum(pos)/length(lnc_list) )
  }
}

res = data.frame( module_lnc = md_lnc, module_lnc_pc = md_pc, n = n, per = per   )

interesting_pc_modules = c('black','green','pink','brown','red','blue','yellow')
interesting_lnc_modules = c('blue','brown','turquoise','black','red')

res = res[res$module_lnc_pc %in% interesting_pc_modules & res$module_lnc %in% interesting_lnc_modules,]
res$module_lnc = factor(res$module_lnc, levels = interesting_lnc_modules)
res$module_lnc_pc = factor(res$module_lnc_pc, levels = rev(interesting_pc_modules))

library(viridis)
g = ggplot(res, aes(x=module_lnc, y=module_lnc_pc, fill = per*100)) + 
  # geom_tile(colour= 'gray60') +
  geom_tile(colour= 'gray20') +
  geom_text(aes(label=sprintf('%.2f%%\nn=%.0f',per*100,n ), colour =(per > 0.4)), size=3) +
  guides(colour='none') + 
  # scale_fill_gradient(low = 'white',high='#2E5B88', name = 'Percentage') +
  scale_fill_gradient2(low = '#F4F9FC',high='#5D91BD', name = 'Percentage') + 
  # scale_fill_gradient2(low = '#white',high='#73bcff', name = 'Percentage') + 
  # scale_fill_gradient(low = 'lightblue2',high='#1C0C5B', name = 'Percentage') + 
  scale_colour_manual(values = c('gray20','white')) + 
  # scale_fill_viridis(name ='Percantage', option = "G") +
  # scale_fill_gradientn(name ='Percantage',colous = paletteer_c("ggthemes::Blue-Green Sequential", 30) ) +
  # scale_colour_manual(values = c('gray60','black')) + 
  xlab('lncRNA modules') + ylab('lncRNA+ProteinCoding genes modules') +
  scale_x_discrete(expand=c(0,0)) +scale_y_discrete(expand=c(0,0))
g

g2 = ggplot(res, aes(x=module_lnc, y=module_lnc_pc, fill = per*100)) + 
  # geom_tile(colour= 'gray60') +
  geom_tile(colour= 'gray20') +
  geom_text(aes(label=sprintf('%.2f%%\nn=%.0f',per*100,n ), colour =(per > 0.4)), size=3) +
  guides(colour='none') + 
  # scale_fill_gradient(low = 'white',high='#2E5B88', name = 'Percentage') +
  # scale_fill_gradient2(low = '#F4F9FC',high='#5D91BD', name = 'Percentage') + 
  # scale_fill_gradient2(low = '#white',high='#73bcff', name = 'Percentage') + 
  # scale_fill_gradient(low = 'lightblue2',high='#1C0C5B', name = 'Percentage') + 
  # scale_colour_manual(values = c('gray20','white')) + 
  scale_fill_viridis(name ='Percentage') +
  # scale_fill_gradientn(name ='Percentage',colous = paletteer_c("ggthemes::Blue-Green Sequential", 30) ) +
  scale_colour_manual(values = c('gray60','black')) + 
  xlab('lncRNA modules') + ylab('lncRNA+ProteinCoding genes modules') +
  scale_x_discrete(expand=c(0,0)) +scale_y_discrete(expand=c(0,0))
g2
# ggsave('plots/module/lncPrtCd_lnc_modules_shared_lncRNAs.pdf',plot = g, device = 'pdf', width = 6,height =4)

res_table = reshape2::dcast(res, module_lnc ~ module_lnc_pc, value.var = 'per', )
