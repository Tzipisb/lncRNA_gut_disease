library(reshape)
library(ggplot2)

data_file = 'res/WGCNA/Protect_gene_module_table_lncPrtCd.txt'

df = read.table(file = data_file, header = T, sep = '\t')

modules = unique(df$module_color)
modules = c( 'All', modules )

res = data.frame(module = modules, lncRNA_n = vector(mode = 'numeric', length = length(modules)), 
                                   proteinCoding_n = vector(mode = 'numeric', length = length(modules)),
                                   gene_n = vector(mode = 'numeric', length = length(modules)))

res$lncRNA_n[1] = sum(df$Gene_type == 'lncRNA')
res$proteinCoding_n[1] = sum(df$Gene_type == 'Protein_coding')
res$gene_n[1] = length(df$Gene_type)

for ( i in 2:dim(res)[1] )
{
  pos = df$module_color == res$module[i]
  res$lncRNA_n[i] = sum(df$Gene_type[pos] == 'lncRNA')
  res$proteinCoding_n[i] = sum(df$Gene_type[pos] == 'Protein_coding')
  res$gene_n[i] = length(df$Gene_type[pos])
  
  # temp = data.frame(lncRNA_n = res$lncRNA_n[c(1, i)], proteinCoding_n = res$proteinCoding_n[c(1, i)] )
  # print(res$module[i])
  # print(chisq.test(temp))
}

res$fc_pc_lnc = res$proteinCoding_n / res$lncRNA_n
res2 = reshape::melt(res[1:4], id.vars = c('module','gene_n'))

res2$per = vector(mode = 'numeric',length = length(res2$value))
for ( i in 1:dim(res2)[1] )
{
  res2$per[i] = res2$value[i] / ( res2$value[i] + 
                                    res2$value[ res2$module == res2$module[i] & 
                                                  res2$variable != res2$variable[i] ] )
  res2$per[i] = res2$per[i]*100
}

interesting_modules = c('All','black','green','pink','brown','red','blue','yellow')

res2 = res2[res2$module %in% interesting_modules, ]


res2$Module = sprintf('%s (n=%s)',res2$module, res2$gene_n)
res2$Module = factor(res2$Module, levels = rev(c('All (n=16424)','black (n=237)','green (n=482)','pink (n=167)',
                                             'brown (n=2546)','red (n=239)',
                                             'blue (n=3034)','yellow (n=930)')))
res2$variable = gsub('_n','',res2$variable)

res2$module = factor(res2$module, levels = rev(c('black','green','pink',
                                                 'brown','red',
                                                 'blue','yellow','All')))
# g = ggplot(data=res2, aes(y=Module, x=per, fill=variable)) +
#   geom_bar(stat="identity")+
#   # geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
#   #           color="white", size=3.5)+
#   scale_fill_manual(values = c('gray75','gray35'), name = '') + 
#   # scale_fill_brewer(palette ='Paired', name = '') + 
#   xlab('Percantage') + 
#   scale_x_continuous(expand = c(0, 0)) +
#   theme_bw() 
# g

g = ggplot(data=res2, aes(y=module, x=value, fill=variable,  label = value)) +
  geom_bar(stat="identity", position = 'fill')+
  geom_text(size = 3, position = position_fill(vjust = 0.5), color = 'gray10') + 
  # geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
  #           color="white", size=3.5)+
  # scale_fill_manual(values = c('gray80','gray50'), name = '') + 
  scale_fill_manual(values = c('#DFEDF6','#81B2D4'), name = '') + 
  # scale_fill_brewer(palette ='Paired', name = '') + 
  xlab('Fraction') + 
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw() 
g
# ggsave('plots/module/lncPrtCd_module_coutns.pdf',plot = g, device = 'pdf', width = 6,height =2)

# r = res[2:dim(res)[1],2:3]
# row.names(r) = res$module[ -1 ]
# chisq <- chisq.test(r)
# 
# chisq$observed
# round(chisq$expected,2)
# 
# library(corrplot)
# corrplot(chisq$residuals, is.cor = FALSE)
