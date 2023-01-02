library(ggplot2)


data_file = 'modules/data/lncRNA_modules_fdr_e8_yaeEdited.txt'

df = read.table(data_file, sep = '\t',header = T, quote = '')
df = df[df$Include == 1,]

df$edited_new_short_name[df$ID == 'GSM854276_500'] = 'Myeloid_Cells,_DC,_Small_Intestine'
df$edited_new_short_name[df$ID == 'GSM777037_500'] = 'Stromal_Cells'
df$edited_new_short_name = gsub('DC','Dendritic cells (DC)',df$edited_new_short_name )
df$edited_new_short_name = gsub('GN','Granulocytes (GN)',df$edited_new_short_name )
df$edited_new_short_name = gsub('Myeloid_Cells,_Mo','Myeloid_Cells,_Monocytes (Mo)',df$edited_new_short_name )
df$edited_new_short_name = gsub('MF','Macrophage (MF)',df$edited_new_short_name )

df$Name = df$edited_new_short_name
df = df[order(df$Name),]

df$Name = sprintf('%s  ', df$Name)
df$Name = gsub('_',' ', df$Name)
df$Category = gsub('_',' ', df$Category)

interesting_modules = c('blue','yellow','black','green', 'brown','pink', 'red')

temp = df
# # temp  = df[,c('Name','FDR_BH','module','Category','ID')]
# # temp$Name = gsub(',','_',temp$Name)
# temp = reshape2::dcast(temp, Name + Category + ID ~ module, value.var = 'FDR_BH', )
# # row.names(temp) = temp$Name
# # temp = temp[,-1]

df = df[, !names(df) %in% c('edited_new_short_name','Include')]
df = melt(data = df, id.vars = c('Name','Category','ID'),variable_name = 'module')
names(df)[names(df) == 'value']  = 'FDR_BH'
df$FDR_BH[df$FDR_BH <=0.05] = NA

df = droplevels(df)

# if(all(temp$red==0))
# {
#   pos = which(names(temp)=='red')
#   temp = temp[,-pos]
#   df = df[df$module!='red',]
#   interesting_modules = c('blue','yellow','black','green', 'brown','pink')
# }

# heatmap(as.matrix(temp), scale = 'none')


# ord = hclust( dist(temp, method = "euclidean") )$order
# dff$Name = factor( dff$Name, levels = unique(dff$Name)[ord] )

# ord = order( apply(temp[,c('blue','yellow')], 1, mean )  ) 
# dff$Name = factor( dff$Name, levels = sort(unique(dff$Name))[ord] )

# ord = order( apply(temp[,c('black','green','brown','pink')], 1, mean )  ) 
# dff$Name = factor( dff$Name, levels = sort(unique(dff$Name))[ord] )

ord = order( apply(temp[,c('black','green','brown','pink')], 1, mean, na.rm=T ) - apply(temp[,c('blue','yellow')], 1, mean, na.rm=T ) ) 
df$Name = factor( df$Name, levels = sort(unique(df$Name))[ord] )


df$module = factor(df$module, levels = interesting_modules)

unique_cols = c('#e6194b', '#3cb44b', '#ffe119',
                '#4363d8', '#f58231', '#911eb4', '#46f0f0',
                '#f032e6', '#bcf60c', '#9a6324', '#008080',
                '#e6beff', '#aaffc3', '#800000', '#fabebe',
                '#fffac8', '#808000', '#ffd8b1', '#000075',
                '#808080', '#000000') # '#ffffff' = white

# "blue"   "yellow" "black"  "green"  "brown"  "pink"   "red" 
interesting_modules_cols = c('#4363d8','#F6D860','#212121','#4E9F3D','#B05E27','pink','#e6194B')

# "magenta"
g = ggplot(df) + 
  geom_point(aes(x=module, y = Name, colour = Category), shape=15, size=0.1, alpha=0) +
  geom_point(shape=21, aes(x=module, y = Name, fill = module, size=FDR_BH)) +
  theme_bw() + ylab('') + 
  # facet_grid(Database~., scales='free') + 
  scale_fill_manual(breaks = interesting_modules, values = interesting_modules_cols)+ 
  theme(axis.ticks.y = element_blank()) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  guides(fill="none", color = guide_legend(override.aes = list(size = 5, alpha=1) ) ) +
  labs(size='-log10(FDR)') + xlab('Module') + 
  scale_colour_brewer(palette = 'Set3') #
# scale_colour_manual(values = unique_cols) #

# Build a legend "bar"
leg <- ggplot(df, aes(y = Name, x = 0)) + 
  geom_point(aes(color = Category), shape = 15, size = 3, show.legend = F) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) +
  # scale_colour_manual(values = unique_cols) #
  scale_colour_brewer(palette = 'Set3') #
# annotate hm with the bar
g = g + annotation_custom(ggplotGrob(leg), 
                          xmin = -.1, xmax = .3, 
                          ymin = 0, ymax =43.5) # n + .5 since tiles are 1 in width centered on the value
g


out_path = 'plots/module/'
ggsave(sprintf('%s/module_function_yaelFiltered_v2.pdf', out_path),plot = g, device = 'pdf', width = 7,height = 6)

# write.table(x = temp, file = sprintf('%s/module_function_qvals_fdr_e-%d.txt', out_path,fdr_cut), quote = F, sep = '\t',col.names = T, row.names = F)











# df$FDR_BH = -1*log10(df$FDR_BH)
# 
# temp  = df[,c('Name','FDR_BH','module')]
# temp = reshape2::dcast(temp, Name ~ module, value.var = 'FDR_BH')
# row.names(temp) = temp$Name
# temp = temp[,-1]
# temp[is.na(temp)]=0
# 
# # ord = hclust( dist(temp, method = "euclidean") )$order
# # df$Name = factor( df$Name, levels = unique(df$Name)[ord] )
# 
# # ord = hclust( dist(t(temp), method = "euclidean") )$order
# # df$module = factor( df$module, levels = unique(df$module)[ord] )
# g
# df$module = factor(df$module, levels = c('blue','yellow','black','green', 'brown','pink', 'red'))
# # "magenta"
# ggplot(df, aes(x=module, y = Name, fill = module, size=FDR_BH)) + geom_point(shape=21) +
#   theme_bw() + 
#   # facet_grid(Database~., scales='free') + 
#   scale_fill_manual(breaks = interesting_modules, values = interesting_modules)

