z_score_flag = T
log_norm_flag = F


type = 'lncRNA'
c = 'Protect'
module = 'blue'
top_num = 25

modules = c('blue','brown','turquoise','black','red')

## read data

in_file = sprintf('res/WGCNA/%s_gene_module_table_%s.txt',c, type)

mdf = read.table(in_file, header = T, sep = '\t')

metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v6.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))

TPM_file = sprintf('data/All_cohorts_%s_TPM.txt',type)
tpm = read.table(file = TPM_file, header = T,sep = '\t', stringsAsFactors = F,row.names = 1)

# filter data
metadata_df = metadata_df[order(metadata_df$SampleID),]
tpm = tpm[, make.names(sprintf('%s', metadata_df$SampleID))]

cht =  metadata_df$Cohort
tpm = tpm[,cht %in% c]
metadata_df = metadata_df[cht %in% c,]

colors = c()
genes = c()
# module = 'blue'
for (module in modules)
{
  temp = mdf[mdf$module_color == module, ]
  temp = temp[ rev(order(temp[[ sprintf('MM%s',module) ]])), ]
  gene1 = temp$Gene[1:top_num]
  
  genes =c(genes, gene1)
  colors =c(colors, rep(module, top_num))
}


# module = 'red'
# temp = mdf[mdf$module_color == module, ]
# temp = temp[ rev(order(temp[[ sprintf('MM%s',module) ]])), ]
# gene2 = temp$Gene[1:top_num]

# mdf = mdf[mdf$Gene %in% c(gene1, gene2),]

mdf = mdf[mdf$Gene %in% genes,]

tpm = tpm[mdf$Gene,]

# heatmap(as.matrix(tpm))

if (z_score_flag)
{
  for ( i in 1:dim(tpm)[1])
  {
    t = as.numeric(tpm[i,])
    tpm[i,] = (t - mean(t))/sd(t)
  }
}
if (log_norm_flag)
{
  for ( i in 1:dim(tpm)[1])
  {
    t = as.numeric(tpm[i,])
    t = log2(t+0.1)
    t = t/mean(t)
    tpm[i,] = t
  }
}

old_tpm = tpm
tpm$gene = row.names(tpm)
tpm2 = melt(tpm, id.vars = c('gene'))
temp = data.frame(variable = metadata_df$SampleID, Dx = metadata_df$Dx, pucai = metadata_df$protect_BASELINE_Pucai_Category,
                  colectomy_3y = metadata_df$protect_3years_colectomy, 
                  W52_no_remssion = metadata_df$protect_CSFREE_REMISSION_WK52_KM,
                  W4_no_response = metadata_df$protect_W4R)
tpm3 = merge(x = tpm2, y=temp, by = 'variable', all.x=T)

mycol <- c("navy", "blue2", "deepskyblue1", "lightcyan", "gold", "firebrick2", "red4")

tpm3$gene = factor(tpm3$gene, levels = genes)
tpm3$pucai = as.character(tpm3$pucai)
tpm3$pucai[is.na(tpm3$pucai) ] = 'Control'
tpm3$pucai = factor(tpm3$pucai, levels = c('Control','1','2','3'))

module_gene_df = data.frame(gene = genes, module = colors)
tpm4 = merge(tpm3, module_gene_df , by = 'gene')
library(colorspace)


ord <- hclust( dist(t(old_tpm), method = "euclidean"), method = "ward.D" )$order
tpm4$variable = factor(tpm4$variable, levels = names(tpm)[ord])

tpm4$module[tpm4$module == 'black'] = 'M4'
tpm4$module[tpm4$module == 'blue'] = 'M1'
tpm4$module[tpm4$module == 'brown'] = 'M2'
tpm4$module[tpm4$module == 'red'] = 'M5'
tpm4$module[tpm4$module == 'turquoise'] = 'M3'

g = ggplot(tpm4, aes(x=gene, y=variable, fill = value)) + geom_tile() + 
  scale_fill_gradient2(low = 'blue4',high = 'firebrick3', mid = 'white', midpoint = 0) + 
  # scale_fill_continuous_divergingx(palette = 'Zissou ', mid = 0) + 
  # facet_grid(pucai~module, scale='free', ) + 
  facet_grid(~module, scale='free') + 
  # scale_fill_gradientn(colours = mycol, ) + 
  # scale_y_discrete(limits=names(tpm)[ord], labels = names(tpm)[ord])+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) + 
  xlab('Top 25 genes') + ylab('')
g


# Build a legend "bar"
# leg <- ggplot(tpm3, aes(y = variable, x = 0)) +
#   geom_point(aes(color = Dx), shape = 15, size = 3, show.legend = F) +
#   theme_classic() +
#   theme(axis.title = element_blank(), axis.line = element_blank(),
#         axis.text = element_blank(), axis.ticks = element_blank(),
#         plot.margin = unit(c(0,0,0,0), "cm")) +
#   # scale_colour_manual(values = unique_cols) #
#   scale_colour_manual(values = c('blue','red')) #

leg <- ggplot(tpm4, aes(y = variable, x = 0)) +
  geom_point(aes(color = pucai), shape = 15, size = 3, show.legend = T) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"), legend.position = 'left') +
  # scale_colour_manual(values = unique_cols) #
  scale_colour_manual(values = c('blue','yellow','orange','red')) #

tpm4$colectomy_3y = factor(tpm4$colectomy_3y, levels = unique(tpm4$colectomy_3y))
leg2 <- ggplot(tpm4, aes(y = variable, x = 0)) +
  geom_point(aes(color = colectomy_3y), shape = 15, size = 3, show.legend = T) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"), legend.position = 'left') +
  # scale_colour_manual(values = unique_cols) #
  scale_colour_manual(values = c('peachpuff','peachpuff3')) #

tpm4$W52_no_remssion = factor(tpm4$W52_no_remssion, levels = unique(tpm4$W52_no_remssion))
leg3 <- ggplot(tpm4, aes(y = variable, x = 0)) +
  geom_point(aes(color = W52_no_remssion), shape = 15, size = 3, show.legend = T) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"), legend.position = 'left') +
  # scale_colour_manual(values = unique_cols) #
  scale_colour_manual(values = c('plum1','plum3')) #

tpm4$W4_no_response = factor(tpm4$W4_no_response, levels = unique(tpm4$W4_no_response))
leg4 <- ggplot(tpm4, aes(y = variable, x = 0)) +
  geom_point(aes(color = W4_no_response), shape = 15, size = 3, show.legend = T) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"), legend.position = 'left') +
  # scale_colour_manual(values = unique_cols) #
  scale_colour_manual(values = c('darkseagreen1','darkseagreen3')) #

library(patchwork)
# g2 = leg + g +
#   plot_layout(nrow = 1, widths = c(1, 50))
g2 = leg + leg2 + leg3 +  leg4 + g +
  plot_layout(nrow = 1, widths = c(1, 1, 1, 1, 50))

g3 = (leg +theme(legend.position = 'none') ) + 
  leg2 +theme(legend.position = 'none') + 
  leg3 +theme(legend.position = 'none') + 
  leg4 +theme(legend.position = 'none')+ 
  g +
  plot_layout(nrow = 1, widths = c(1, 1, 1, 1, 50))

out_path = 'plots/WGCNA/'
# out_file = sprintf('%s/top_modules_heatmap_%s_top_%s_%s_pucai.pdf', out_path,module, top_num, type)
out_file = sprintf('%s/top_modules_heatmap_top_%s_%s.pdf', out_path, top_num, type)
out_file2 = sprintf('%s/top_modules_heatmap_top_%s_%s_noLeg.pdf', out_path, top_num, type)
if (z_score_flag)
{
  out_file = gsub('.pdf','_zScore.pdf',out_file)
  out_file2 = gsub('.pdf','_zScore.pdf',out_file2)
}

# ggsave(out_file,plot = g2, device = 'pdf', width = 18,height =9)
# ggsave(out_file2,plot = g3, device = 'pdf', width = 12,height =9)
# 
# write.table(x = t(old_tpm), file = sprintf('%s/top_modules_tpm_top_%s_%s.txt', out_path, top_num, type), row.names = T, sep ='\t', quote = F)
# write.table(x = metadata_df, file = sprintf('%s/top_modules_metadata_top_%s_%s.txt', out_path, top_num, type), row.names = F, sep ='\t', quote = F)
# write.table(x = module_gene_df, file = sprintf('%s/top_modules_gene_module_top_%s_%s.txt', out_path, top_num, type), row.names = F, sep ='\t', quote = F)


# # annotate hm with the bar
# g = g + annotation_custom(ggplotGrob(leg), 
#                           xmin = .2, xmax = .5, 
#                           ymin = 0, ymax =227.5) # n + .5 since tiles are 1 in width centered on the value

# tpm3$gene_dx = sprintf('%s_%s',tpm3$gene,tpm3$Dx)
# 
# tpm4 = aggregate(value~gene_dx, data=tpm3, mean)
# tpm4$gene = gsub('_.*','',tpm4$gene_dx) 
# tpm4$Dx = gsub('.*_','',tpm4$gene_dx)
# tpm4$gene = factor(tpm4$gene, levels = mdf$Gene)
# 
# legend_name = 'Mean TPM'
# if (z_score_flag)
# {
#   legend_name = 'Mean\nZ-score\nTPM' 
# }
# 
# mid_point = ( max(tpm4$value) + min(tpm4$value) )/2
# g2 = ggplot(tpm4, aes(x=gene, y=Dx, fill = value)) + 
#   geom_tile(colour='black', size=0.4) + 
#   # guides(fill=guide_legend(title="Mean TPM"))+
#   scale_fill_gradient2(low = '#3366CC',mid='white',high = '#d53e4f', midpoint = mid_point, name = legend_name) + 
#   # scale_fill_gradient2() + 
#   # scale_fill_fermenter(n.breaks = 11, palette = "PuOr") + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#   xlab('Top 25 genes') + ylab('') +
#   scale_y_discrete(expand=c(0,0))+
#   #define new breaks on x-axis
#   scale_x_discrete(expand=c(0,0)) + coord_flip()
# g2
# 
# out_path = 'plots/WGCNA/'
# out_file = sprintf('%s/module_heatmap_%s_top_%s_%s.pdf', out_path,module, top_num, type)
# if (z_score_flag)
#   out_file = gsub('.pdf','_zScore.pdf',out_file)
# # ggsave(out_file,plot = g2, device = 'pdf', width = 3,height =5)
