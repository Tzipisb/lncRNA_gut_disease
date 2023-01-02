## PCA
# library(car)
# library(devtools)
# library(ggbiplot)

diagnosis_prm = 'Dx'
col_prm = 'Dx'
col_prm = 'Source'

extra_filter_flag = T

# type = 'lncRNA'
# type = 'ProteinCoding'

if (!exists('type'))
  type = 'lncRNA'

# type = 'lncPrtCd'

# cols = c('#343f56','#387c6d','#e9896a')
# cols = c('#282846','#007580','#fed049')

# # summary(as.factor(metadata_df$Dx))
# # c = c('Protect','SEEM','SOURCE','RISK','Celiac_Leonard','Crohn_Yim','Crohn_Yim','blood','caco')
# c = c('Protect','SEEM','SOURCE','RISK','Celiac_Leonard')
# # c = c('Protect','SEEM','SOURCE')
# # to_test = c('Control','UC','CD','Celiac')
# to_test = c('Control')
# test_source = c('Rectal','terminal_ileum','ileal','Duodenum')


# type = 'ProteinCoding'
c = c('Protect','SEEM','SOURCE','RISK','Celiac_Leonard')
to_test = c('Control')
test_source = c('Rectal','terminal_ileum','ileal','Duodenum')

# type = 'lncRNA'
# c = c('Protect','SEEM','SOURCE','RISK','Celiac_Leonard')
# to_test = c('Control')
# test_source = c('Rectal','terminal_ileum','ileal','Duodenum')


# genes_file = sprintf('res/gene_lists/main3_TPMfilterd_%s_list.txt',type)
# g_df = read.table(file = genes_file, header = F,sep = '\t')
# genes = g_df$V1


tpm_file =  sprintf('data/All_cohorts_%s_TPM_v8.txt',type)
tpm = read.table(file = tpm_file, header = T, row.names = 1)
# tpm = tpm[as.character(genes),]
names(tpm) = make.names(names(tpm))

metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v8.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)
metadata_df = metadata_df[order(metadata_df$SampleID),]
tpm = tpm[, make.names( metadata_df$SampleID )]

cht =  metadata_df$Cohort
tpm = tpm[,cht %in% c]
metadata_df = metadata_df[cht %in% c,]

if (extra_filter_flag)
{
  # extra filteres
  pos = metadata_df$Source %in% test_source
  pos2 = metadata_df[[diagnosis_prm]] %in% to_test
  pos = pos & pos2 
  # pos = pos2
  tpm = tpm[,pos]
  metadata_df = metadata_df[pos,]
}


## filter TPM to 1 TPM>0.2%
tpm_genes = row.names(tpm)
wanted_genes_pos2 = vector(mode = 'logical',length = length(tpm_genes))
for ( i in 1:length(tpm_genes) )
{
  per = sum(tpm[i,] > 1) / dim(tpm)[2]
  wanted_genes_pos2[i] = per > 0.2
}
# wanted_genes_pos = wanted_genes_pos1 & wanted_genes_pos2
wanted_genes_pos = wanted_genes_pos2

tpm = tpm[wanted_genes_pos, ]

# tpm2 = tpm[rowSums(is.na(tpm)) != ncol(tpm),] # remove NAs lines
tpm2 = tpm[rowSums(is.na(tpm)) == 0,] # remove any line with >1 NA
tpm2 = tpm2[rowSums(tpm2) != 0,] # remove any gene with 0 expression
input_data = t(tpm2)

df.pca = prcomp(input_data, center = T, scale = T)
imp = summary(df.pca)$importance[2,]

data = as.data.frame(df.pca$x)

new_chr = metadata_df$Cohort
new_chr[metadata_df$Cohort == 'RISK' &  
          metadata_df$Source == 'Rectal' & 
          metadata_df$Dx %in% c('UC','Control')] = 'RISK Rectal'
new_chr[metadata_df$Cohort == 'RISK' &  
          metadata_df$Source == 'ileal' & 
          metadata_df$Dx_specific %in% c('iCD','Control')] = 'RISK Ileal'
metadata_df$new_chr = new_chr

# g = ggplot(data, aes(x=PC1, y=PC2)) +
#   geom_point(aes( colour = metadata_df[[col_prm]], shape = metadata_df$Cohort )) + 
#   guides(colour=guide_legend(title=col_prm))# + scale_colour_brewer(palette = 'Paired')

if (type == 'ProteinCoding')
  data$PC1 = -1*data$PC1

Location = metadata_df$Source
Location[Location %in% c('terminal_ileum','ileal')] = 'Ileal'

metadata_df$new_chr[metadata_df$new_chr=='Protect'] = 'PROTECT'

cols = c('#00a6ed','#7fb800','#ffb400')
cols2 = c('#000099','#336600','#CC6633')

Pipeline = metadata_df$Library_Prep
Pipeline[metadata_df$Cohort %in% c('Protect','SOURCE','SEEM')] = 'Paired_polyA'

shapes = c(21,22,24)
# shapes = c(0,1,2)
g = ggplot(data, aes(x=PC1, y=PC2)) +
  geom_point(aes( fill = metadata_df$new_chr, shape = Location), size=3 ) + 
  theme_bw() + 
  scale_shape_manual(values = shapes) + 
  scale_fill_manual(# values = c(cols, 'black','gray','gray45'),
    values = c(cols, cols2),
    # values = c('blue','springgreen3','red','cyan','orange'), 
    # breaks = c('PROTECT','SOURCE','SEEM','RISK','Celiac_Leonard')) +
    breaks = c('PROTECT','SOURCE','SEEM','RISK Rectal','RISK Ileal','Celiac_Leonard')) +
  guides(fill=guide_legend(title='Cohort', override.aes = list(shape = 21, size=3)))

# g = g + # stat_ellipse(aes(colour = metadata_df$Cohort %in% c('Protect','SOURCE','SEEM')), linetype = 2 ) +
#   stat_ellipse(aes(colour = Pipeline), linetype = 2 ) +
#   # scale_colour_manual(values = c(NA,'#f6511d')) +
#   scale_colour_manual(values = c('gray','gray45','gray55')) +
#   guides(colour = "none")

# g2 = ggplot(data, aes(x=PC1, y=PC2)) +
#   geom_point(aes( fill = metadata_df$Cohort, shape = Location, colour=metadata_df$Dx == 'Control'), size=3 ) + 
#   theme_bw() + 
#   scale_shape_manual(values = c(21,22,24)) + 
#   scale_colour_manual(values = c('red','black')) + 
#   scale_fill_manual(values = c('blue','springgreen3','red','cyan','orange'), 
#                     breaks = c('Protect','SOURCE','SEEM','RISK','Celiac_Leonard')) +
#   guides(fill=guide_legend(title='Cohort', override.aes = list(shape = 21, size=3)))

out_path = 'plots/PCA_fitTPMfilter/'
if ( extra_filter_flag )
{
  out_file = sprintf('%s/PCA_cohorts_%s_dx_%s_%s_Cohort.pdf', out_path, paste(c,collapse = '_'),paste(to_test,collapse = '_'), type)
} else
  out_file = sprintf('%s/PCA_cohorts_%s_%s.pdf', out_path, paste(c,collapse = '_'), type)
dir.create(out_path)
# ggsave(filename = out_file, plot = g, device = 'pdf', width = 6,height = 4)



g_location = ggplot(data, aes(x=PC1, y=PC2)) +
  geom_point(aes( fill = Location), size=3, shape =21 ) + 
  theme_bw() + 
  xlab(sprintf('PC1 (%.2f%%)', imp[1]*100)) + ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) +
  # scale_shape_manual(values = shapes) + 
  scale_fill_manual(# values = c(cols, 'black','gray','gray45'),
    values = c(cols))
# values = c('blue','springgreen3','red','cyan','orange'), 
# breaks = c('PROTECT','SOURCE','SEEM','RISK','Celiac_Leonard')) +
# breaks = c('PROTECT','SOURCE','SEEM','RISK Rectal','RISK Ileal','Celiac_Leonard')) +
# guides(fill=guide_legend(title='Cohort', override.aes = list(shape = 21, size=3)))
out_file = sprintf('%s/PCA_cohorts_%s_%s_Location.pdf', out_path, paste(c,collapse = '_'), type)
# ggsave(filename = out_file, plot = g, device = 'pdf', width = 6,height = 4)


b = ggplot(data, aes(x=Location, y=PC1)) + geom_boxplot( aes( fill = Location) )
res = add_significant_asterix_to_plot_BH_v2(p = b, var = as.factor(Location), val = data$PC1, test_type = 'wilcox')

b = ggplot(data, aes(x=Location, y=PC2)) + geom_boxplot( aes( fill = Location) )
res = add_significant_asterix_to_plot_BH_v2(p = b, var = as.factor(Location), val = data$PC2, test_type = 'wilcox')

g_Platform = ggplot(data, aes(x=PC1, y=PC2)) +
  geom_point(aes( fill = Pipeline), size=3, shape =21 ) + 
  theme_bw() + 
  xlab(sprintf('PC1 (%.2f%%)', imp[1]*100)) + ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) +
  # scale_shape_manual(values = shapes) + 
  scale_fill_manual(# values = c(cols, 'black','gray','gray45'),
    values = c(cols), name = 'Platform')
# values = c('blue','springgreen3','red','cyan','orange'), 
# breaks = c('PROTECT','SOURCE','SEEM','RISK','Celiac_Leonard')) +
# breaks = c('PROTECT','SOURCE','SEEM','RISK Rectal','RISK Ileal','Celiac_Leonard')) +
# guides(fill=guide_legend(title='Cohort', override.aes = list(shape = 21, size=3)))
out_file = sprintf('%s/PCA_cohorts_%s_%s_pipeline.pdf', out_path, paste(c,collapse = '_'), type)
# ggsave(filename = out_file, plot = g, device = 'pdf', width = 6,height = 4)

Use = ifelse(metadata_df$Cohort %in%c('Protect','SOURCE','SEEM'), yes = 'Main',no = 'Validation' )
g = ggplot(data, aes(x=PC1, y=PC2)) +
  geom_point(aes( fill = Use), size=3, shape =21 ) + 
  theme_bw() # + 
# scale_shape_manual(values = shapes) + 
# scale_fill_manual(# values = c(cols, 'black','gray','gray45'),
#   values = c(cols))
# values = c('blue','springgreen3','red','cyan','orange'), 
# breaks = c('PROTECT','SOURCE','SEEM','RISK','Celiac_Leonard')) +
# breaks = c('PROTECT','SOURCE','SEEM','RISK Rectal','RISK Ileal','Celiac_Leonard')) +
# guides(fill=guide_legend(title='Cohort', override.aes = list(shape = 21, size=3)))
out_file = sprintf('%s/PCA_cohorts_%s_%s_Use.pdf', out_path, paste(c,collapse = '_'), type)
# ggsave(filename = out_file, plot = g, device = 'pdf', width = 6,height = 4)



