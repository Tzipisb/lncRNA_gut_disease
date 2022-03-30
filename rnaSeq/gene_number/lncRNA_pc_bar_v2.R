library(ggplot2)

# cols = c('#343f56','#387c6d','#e9896a')
# cols = c('#282846','#007580','#fed049')
# cols = c('#00a6ed','#7fb800','#ffb400')

## long
# source('expression_plos/check_expression_main_controls.R')
pers_df2 = pers_df

pers_df2 = pers_df2[pers_df2$pers1 > 0.2,]
pers = pers_df2[,c('pers1','cohort','Type')]

pers$cohort_type = sprintf('%s %s',pers$cohort, pers$Type)
df = as.data.frame(summary(as.factor(pers$cohort_type)))
names(df)[1] = 'n'

df$Cohort = gsub(pattern = ' ProteinCoding','',row.names(df))
df$Cohort = gsub(pattern = ' lncRNA','',df$Cohort)
df$type = gsub(pattern = '.* ','',row.names(df))

# lnc_file = sprintf('res/gene_count_over_pipeline/%s_cohorts_cnt.txt','lncRNA')
# prt_file = sprintf('res/gene_count_over_pipeline/%s_cohorts_cnt.txt','ProteinCoding')
# 
# lnc = read.table(file = lnc_file, header = T, sep = '\t')
# prt = read.table(file = prt_file, header = T, sep = '\t')
# 
# lnc$type = 'lncRNA'
# prt$type = 'ProteinCoding'
# 
# df = rbind(lnc, prt)
# 
# df = df[df$Cohort %in% c('Protect','SOURCE','SEEM'),]
# df =df[,names(df)!= 'TPM_filter_control_genes']

library(reshape2)
df = melt(data = df, id.vars = c('Cohort','type'))
df$variable = as.character(df$variable)
# df$variable[df$variable == 'DE_genes'] = 'DE'
# df$variable[df$variable == 'TPM_filter_genes'] = 'TPM filtered'
# df$variable[df$variable == 'TPM_filter_control_genes'] = 'TPM filtered control'


# df$variable = factor(df$variable,levels = c('TPM filtered','TPM filtered control','DE') )
df$type = factor(df$type,levels = c('ProteinCoding','lncRNA') )

# df$Cohort[df$Cohort == 'Protect'] = 'PROTECT'
# df$Cohort = factor(df$Cohort,levels = rev(c('PROTECT','SOURCE','SEEM')) )

# counts = ggplot(df, aes(x=df$type, y=df$value, group = df$cohort, fill = Cohort, alpha = df$variable)) + 
#   geom_bar(position="dodge2", stat="identity", color = 'black') + 
#   scale_alpha_discrete(range=c(1, 0.2,0.15)) + xlab('') + ylab('') + theme_classic() +
#   labs(alpha = "") + scale_y_continuous(expand=c(0,0)) 

# df2$Cohort = sprintf('%s\n(lncRNA=%s, ProteinCoding=%s)',res2$Cohort, res2$gene_n)


# df2 = df[df$variable == 'TPM filtered control',]
df2 =df

df2$Cohort[df2$Cohort == 'Celiac_Leonard'] = 'PRJNA52875'
df2$Cohort = factor(df2$Cohort, levels = rev(c('PROTECT','SOURCE','SEEM',
                                               'RISK Rectal','RISK Ileal','PRJNA52875')))
df2$type = factor(df2$type, levels = unique(df2$type))

g = ggplot(data=df2, aes(y=Cohort, x=value, fill=type,  label = value)) +
  geom_bar(stat="identity", position="fill")+
  geom_text(size = 4, position = position_fill(vjust = 0.5)) + 
  # geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
  #           color="white", size=3.5)+
  scale_fill_manual(values = c('#DFEDF6','#81B2D4'), name = '') + 
  # scale_fill_brewer(palette ='Paired', name = '') + 
  xlab('Fraction') + 
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw()

out_path = 'plots/gene_number/'
# dir.create(out_path)
# ggsave(filename = sprintf('%s/TPM_all_ctrl_DE_counts.tiff', out_path), plot = counts, device = 'tiff', width = 6,height = 4, compress = 'lzw')
# ggsave(filename = sprintf('%s/TPM_ctrl_lnc_pc_bar_allCohorts.pdf', out_path), plot = g, device = 'pdf', width = 6,height = 3)

