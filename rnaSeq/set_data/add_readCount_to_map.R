# install.packages("rjson")
library("rjson")

metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v7.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))

metadata_df$read_count = vector(mode = 'numeric',length = dim(metadata_df)[1] )
for ( i in 1:dim(metadata_df)[1] )
{
  file_path = sprintf('../%s/kallisto/%s/run_info.json', metadata_df$Cohort[i], metadata_df$SampleID_noCohort[i])
  json_data <- fromJSON(paste(readLines(file_path), collapse=""))
  metadata_df$read_count[i] = json_data$n_processed
}

out_file = '../metadata/lncRNA_meta_analysis_megamap_v7_withReadCount.txt' 
write.table(x = metadata_df, file = out_file, quote = F, sep='\t', row.names = F)