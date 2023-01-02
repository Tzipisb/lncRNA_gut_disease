library(tximport)

# library('DESeq2')
# source('/pita/users/tzipi/code/rnaSeq/DESeq2_DE_funcs.r')

# metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v7.txt' 
metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v8.txt' 

# only_control_flag = TRUE
only_control_flag = FALSE

# ony_ileum_flag = TRUE
ony_ileum_flag = FALSE

filter_by_TPM_flag = TRUE
# filter_by_TPM_flag = FALSE

# type = 'lncRNA'
# gene_file = '/pita/users/tzipi/projects/rnaSeq/lncRNA/tables/from_yael/lncRNA_gencode24_forR.txt'

# type = 'ProteinCoding'
# gene_file = '/pita/users/tzipi/projects/rnaSeq/lncRNA/tables/from_yael/19815_Protein_coding_genes_GencodeV24_fixed.txt'

type = 'lncPrtCd'
gene_file = '/pita/users/tzipi/projects/rnaSeq/lncRNA/tables/from_yael/lncRNA_proteinCoding_gencode24.txt'


# cohorts = c('Celiac_Nurit','Celiac_Leonard','Crohn_Yim','IBD_Howell','Protect','RISK','SOURCE','SEEM','caco','blood')
cohorts = c('Celiac_Nurit','Celiac_Leonard','Crohn_Yim','IBD_Howell','Protect','RISK','SOURCE_v2','SEEM','caco','blood')
exmp_files = c('celiac_7002','SRR8774197','SRR5659423','ERR2270963','BR183R','01-001','1-TI2-A002','2540002-9','caco_UN_1','SRR1740034')

# cohort = 'Celiac_Nurit'
# exmp_file = 'celiac_7002'

# cohort = 'Celiac_Leonard'
# exmp_file = 'SRR8774197'

# cohort = 'Crohn_Yim'
# exmp_file = 'SRR5659423'

# cohort = 'IBD_Howell'
# exmp_file = 'ERR2270963'

# cohort = 'Protect'
# exmp_file = 'BR183R'

# cohort = 'RISK'
# exmp_file = '01-001'

cohort = 'SOURCE_v2'
exmp_file = '1-TI2-A002'

# cohort = 'blood'
# exmp_file = 'SRR1740034'

# cohort = 'UC_PRJNA295143'
# exmp_file = 'SRR2313973'
# metadata_file = '../UC_PRJNA295143/data/UC_PRJNA295143_map.txt'

# for ( i in 1:length(cohorts) )
{
  # cohort = cohorts[i]
  # exmp_file = exmp_files[i]

  name = cohort
  input_dir = sprintf('../%s/kallisto',cohort)
  
  cohort2 = gsub('_v2','',cohort)
  
  cohort_specific = cohort2
  if (only_control_flag)
    cohort_specific = sprintf('%s_Control', cohort_specific)

  if (ony_ileum_flag)
    cohort_specific = sprintf('%s_ileal', cohort_specific)
  
  # diagnosis_prm = 'Diagnosis'
  # to_test = c('CELIAC','CONTROL')
  
  metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)
  metadata_df = metadata_df[order(metadata_df$SampleID_noCohort),]
  
  metadata_df = metadata_df[metadata_df$Cohort == cohort2,]
  
  if (only_control_flag)
    metadata_df = metadata_df[metadata_df$Dx == 'Control',]

  if (ony_ileum_flag)
    metadata_df = metadata_df[metadata_df$Source == 'ileal',]

  
  # make a transcript to gene invertion table based on a random kallisto-output file
  example_file = sprintf('%s/%s/abundance.tsv', input_dir, exmp_file) 
  exp = read.table(example_file,sep="\t", header=T)
  trans_names = exp$target_id
  gene_name = as.character(trans_names)
  if ( cohort2 != 'Protect' & cohort2 != 'SEEM')#  & cohort2 != 'SOURCE' ) # those cohort were sent as kallisto and not fastq -> different format
  { 
    for ( i in 1:length(gene_name) )  {   gene_name[i] = strsplit(gene_name[i], '|', fixed = TRUE)[[1]][6] }
  }
  gene_name = gsub(pattern = '_.*',replacement = '',x = gene_name)
  tx2gene = data.frame(TXNAME = trans_names, GENEID = gene_name)
  
  # get all files
  files = metadata_df$SampleID_noCohort
  for( i in 1:length(files)) 
  {
    temp = sprintf('%s/%s/abundance.h5', input_dir, metadata_df$SampleID_noCohort[i])
    # temp = sprintf('%s/*/%s/abundance.h5', input_dir, metadata_df$RNAseq_SampleID[i])
    files[i] = Sys.glob( temp )
    # print(files[i])
    print(Sys.glob( temp ))
  }
  
  # run tximport to summarizr transcript to gene level
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
  
  # filtering tximport results to wanted genes (lncRNA from the table yael send)
  
  gene_df = read.table(gene_file, header = T)
  wanted_genes = gene_df$Gene
  txi_genes = row.names(txi$abundance)
  txi_genes = row.names(txi$abundance)
  wanted_genes_pos = txi_genes %in% wanted_genes
  
  txi$abundance = txi$abundance[wanted_genes_pos, ]
  txi$counts = txi$counts[wanted_genes_pos, ]
  txi$length = txi$length[wanted_genes_pos, ]
  
  # name = sprintf('%s_lncRNA_Filtered', name)
  name = sprintf('%s_%s_Filtered', name, type)
  
  
  if ( filter_by_TPM_flag )
  {
    txi_genes = row.names(txi$abundance)
    wanted_genes_pos = c()
    for ( i in 1:length(txi_genes) )
    {
      per = sum(txi$abundance[i,] > 1) / dim(txi$counts)[2]
      if ( per > 0.2 ) {wanted_genes_pos = c(wanted_genes_pos, i) }
    }
    txi$abundance = txi$abundance[wanted_genes_pos, ]
    txi$counts = txi$counts[wanted_genes_pos, ]
    txi$length = txi$length[wanted_genes_pos, ]
    name = sprintf('%s_TPM_Filtered', name)
  }
  f = sprintf('../%s/%s_txi_%s_TPM_filtered.rds',cohort,cohort_specific, type)
  if ( !filter_by_TPM_flag )
    f = sprintf('../%s/%s_txi_%s.rds',cohort,cohort_specific, type)
  saveRDS(txi, file = f)
  
  TPM_df = txi$abundance
  TPM_df = data.frame(Gene = row.names(TPM_df), TPM_df, row.names=NULL )
  colnames(TPM_df) = c( 'Gene',sprintf('%s__%s',cohort, metadata_df$SampleID_noCohort ) )
  
  # write.table(x = TPM_df, file = sprintf('../%s/%s_txi_lncPrtCd_TPM_filtered_TPM.txt',cohort,cohort), quote = F,sep = '\t',row.names = F )
  # write.table(x = TPM_df, file = sprintf('../%s/%s_txi_lncPrtCd_filtered_TPM.txt',cohort,cohort), quote = F,sep = '\t',row.names = F )
  # write.table(x = TPM_df, file = sprintf('../%s/%s_txi_lncRNA_TPM_filtered_TPM.txt',cohort,cohort), quote = F,sep = '\t',row.names = F )
  # write.table(x = TPM_df, file = sprintf('../%s/%s_txi_ProteinCoding_TPM_filtered_TPM.txt',cohort,cohort), quote = F,sep = '\t',row.names = F )
  
  f = sprintf('../%s/%s_txi_%s_TPM_filtered_TPM.txt',cohort,cohort_specific, type)
  if ( !filter_by_TPM_flag )
    f = sprintf('../%s/%s_txi_%s_TPM.txt',cohort,cohort_specific, type)
  write.table(x = TPM_df, file = f, quote = F,sep = '\t',row.names = F )
}
