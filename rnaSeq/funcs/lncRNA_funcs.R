
# gets genes list, metadata table and deseq2 DE output text file, TPM talbe and metadata to use
# and gets a dataframe with DE FC+qvalue and average TPM values per DE group for each gene in the list.
# NA for each gene that is not in the DE table.
cohort_FC_pqval_avgTPM_table_by_geneList = function(genes, DE, TPM, metadata_df, Cohort, cohort_name, source, diagnosis_prm, to_test )
{
  FC = vector(mode = 'numeric', length = length(genes))
  pval = vector(mode = 'numeric', length = length(genes))
  qval = vector(mode = 'numeric', length = length(genes))
  t1_mean = vector(mode = 'numeric', length = length(genes))
  t2_mean = vector(mode = 'numeric', length = length(genes))
  for ( i in 1:length(genes) )
  {
    gen = genes[i]
    pos = which(DE$Gene == gen)
    
    if ( length(pos) == 0)
    {
      FC[i] = NA
      pval[i] = NA
      qval[i] = NA
      t1_mean[i] = NA
      t2_mean[i] = NA
    } else if ( length(pos) >1 )
    {
      print(sprintf('Problem with %s', gen))
    } else if ( length(pos) == 1 )
    {
      FC[i] = DE$log2FoldChange[pos]
      pval[i] = DE$pvalue[pos]
      qval[i] = DE$padj[pos]
      t1_samps = metadata_df$SampleID[ metadata_df$Cohort == cohort & 
                                         metadata_df$Source == source & 
                                         metadata_df[[diagnosis_prm]] == to_test[1] ]
      t1_mean[i] = mean( as.numeric(TPM[TPM$Gene == gen, make.names(t1_samps)] ) )
      t2_samps = metadata_df$SampleID[ metadata_df$Cohort == cohort & 
                                         metadata_df$Source == source & 
                                         metadata_df[[diagnosis_prm]] == to_test[2] ]
      t2_mean[i] = mean( as.numeric(TPM[TPM$Gene == gen, make.names(t2_samps)] ) )
    }
  }
  df = data.frame(Gene = genes, FC= FC, pval = pval, qval= qval, t1_mean = t1_mean,t2_mean = t2_mean)
  names(df) = c('Gene', sprintf('%s__logFC',cohort_name), 
                sprintf('%s__Pvalue',cohort_name),
                sprintf('%s__Qvalue',cohort_name),
                sprintf('%s__%s_averageTPM',cohort_name,to_test[1]),
                sprintf('%s__%s_averageTPM',cohort_name,to_test[2]))
  return(df)
}


# calculating score of fit = number of NA, tells if akk FC are in the same direction.
# v2 = add FC cutoff
calc_score_v2 = function(df, names, FC_cutoff = log2(1.2) )
{
  score = vector(mode = 'character',length = dim(df)[1] )
  for ( i in 1:dim(df)[1] )
  {
    temp = as.numeric( df[i, names] )
    na_num = sum( is.na( temp ) )
    temp = temp[ !is.na(temp) ]
    if ( all(temp>0) |  all(temp<0) ) 
    {
      if ( na_num == 0 )
      {
        low_FC_num = sum( abs(temp) < FC_cutoff )
        if ( low_FC_num == 0 )
        {
          score[i] = 'fit'
        } else
        {
          score[i] = sprintf('fit_%d_lowFC', low_FC_num)
        }
      } else
      {
        score[i] = sprintf('fit_%d_NAs', na_num)
      } 
    } else
    {
      score[i] = 'no_fit'
    }
  }
  return(score)
}

# calculating score of fit = number of NA, tells if FC are in the same direction.
# v3 = add FC cutoff + Q val
calc_score_v3 = function(df, cohorts, FC_cutoff = log2(1.2),  Qval_cutoff = 0.1 )
{
  FC_names = c( sprintf('%s__log2FoldChange',cohorts), sprintf('%s__logFC',cohorts) )
  FC_names = FC_names[FC_names %in% names(df) ]
  Q_names = c( sprintf('%s__padj',cohorts), sprintf('%s__Qvalue',cohorts) )
  Q_names = Q_names[Q_names %in% names(df) ]
  score = vector(mode = 'character',length = dim(df)[1] )
  for ( i in 1:dim(df)[1] )
  {
    FC_vals = as.numeric( df[i, FC_names] )
    na_num = sum( is.na( FC_vals ) )
    FC_vals = FC_vals[ !is.na(FC_vals) ]
    low_FC_num = sum( abs(FC_vals) <= FC_cutoff )
    
    Q_vals = as.numeric( df[i, Q_names] )
    Q_vals = Q_vals[ !is.na(Q_vals) ]
    high_Q_num = sum( Q_vals >= Qval_cutoff )
    if ( all(FC_vals>0) |  all(FC_vals<0) ) 
    {
      score[i] = 'fit'
    } else
    {
      score[i] = 'no_fit'
    }
    score[i] = sprintf('%s_%d_NAs', score[i], na_num)
    score[i] = sprintf('%s_%d_lowFC', score[i], low_FC_num)
    score[i] = sprintf('%s_%d_highFDR', score[i], high_Q_num)
  }
  score = gsub(pattern = '_0_NAs', replacement = '',x = score, fixed = T)
  score = gsub(pattern = '_0_lowFC', replacement = '',x = score, fixed = T)
  score = gsub(pattern = '_0_highFDR', replacement = '',x = score, fixed = T)
  return(score)
}

# Score for specificity in the main cohorts = what is fit in them and not in the other cohorts.
calc_score_specific = function(df, main_cohorts, other_cohorts , FC_cutoff = log2(1.2),  Qval_cutoff = 0.1 )
{
  main_score = calc_score_v3(df, main_cohorts, FC_cutoff = log2(1.2),  Qval_cutoff = 0.1 )
  all_score = calc_score_v3(df, c(main_cohorts, other_cohorts), FC_cutoff = log2(1.2),  Qval_cutoff = 0.1 )
  
  all_score[main_score != 'fit'] = 'Main_not_fit'
  all_score[grepl(pattern = 'no_fit_',x = all_score)] = 'no_fit'
  return(all_score)
}
