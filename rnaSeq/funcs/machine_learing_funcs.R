source('/pita/users/tzipi/code/R_figs/lay_out.R')

# makes sure all groups in check_val have an identical number of samples in taxa (thus fixing sample size bias)
fix_group_size_bias = function(taxa, check_val)
{
  taxa[[check_val]] = as.factor(taxa[[check_val]])
  taxa = droplevels(taxa)
  groups = levels( taxa[[check_val]] )
  min_group_size = min( summary(taxa[[check_val]]) )
  
  pos = c()
  for (g in groups)
  {
    temp =  which( taxa[[check_val]] == g )
    pos = c(pos, temp[1:min_group_size])
  }
  taxa = taxa[pos,]
  return(taxa)
}

# makes sure all groups in check_val have an identical number of samples in taxa (thus fixing sample size bias)
fix_group_size_bias_rnd = function(taxa, check_val)
{
  taxa[[check_val]] = as.factor(taxa[[check_val]])
  taxa = droplevels(taxa)
  groups = levels( taxa[[check_val]] )
  min_group_size = min( summary(taxa[[check_val]]) )
  
  pos = c()
  for (g in groups)
  {
    temp =  which( taxa[[check_val]] == g )
    rnd_group_pos = sample(x=temp, size=min_group_size, replace = FALSE, prob = NULL)
    
    pos = c(pos, rnd_group_pos)
  }
  taxa = taxa[pos,]
  return(taxa)
}

# log transform the taxa in study taxa (numeric only), add epsolin to avoid 0 problem.
log_transform_taxa = function(taxa)
{
  for ( i in 1:dim(taxa)[1] )
    for ( j in 1:dim(taxa)[2] )
      taxa[i,j] = log2( taxa[i,j] + .Machine$double.eps)
    return(taxa)
}

# remove taxa that appear in only less than cutoff samples
remove_low_appear_taxa = function(taxa, high_cutoff = 5)
{
  for ( i in dim(taxa)[2]:1 )
  {
    no_zero_num = sum(taxa[,i]>0)
    if ( no_zero_num < high_cutoff )
      taxa[,i] = NULL
  }
  return(taxa)
}

plot_pred_type_distribution <- function(df, threshold, rev_level_flag = F) 
{
  # df$pred = df$adjacent_samples_bDiv
  # df$survived = numeric(length = dim(df)[[1]] )
  
  lvls = levels( df$survived )
  if ( rev_level_flag )
  {
    lvls = rev(lvls)
  }
  # df$survived[df$Dx == 'flared'] = 1
  # df$survived[df$Dx == 'never_flared'] = 0
  
  
  v <- rep(NA, nrow(df))
  v <- ifelse(df$pred >= threshold & df$survived == lvls[[1]], "TP", v)
  v <- ifelse(df$pred >= threshold & df$survived == lvls[[2]], "FP", v)
  v <- ifelse(df$pred < threshold & df$survived == lvls[[1]], "FN", v)
  v <- ifelse(df$pred < threshold & df$survived == lvls[[2]], "TN", v)
  
  
  
  # v <- ifelse(df$pred >= threshold & df$survived == 1, "TP", v)
  # v <- ifelse(df$pred >= threshold & df$survived == 0, "FP", v)
  # v <- ifelse(df$pred < threshold & df$survived == 1, "FN", v)
  # v <- ifelse(df$pred < threshold & df$survived == 0, "TN", v)
  
  # df$survived = df$Dx
  
  df$pred_type <- v
  
  p = ggplot(data=df, aes(x=survived, y=pred)) + 
    geom_violin(fill=rgb(1,1,1,alpha=0.6), color=NA) + 
    geom_jitter(aes(color=pred_type), alpha=0.6) +
    geom_hline(yintercept=threshold, color="red", alpha=0.6) +
    scale_color_discrete(name = "type") +
    labs(title=sprintf("Threshold at %.4f", threshold)) + xlab('Group') + ylab('Score')
  
  return(list(p, v))
}

calc_youden <- function(df, rev_level_flag = F) 
{
  # df$pred = df$adjacent_samples_bDiv
  # df$survived = numeric(length = dim(df)[[1]] )
  
  lvls = levels( df$survived )
  if ( rev_level_flag )
  {
    lvls = rev(lvls)
  }
  
  youden_index = vector(mode = 'numeric', length = 100)
  for ( i in 1:100 )
  {
    threshold = i/100
    
    v <- rep(NA, nrow(df))
    v <- ifelse(df$pred >= threshold & df$survived == lvls[[1]], "TP", v)
    v <- ifelse(df$pred >= threshold & df$survived == lvls[[2]], "FP", v)
    v <- ifelse(df$pred < threshold & df$survived == lvls[[1]], "FN", v)
    v <- ifelse(df$pred < threshold & df$survived == lvls[[2]], "TN", v)
    
    sensitivity = sum(v=='TP') / (sum(v=='TP') + sum(v=='FN'))
    specificity = sum(v=='TN') / (sum(v=='TN') + sum(v=='FP')) 
    youden_index[i] = sensitivity + specificity - 1 
  }
  pos = which(youden_index == max(youden_index) ) 

  return( pos[1]/100 )
}

calculate_roc <- function(df, cost_of_fp, cost_of_fn, n=100) 
{
  tpr <- function(df, threshold) {
    sum(df$pred >= threshold & df$survived == levels( df$survived )[[1]]) / sum(df$survived == levels( df$survived )[[1]])
  }
  
  fpr <- function(df, threshold) {
    sum(df$pred >= threshold & df$survived == levels( df$survived )[[2]]) / sum(df$survived == levels( df$survived )[[2]])
  }
  
  cost <- function(df, threshold, cost_of_fp, cost_of_fn) {
    sum(df$pred >= threshold & df$survived == levels( df$survived )[[2]]) * cost_of_fp + 
      sum(df$pred < threshold & df$survived == levels( df$survived )[[1]]) * cost_of_fn
  }
  
  
  roc <- data.frame(threshold = seq(0,1,length.out=n), tpr=NA, fpr=NA)
  roc$tpr <- sapply(roc$threshold, function(th) tpr(df, th))
  roc$fpr <- sapply(roc$threshold, function(th) fpr(df, th))
  roc$cost <- sapply(roc$threshold, function(th) cost(df, th, cost_of_fp, cost_of_fn))
  
  return(roc)
}

plot_roc <- function(roc, threshold, cost_of_fp, cost_of_fn) 
{
  library(gridExtra)
  norm_vec <- function(v) (v - min(v))/diff(range(v))
  
  idx_threshold = which.min(abs(roc$threshold-threshold))
  
  col_ramp <- colorRampPalette(c("green","orange","red","black"))(100)
  col_by_cost <- col_ramp[ceiling(norm_vec(roc$cost)*99)+1]
  p_roc <- ggplot(roc, aes(fpr,tpr)) + 
    geom_line(color=rgb(0,0,1,alpha=0.3)) +
    geom_point(color=col_by_cost, size=4, alpha=0.5) +
    # coord_fixed() +
    geom_line(aes(threshold,threshold), color=rgb(0,0,1,alpha=0.5)) +
    labs(title = sprintf("ROC")) + xlab("FPR") + ylab("TPR") +
    geom_hline(yintercept=roc[idx_threshold,"tpr"], alpha=0.5, linetype="dashed") +
    geom_vline(xintercept=roc[idx_threshold,"fpr"], alpha=0.5, linetype="dashed")
  
  p_cost <- ggplot(roc, aes(threshold, cost)) +
    geom_line(color=rgb(0,0,1,alpha=0.3)) +
    geom_point(color=col_by_cost, size=4, alpha=0.5) +
    labs(title = sprintf("cost function")) +
    geom_vline(xintercept=threshold, alpha=0.5, linetype="dashed")
  
  sub_title <- sprintf("threshold at %.2f - cost of FP = %d, cost of FN = %d", threshold, cost_of_fp, cost_of_fn)
  
  lay_out(list(p_roc , 1, 1),list(p_cost, 1, 2))
}


plot_roc_simple <- function(roc, threshold, cost_of_fp, cost_of_fn) 
{
  library(gridExtra)
  norm_vec <- function(v) (v - min(v))/diff(range(v))
  
  idx_threshold = which.min(abs(roc$threshold-threshold))
  
  col_ramp <- colorRampPalette(c("green","orange","red","black"))(100)
  col_by_cost <- col_ramp[ceiling(norm_vec(roc$cost)*99)+1]
  p_roc <- ggplot(roc, aes(rev(fpr),rev(tpr))) + 
    geom_line(color=rgb(0,0,1)) +
    # geom_point(color=col_by_cost, size=4, alpha=0.5) +
    # coord_fixed() +
    geom_line(aes(threshold,threshold)) +
    labs(title = sprintf("ROC")) + xlab("FPR") + ylab("TPR") # +
    # geom_hline(yintercept=roc[idx_threshold,"tpr"], alpha=0.5, linetype="dashed") +
    # geom_vline(xintercept=roc[idx_threshold,"fpr"], alpha=0.5, linetype="dashed")
  
  # p_cost <- ggplot(roc, aes(threshold, cost)) +
  #   geom_line(color=rgb(0,0,1,alpha=0.3)) +
  #   geom_point(color=col_by_cost, size=4, alpha=0.5) +
  #   labs(title = sprintf("cost function")) +
  #   geom_vline(xintercept=threshold, alpha=0.5, linetype="dashed")
  # 
  # sub_title <- sprintf("threshold at %.2f - cost of FP = %d, cost of FN = %d", threshold, cost_of_fp, cost_of_fn)
  # 
  print(p_roc)
  # lay_out(list(p_roc , 1, 1),list(p_cost, 1, 2))
  return(p_roc)
}

# z score normalize the taxa in study taxa (numeric only)
z_score_norm_taxa = function(taxa)
{
  for ( i in 1:dim(taxa)[2] )
    if ( is.numeric( taxa[, i] ))
      taxa[, i] = ( taxa[, i] - mean( taxa[, i]) ) / sd( taxa[, i] ) 
  return(taxa)
}

z_score_norm = function(val)
{
  return( ( val - mean( val ) ) / sd( val )  )
}


# run random forest on the average value per patients of 'check_vars' variables. 
RF_perPatient_avg = function(taxa, check_vars, b_file = NA, samples_number = 2, last_sample_flag = F)
{
  temp_df = in_patient_var_get(taxa, check_var = check_vars[[1]], samples_number = samples_number, b_file = b_file, last_sample_flag = last_sample_flag)
  dists_df = data.frame(patient = temp_df$patient, Dx = temp_df$Dx)
  
  for ( check_var in check_vars )
  {
    temp_df = in_patient_var_get(taxa, check_var = check_var, samples_number = samples_number, b_file = b_file, last_sample_flag = last_sample_flag)
    dists_df[[check_var]] = temp_df$dists
    # dists_df[[check_var]] = z_score_norm( temp_df$dists )
  }
  if ( 'highest_lewis_score' %in% check_vars )
    dists_df = dists_df[!is.na(dists_df$highest_lewis_score),]
  
  study_data = dists_df[ check_vars ]
  # study_data = z_score_norm_taxa(study_data)
  
  output.forest <- randomForest(dists_df$Dx ~  . , ntree  = 5000, data = study_data)#, cutoff=c(0.7,1-0.7))# , na.action = na.omit)
  
  return(output.forest)
}

# gets the average value per patients of 'check_vars' variable - when you only want the values for a single variable. 
oneVar_perPatient_avg = function(taxa, check_vars, b_file = NA, samples_number = 2, last_sample_flag = F)
{
  temp_df = in_patient_var_get(taxa, check_var = check_vars[[1]], samples_number = samples_number, b_file = b_file, last_sample_flag = last_sample_flag)
  dists_df = data.frame(patient = temp_df$patient, Dx = temp_df$Dx)
  
  for ( check_var in check_vars )
  {
    temp_df = in_patient_var_get(taxa, check_var = check_var, samples_number = samples_number, b_file = b_file, last_sample_flag = last_sample_flag)
    dists_df[[check_var]] = temp_df$dists
    # dists_df[[check_var]] = z_score_norm( temp_df$dists )
  }
  if ( 'highest_lewis_score' %in% check_vars )
    dists_df = dists_df[!is.na(dists_df$highest_lewis_score),]
  
  # study_data = dists_df[ check_vars ]
  # study_data = z_score_norm_taxa(study_data)
  
  # output.forest <- randomForest(dists_df$Dx ~  . , ntree  = 5000, data = study_data)#, cutoff=c(0.7,1-0.7))# , na.action = na.omit)
  
  return(dists_df)
}

# calculate roc
calc_roc = function( pred, real_res )
{
  res = roc( pred, factor( 1 * (real_res==levels(real_res)[1] ) ) )
  return( res )
}

# calculate AUC
calc_auc = function( pred, real_res )
{
  res = calc_roc( pred, real_res )
  auc_res = auc(res,min = 0, max = 1)
  return( auc_res )
}

# ploting the MeanDecreaseGini of random forest result, using wanted parameters number, in ggplot
plot_MeanDecreaseGini = function( output.forest, var_num, sum_name_flag = T )
{
  imp = as.data.frame(importance(output.forest,type = 2))
  imp = data.frame( taxa = rownames(imp), MeanDecreaseGini = imp$MeanDecreaseGini)
  imp = imp[ rev(order(imp$MeanDecreaseGini)) ,]
  
  
  imp = imp[1:var_num, ]
  imp = droplevels(imp)
  if (sum_name_flag) 
  {
    imp$Taxa = taxa_full_name_to_last_level_many(imp$taxa)
  } else
    imp$Taxa = imp$taxa
  imp$Taxa = factor(imp$Taxa, levels = imp$Taxa[ order(imp$MeanDecreaseGini) ] )
  p = ggplot(imp, aes(x=MeanDecreaseGini,y=Taxa)) + geom_point(size=2) + theme_bw() + 
    xlab('Mean decrease\ngini') + ylab('') # + coord_flip()
  return(p)
}