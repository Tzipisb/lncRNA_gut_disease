source('/pita/users/tzipi/code/R_figs/machine_learing_funcs.R')
library(ggplot2)
library(AUC)

set.seed(4)

fix_group_size_bias2 = function(data, check_val)
{
  groups = levels( check_val )
  min_group_size = min( summary( check_val ) )
  
  pos = c()
  for (g in groups)
  {
    temp =  which( check_val == g )
    pos = c( pos, temp[1:min_group_size] )
  }
  return(pos)
}

if (!exists('cht'))
  cht = 'Protect'

out_path = 'code_for_paper/rnaSeq/random_forest/'
# auc_file = sprintf('%s/AUC_accuracy_RF_v2.txt',out_path)
# cat("Train cohort\tTest cohort\tTrain or Test\tType\tGenes\tAUC\taccuracy\tYouden_point",file=auc_file,sep="\n")

res_df2 = data.frame()

cohorts = c('Protect','SEEM','SOURCE')
# cht = c('Protect')
# for ( cht in cohorts )
{
  if ( cht == 'Protect')
  {
    cohort = 'Protect'
    c = c('Protect','RISK')
    source = 'Rectal'
    test_source = source
    diagnosis_prm = 'Dx'
    to_test = c('UC','Control')
  }
  
  if ( cht == 'SEEM')
  {
    c = c('SEEM','Celiac_Leonard')
    cohort = 'SEEM'
    diagnosis_prm = 'Dx_specific'
    source = 'Duodenum'
    test_source = source
    to_test = c('Celiac','Control','Celiac_active')
  }
  
  if ( cht == 'SOURCE')
  {
    c= c('SOURCE', 'RISK')
    cohort = 'SOURCE'
    diagnosis_prm = 'Dx_specific'
    source = 'terminal_ileum'
    test_source = c('terminal_ileum','ileal')
    to_test = c('CD','Control','iCD')
  }
  
  DE_flag = F
  # for (DE_flag in c(T,F))
  {
    genes_used = ifelse(DE_flag,'DEgenes',sprintf('%sGenes',c[[1]]) )
    
    type = 'lncRNA'
    for ( type in c('lncRNA','ProteinCoding') )
    {
      if ( DE_flag )
      {
        DE_name = sprintf('%s_%s_%s_%s_vs_%s', cohort, type, source, to_test[1],to_test[2])
        # if ( group == 'Train') { DE_name = sprintf('%s_Train',DE_name)}
        DE_file = sprintf('../%s/res/%s_deseq2_res_filtered.tsv',cohort, DE_name)
        DE_df = read.table(file = DE_file, header = T,sep = '\t')
        genes = DE_df$Gene   
      } else
      {
        # genes_file = sprintf('res/gene_lists/main3_TPMfilterd_%s_list.txt',type)
        # genes_file = sprintf( 'plots/venns/lists/%s_%s_TPM_filter.txt', c[1],  type )
        # g_df = read.table(file = genes_file, header = F,sep = '\t')
        cht2 = gsub('SOURCE','SOURCE_v2',cht)
        TPM_filtered_files = sprintf('../%s/%s_txi_%s_TPM_filtered_TPM.txt',cht2, cht, type)
        g_df = read.table(file = TPM_filtered_files, header = T,sep = '\t')
        genes = g_df$Gene
      }
      
      tpm_file =  sprintf('data/All_cohorts_%s_TPM_v8.txt',type)
      tpm = read.table(file = tpm_file, header = T, row.names = 1)
      names(tpm) = make.names(names(tpm))
      tpm = tpm[as.character(genes),]
      
      
      metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v8.txt' 
      metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F)
      metadata_df = metadata_df[order(metadata_df$SampleID),]
      tpm = tpm[, make.names( metadata_df$SampleID )]
      
      # cht = gsub(pattern = '__.*','',names(tpm))=='RISK' ## ,make sure its the same in both versions = no mistakes
      met_cht =  metadata_df$Cohort
      
      #             #### choose data to look at (usign DE genes)
      # # filter to wanted cohort
      tpm = tpm[,met_cht %in% c]
      metadata_df = metadata_df[met_cht %in% c,]
      
      
      # extra filteres
      pos = metadata_df$Source %in% test_source
      pos2 = metadata_df[[diagnosis_prm]] %in% to_test
      pos = pos & pos2 
      # pos = pos2
      tpm = tpm[,pos]
      metadata_df = metadata_df[pos,]
      
      
      # tpm2 = tpm[rowSums(is.na(tpm)) != ncol(tpm),] # remove NAs lines
      tpm2 = tpm[rowSums(is.na(tpm)) == 0,] # remove any line with >1 NA
      tpm2 = tpm2[rowSums(tpm2) != 0,] # remove any gene with 0 expression
      # input_data = t(tpm2)
      tpm2 = as.data.frame(t(tpm2))
      names(tpm2) = make.names(names(tpm2))
      
      
      metadata_df_train = metadata_df[metadata_df$Cohort == c[[1]],]
      tpm_train = tpm2[row.names(tpm2) %in% make.names(metadata_df_train$SampleID) ,]
      
      pos = fix_group_size_bias2(tpm2, as.factor(metadata_df_train$Dx))
      # pos = 1:length(metadata_df_train$Dx)
      metadata_df_train_sf = metadata_df_train[pos,]
      tpm_train_sf = tpm_train[pos,]
      
      metadata_df_test = metadata_df[metadata_df$Cohort == c[[2]],]
      tpm_test = tpm2[row.names(tpm2) %in% make.names(metadata_df_test$SampleID), ]
      
      library(randomForest)
      
      # rf <- randomForest(as.factor(metadata_df_train$Dx) ~ .,data=tpm_train)
      rf <- randomForest(as.factor(metadata_df_train_sf$Dx) ~ .,data=tpm_train_sf)
      
      # pred = predict(rf, newdata=tpm_test)
      prob = predict(rf, newdata=tpm_test, type = "prob")
      df = data.frame(pred = prob[,1], survived = as.factor(metadata_df_test$Dx) )
      youden = calc_youden(df)
      cutoff = youden
      cm = table(metadata_df_test$Dx, prob[,1]<cutoff)
      
      # imp = as.data.frame(importance(rf,type = 2))
      # imp = data.frame( gene = rownames(imp), MeanDecreaseGini = imp$MeanDecreaseGini)
      # imp = imp[ rev(order(imp$MeanDecreaseGini)) ,]
      # 
      # imp_num = 30
      # varImpPlot(rf,  sort = T, n.var=imp_num, main=sprintf("Top %i - Variable Importance",imp_num))#, class = taxa[[check_val]]) 
      
      cal_roc_res = calculate_roc(df, 1, 2, n=100) 
      plot_roc(cal_roc_res, cutoff, 1, 2)
      res2 = plot_pred_type_distribution(df, cutoff) 
      
      
      df = data.frame(pred = prob[,1], survived = as.factor(metadata_df_test$Dx) )
      # cal_roc_res = calculate_roc(df, 1, 2, n=100) 
      
      if(c[1]== 'Protect')
        c[1] = 'PROTECT'
      
      res_main = roc(rf$votes[,1],factor(1 * (rf$y==levels(rf$y)[1] )))
      res_main_df = data.frame(fpr = res_main$fpr, tpr = res_main$tpr, Cohort = sprintf('%s train',c[1]))
      
      res = roc(prob[,1],factor(1 * (df$survived==levels(df$survived)[1] )))
      res_df = data.frame(fpr = res$fpr, tpr = res$tpr, Cohort = sprintf('%s test',c[2]))
      
      res_df = rbind(res_df, res_main_df)
      p_roc <- ggplot(res_df) + 
        geom_line(aes(fpr,tpr, colour = Cohort)) + xlab("FPR") + ylab("TPR") + theme_bw() + 
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                     color="darkgrey", linetype="dashed")
      
      out_file = sprintf('%s/RF_%s_train_%s_test_%s_%s_ROC.pdf',out_path, c[[1]],c[[2]],type, genes_used)
      # ggsave(out_file, plot = p_roc, device = 'pdf', width = 4,height = 3)
      
      res2 = plot_pred_type_distribution(df, youden) 
      # out_file = sprintf('%s/RF_%s_train_%s_test_%s_%s_score.jpg',out_path, c[[1]],c[[2]],type, genes_used)
      # ggsave(filename = out_file, plot = res2[[1]], device = 'jpg', width = 6,height = 4)
      
      cm = res2[[2]]
      cm = factor(cm, c('TP','TN','FN','FP'))
      p = ggplot(data=df, aes(x=survived, y=pred)) + 
        geom_boxplot() + 
        geom_jitter(aes(color=cm)) +
        geom_hline(yintercept=youden, color="red", alpha=0.6) +
        theme_bw() + 
        scale_color_manual(name = "", values = c('blue','red','#33CCFF','#FF9966')) + 
        xlab('') + ylab('Random Forest Score')
      out_file = sprintf('%s/RF_%s_train_%s_test_%s_%s_score.pdf',out_path, c[[1]],c[[2]],type, genes_used)
      # ggsave(out_file,plot = p, device = 'pdf', width = 3,height =3)
      
      
      library(AUC)
      res = roc(df$pred, factor(1 * (df$survived==levels(df$survived)[1])))
      auc_res = print( auc(res,min = 0, max = 1) )
      TP = sum ( ( df$pred >= cutoff ) & ( df$survived == levels(df$survived)[1] ) )
      TN = sum ( ( df$pred < cutoff ) & ( df$survived == levels(df$survived)[2] ) )
      accuracy = (TP+TN) / length(df$survived)
      
      # cat(sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', c[1], c[2], 'Test',type, genes_used, auc_res, accuracy, youden),
      #     file=auc_file,append=TRUE)
      
      auc_res_main = print( auc(res_main,min = 0, max = 1) )
      df_main = data.frame(pred = rf$votes[,1], survived = as.factor(metadata_df_train_sf$Dx) )
      cutoff_main = calc_youden(df_main)
      TP_main = sum ( ( df_main$pred >= cutoff_main ) & ( df_main$survived == levels(df_main$survived)[1] ) )
      TN_main = sum ( ( df_main$pred < cutoff_main ) & ( df_main$survived == levels(df_main$survived)[2] ) )
      accuracy_main = (TP_main+TN_main) / length(df_main$survived)
      
      # cat(sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', c[1], c[2], 'Train',type, genes_used, auc_res_main, accuracy_main, cutoff_main),
      #     file=auc_file,append=TRUE)
      
      if(c[1]== 'PROTECT')
        c[1] = 'Protect'
      
      # dir.create(out_path)
      # ggsave(filename = out_file, plot = g, device = 'jpg', width = 6,height = 4) 
      res_df$Type = type
      
      res_df2  = rbind(res_df2, res_df)
    }
    c2 = c
    c2[c2 == 'Protect'] = 'PROTECT'
    c2[c2 == 'Celiac_Leonard'] = 'PRJNA52875'
    res_df2$Cohort = gsub(' .*','', res_df2$Cohort)
    res_df2$Cohort = gsub('Celiac_Leonard','PRJNA52875', res_df2$Cohort)
    res_df2$Cohort = factor(res_df2$Cohort, levels = c2)
    p_roc2 <- ggplot(res_df2) +
      geom_line(aes(fpr,tpr, colour = Cohort), size=1.5) + xlab("FPR") + ylab("TPR") + theme_bw() +
      facet_wrap(~Type,) +  
      # scale_colour_brewer(palette = 'Set2') + 
      scale_colour_manual(values = c('orange','#42d4f4'), name = element_blank())  + 
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                   color="darkgrey", linetype="dashed")
  }
}






