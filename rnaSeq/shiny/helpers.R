
# adds significant asterix by p. value to plot, 
# by Benjamini–Hochberg procedure (not Bonferroni correncted here) and wilcox test 
# p = the plot to add on. has to be built as ggplot(data = ?), with aes in geom_boxplot etc.
# var = the variable to divide by. needs to be factor
# val = the values
# p_val_cutoff = cuoff of p.value to show differnce form. default is 0.05
# line_width = width of the line. astetic thing
add_significant_asterix_to_plot_BH = function(p, var, val, p_val_cutoff = 0.05, line_width = 1, print_pvals = F, log_scale_flag=F, test_type = 't_test',label_size = 8, scale_var = 6, asterix_scale_var = 1.3, manual_label = '')
{
  all_vars = levels(var)
  
  # calculate the number of comparisons
  l = length( all_vars )
  
  scale = sd(val, na.rm = T)/scale_var
  # if ( log_scale_flag )
  #   scale = 2^scale
  
  max_pos = max(val, na.rm = T) + scale
  # going over all comparisons between different variable (no repeats)
  # this repeat is to calculate p vals
  
  p_vals = c()
  vars_1 = c()
  vars_2 = c()
  is = c()
  js = c()
  for ( i in 1:l )
  {
    if ( i-1 != 0 ) # because 1:0 is 1 0 and not nothing in R
    {
      for ( j in 1:(i-1) )
      {
        var_1 = all_vars[i]
        var_2 = all_vars[j]
        
        pos1 = var == var_1 & !is.na(var)
        pos2 = var == var_2 & !is.na(var)
        
        # if there is enough data to check
        res = list()
        if ( sum(pos1)>1 & sum(pos2)>1 )
        {
          # calculate p.vlaue
          if ( test_type == 't_test') {
            ##                            big problem!!!!!!!!!!!
            # tryCatch( res = t.test(x = val[pos1], y = val[pos2], na.rm = T),
            #           error = function(e){print(sprintf('t_test error %s',e$message) )}, 
            #           finally = {res$p.value = 1})
            res = t.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else if ( test_type == 'wilcox' ) {
            res = wilcox.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else {print ('test type does not fit to "t_test" or "wilcox')}
        } else {res$p.value = NA}
        
        p_val = res$p.value 
        if ( !is.na(p_val) ) # if a legit p value (enough data to make a comparison, etc.), update it
        {
          p_vals = c(p_vals, p_val)
          vars_1 = c(vars_1, var_1)
          vars_2 = c(vars_2, var_2)
          is = c(is, i)
          js = c(js, j)
        }
      }
    }
  }
  # calculating Benjamini–Hochberg
  BH = p.adjust(p_vals,method="BH") 
  
  # going over all comparisons between differnt variable (no repeats)
  for ( i in 1:l )
  {
    if ( i-1 != 0 ) # because 1:0 is 1 0 and not nothing in R
    {
      for ( j in 1:(i-1) )
      {
        var_1 = all_vars[i]
        var_2 = all_vars[j]
        
        pos = which( is == i & js == j )
        
        if ( length(pos) >1 ) { print('ERROR! check the p value assinged'); return() }
        if ( length(pos) != 0 ) # if there is enough data to check
        {
          if (print_pvals)
          {
            print(sprintf('%s vs %s p.val = %f , BH = %f', var_1, var_2, p_vals[pos], BH[pos] ))
          }
          if ( BH[pos] < p_val_cutoff ) # if p.value is significant enough, add asterix
          {
            # paramters
            astx_pos = max_pos
            start_pos = i
            end_pos = j
            if (manual_label == '')
              lbl = find_pVal_asterix(BH[pos])
            else 
              lbl = manual_label
            
            # building neccecary data frame
            x_pos = c(start_pos, start_pos,end_pos,end_pos)
            y_pos = c(astx_pos, astx_pos+scale*1, astx_pos+scale*1, astx_pos)
            df2 <- data.frame(x_pos = x_pos, y_pos = y_pos )
            
            # adding the asterix to plot
            p = p + geom_path(data = df2, aes(x = x_pos, y = y_pos),size = line_width) + 
              annotate("text", x = mean( c(start_pos, end_pos) ), y = astx_pos+scale*asterix_scale_var, label = lbl, size = label_size) 
            
            # updating maximal position on plot for next time
            max_pos = max_pos + scale*2.5
            if (log_scale_flag)
              max_pos = max_pos + scale*max_pos
          }
        }
      }
    }
  }
  # creating a dataframe with the p value the BH corrected values per varibles 
  pvals_df = data.frame( var_1 = vars_1, var_2 = vars_2, p_val = p_vals, BH = BH )
  return( list(p, pvals_df) )
}

# find fitting asterix to represent the p_value (usually in plots)
find_pVal_asterix = function(p_val)
{
  if (p_val <= 0.001)
  { return('***') }
  if (p_val <= 0.01)
  { return('**') }
  if (p_val <= 0.05)
  { return('*') }
  return('ns')
}

make_disease_boxplot = function(tpm3, lnc)
{
  temp = unique(tpm3$Dx)
  pos = which(temp == 'Control')
  tpm3$Dx = factor(tpm3$Dx, levels = c(temp[temp == 'Control'],temp[temp != 'Control']) ) 
  
  g = ggplot(tpm3) + geom_boxplot(aes(x=Dx, y=tpm3[,lnc]), fill = 'gray70')  + 
    # xlab('Diagnosis') +
    xlab('') +
     ylab(sprintf('%s (TPM) in\n%s',names(tpm3)[lnc], unique(tpm3$Source) )) + 
    theme_bw() + facet_wrap(~Cohort)
  res = add_significant_asterix_to_plot_BH(p = g, var = tpm3$Dx, val = tpm3[,lnc], test_type = 'wilcox')
  if ( ( sum(tpm3[,lnc] > 1 ) / length(tpm3[,lnc]) ) < 0.2 )
    res[[1]] = res[[1]] + ggtitle('Note that expression\nis very low here,\nand should be ignored.', ) +
    theme(plot.title = element_text(color="red", size=14) )
  res[[1]] = res[[1]] + theme(axis.title = element_text(size=18) , axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), strip.text.x = element_text(size = 14))
  return(res[[1]])
  return(g)
}

