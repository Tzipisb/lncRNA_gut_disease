# lncRNA_gut_disease

This is the code associated with the paper : 
lncRNA atlas of inflammatory diseases along the GI tract highlights regulatory mitochondrial metabolic epithelial functions of GATA6-AS1  

Here you can find the code itself. rnaSeq data is publicly available at the following:


Rectal UC: PROTECT (GSE10914), RISK rectal (GSE117993)

Ileal Crohn's disease: SOURCE (GSE199906), RISK ileal (GSE101794)

Duodenum Celiac: SEEM (GSE159495), leonard et al. (PRJNA5287557)


16S PROTECT data is publicly available at PRJNA436359.





R sessionInfo:

R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] tximport_1.22.0             rjson_0.2.21                DESeq2_1.34.0              
 [4] SummarizedExperiment_1.24.0 Biobase_2.54.0              MatrixGenerics_1.6.0       
 [7] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1         IRanges_2.28.0             
[10] S4Vectors_0.32.4            BiocGenerics_0.40.0         ggrepel_0.9.2              
[13] randomForest_4.7-1.1        AUC_0.3.2                   WGCNA_1.72-1               
[16] fastcluster_1.2.3           dynamicTreeCut_1.63-1       Maaslin2_1.8.0             
[19] vegan_2.6-4                 lattice_0.20-45             permute_0.9-7              
[22] matrixStats_0.63.0          stringr_1.5.0               reshape_0.8.9              
[25] patchwork_1.1.2             ggplot2_3.4.0               shinysky_0.1.3             
[28] shinythemes_1.2.0           shiny_1.7.4                

loaded via a namespace (and not attached):
  [1] backports_1.4.1        Hmisc_4.7-2            plyr_1.8.8            
  [4] splines_4.1.2          BiocParallel_1.28.3    lpsymphony_1.22.0     
  [7] digest_0.6.31          foreach_1.5.2          htmltools_0.5.4       
 [10] GO.db_3.14.0           rsconnect_0.8.29       fansi_1.0.3           
 [13] magrittr_2.0.3         checkmate_2.1.0        memoise_2.0.1         
 [16] cluster_2.1.2          doParallel_1.0.17      Biostrings_2.62.0     
 [19] annotate_1.72.0        askpass_1.1            jpeg_0.1-10           
 [22] colorspace_2.0-3       blob_1.2.3             xfun_0.36             
 [25] dplyr_1.0.10           crayon_1.5.2           RCurl_1.98-1.9        
 [28] jsonlite_1.8.4         genefilter_1.76.0      biglm_0.9-2.1         
 [31] impute_1.68.0          survival_3.2-13        iterators_1.0.14      
 [34] glue_1.6.2             gtable_0.3.1           zlibbioc_1.40.0       
 [37] XVector_0.34.0         DelayedArray_0.20.0    DEoptimR_1.0-11       
 [40] scales_1.2.1           mvtnorm_1.1-3          DBI_1.1.3             
 [43] Rcpp_1.0.10            xtable_1.8-4           htmlTable_2.4.1       
 [46] foreign_0.8-81         bit_4.0.5              preprocessCore_1.56.0 
 [49] Formula_1.2-4          htmlwidgets_1.6.1      httr_1.4.4            
 [52] getopt_1.20.3          RColorBrewer_1.1-3     ellipsis_0.3.2        
 [55] pkgconfig_2.0.3        XML_3.99-0.13          nnet_7.3-16           
 [58] deldir_1.0-6           locfit_1.5-9.7         utf8_1.2.2            
 [61] RJSONIO_1.3-1.7        reshape2_1.4.4         tidyselect_1.2.0      
 [64] rlang_1.0.6            later_1.3.0            AnnotationDbi_1.56.2  
 [67] munsell_0.5.0          tools_4.1.2            cachem_1.0.6          
 [70] cli_3.6.0              generics_0.1.3         RSQLite_2.2.20        
 [73] fastmap_1.1.0          knitr_1.41             bit64_4.0.5           
 [76] robustbase_0.95-0      KEGGREST_1.34.0        nlme_3.1-153          
 [79] mime_0.12              compiler_4.1.2         rstudioapi_0.14       
 [82] curl_5.0.0             png_0.1-8              tibble_3.1.8          
 [85] geneplotter_1.72.0     pcaPP_2.0-3            stringi_1.7.12        
 [88] Matrix_1.5-3           vctrs_0.5.1            pillar_1.8.1          
 [91] lifecycle_1.0.3        optparse_1.7.3         data.table_1.14.6     
 [94] bitops_1.0-7           httpuv_1.6.8           R6_2.5.1              
 [97] latticeExtra_0.6-30    promises_1.2.0.1       gridExtra_2.3         
[100] codetools_0.2-18       MASS_7.3-54            assertthat_0.2.1      
[103] openssl_2.0.5          withr_2.5.0            GenomeInfoDbData_1.2.7
[106] mgcv_1.8-38            parallel_4.1.2         grid_4.1.2            
[109] rpart_4.1-15           base64enc_0.1-3        interp_1.1-3 
