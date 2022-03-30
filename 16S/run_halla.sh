INPATH=$'res_v2/'

TYPE=$'stool' ## TYPE=$'biopsy' 

halla \
      -x $INPATH'/'$TYPE/'taxa_pucai_DE.txt' \
      -y $INPATH'/'$TYPE/'taxa_pucai_DE.txt' \
      -o $INPATH'/'$TYPE/'tpmDE_taxaDE_res_a25/' \
      -m spearman \
      --fdr_alpha 0.1 &


