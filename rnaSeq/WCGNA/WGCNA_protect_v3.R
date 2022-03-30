# install.packages(c("matrixStats", "Hmisc", "splines",
#                    "foreach", "doParallel", "fastcluster",
#                    "dynamicTreeCut", "survival"))
# # source("http://bioconductor.org/biocLite.R")
# BiocManager::install("GO.db" )
# BiocManager::install("preprocessCore" )
# BiocManager::install("impute" )
# # 
# install.packages('/pita/users/tzipi/bin/WGCNA_1.69-81.tar.gz', repos = NULL, lib=.Library)             
library(WGCNA)
options(stringsAsFactors = FALSE);

source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/lncRNA_funcs.R')

type = 'lncRNA'
# type = 'ProteinCoding'
# type = 'lncPrtCd'

# c = c('Protect','SEEM','SOURCE')
c = c('Protect')

if (c == 'Protect' & type %in% c('lncPrtCd', 'ProteinCoding') )
{
  cut_hight = 50000 # 100000 # protet
  pwr = 12 # 6 # protect
}
if (c == 'Protect' & type == 'lncRNA')
{
  cut_hight = 50000 # 100000 # protet
  pwr = 5
}

# set needed tables
metadata_file =  '../metadata/lncRNA_meta_analysis_megamap_v6.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))

TPM_file = sprintf('data/All_cohorts_%s_TPM.txt',type)
tpm = read.table(file = TPM_file, header = T,sep = '\t', stringsAsFactors = F,row.names = 1)

# genes_file = 'res/supps/main3_validated_genes_cohorts_ProteinCoding_supp.txt'
# genes_df = read.table(file = genes_file, header = T,sep = '\t', stringsAsFactors = F,row.names = 1)
# 
# lnc_file = 'res/supps/main3_validated_genes_cohorts_lncRNA_supp.txt'
# lnc_df = read.table(file = lnc_file, header = T,sep = '\t', stringsAsFactors = F,row.names = 1)
# 
# genes_list = row.names(genes_df)[genes_df$Score_UC_CD_Celiac == 'fit']
# lnc_list = row.names(lnc_df)[lnc_df$Score_UC_CD_Celiac == 'fit']
# gene = c(genes_list, lnc_list)
# 
# tpm = tpm[row.names(tpm) %in% gene,]

metadata_df = metadata_df[order(metadata_df$SampleID),]
tpm = tpm[, make.names(sprintf('%s', metadata_df$SampleID))]

cht =  metadata_df$Cohort
tpm = tpm[,cht %in% c]
metadata_df = metadata_df[cht %in% c,]

## filter TPM to 1 TPM>0.2%
tpm_genes = row.names(tpm)
wanted_genes_pos = c()
for ( i in 1:length(tpm_genes) )
{
  per = sum(tpm[i,] > 1) / dim(tpm)[2]
  if ( per > 0.2 ) {wanted_genes_pos = c(wanted_genes_pos, i) }
}
tpm = tpm[wanted_genes_pos, ]

### data input

datExpr0 = as.data.frame(t(tpm));
names(datExpr0) = row.names(tpm);
rownames(datExpr0) = names(tpm);


gsg = goodSamplesGenes(datExpr0, verbose = 3, );
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# cut_hight = 50000 # 100000 # protet
# # cut_hight = 40000 # SEEM
# # cut_hight = 58000 # SOURCE
# Plot a line to show the cut
abline(h = cut_hight, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut_hight, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
metadata_df$SampleID[!keepSamples]
bad_samples = metadata_df$SampleID[make.names(metadata_df$SampleID) %in% row.names(datExpr0)[clust==0] ] 

collectGarbage();

row.names(datExpr) = make.names( row.names(datExpr) )
metadata_df$SampleID = make.names( metadata_df$SampleID )
metadata_df = metadata_df[metadata_df$SampleID %in% row.names(datExpr), ]

metadata_df$healthy = ifelse(metadata_df$Dx=="Control", 0, 1)

metadata_df$Gender = ifelse(metadata_df$Gender=="Female", 2, 1)

metadata_df$Age_years = as.numeric(metadata_df$Age_years)
# # datTraits = metadata_df[,c('Cohort','Dx','Dx_specific','Source','Gender','Age_years','Age_group')]
# datTraits = data.frame(row.names = metadata_df$SampleID,
#                        Age_years = as.numeric(metadata_df$Age_years),
#                        # Dx = metadata_df$Dx,
#                        # Dx_specific = metadata_df$Dx_specific,
#                        Gender = as.numeric(metadata_df$Gender),
#                        healthy = as.numeric(healthy))

protect_metadata = c('Age_years','Gender','protect_BL_BMIZ',
                     'healthy','protect_BASELINE_Pucai',
                     'protect_BASELINE_Pucai_Category',
                     'protect_ENDO_MAYO',
                     'protect_TOTAL_MAYO_EEF',
                     'protect_TOTAL_MAYO_C3',
                     'protect_PARIS_CLASS',
                     # 'protect_Nancy_total_histology',
                     # 'protect_CBH_VILLIFORM',
                     # 'protect_CBH_INFLAM',
                     # 'protect_CBH_EOSCOUNT',
                     # 'protect_EOSGRADE_C2',
                     'protect_ALB_C2',
                     'protect_W4R',
                     # 'protect_W12_CSfreeR',
                     'protect_CSFREE_REMISSION_WK52_KM',
                     'protect_3years_colectomy')
                     
                     # 'protect_Calprotectin',
                     # 'protect_MONTREAL',
                     # 'protect_Overall_WEEK_4_Pucai_Score',
                     # 'protect_Overall_WEEK_4_Pucai_Category','protect_Value_of_week4_PUCAI',
                     # 'protect_CBH_ULCER','protect_CBH_CRYPT',
                     # 'protect_CBH_BPLASMAC','protect_CBH_BLYMPHOID',
                     # 'protect_CBH_PANETH',
                     # 'protect_CBH_EOSGRADE',
                     # 'protect_OMPC','protect_OMPC_over12',
                     # 'protect_INITIAL_TRT_C4')
datTraits = metadata_df[,protect_metadata]
row.names(datTraits) = metadata_df$SampleID

datTraits$protect_PARIS_CLASS = gsub(pattern = 'E',replacement = '',x = datTraits$protect_PARIS_CLASS)
datTraits$protect_PARIS_CLASS = as.numeric(datTraits$protect_PARIS_CLASS)

# datTraits$protect_MONTREAL = gsub(pattern = 'E',replacement = '',x = datTraits$protect_MONTREAL)
# datTraits$protect_MONTREAL = as.numeric(datTraits$protect_MONTREAL)

# datTraits$protect_INITIAL_TRT_C4 = ifelse(datTraits$protect_INITIAL_TRT_C4=="5ASA", 0, 1)

datTraits$protect_CSFREE_REMISSION_WK52_KM = 1 - datTraits$protect_CSFREE_REMISSION_WK52_KM
datTraits$protect_W4R = 1 - datTraits$protect_W4R

# YN_prms = c('protect_CBH_ULCER','protect_CBH_CRYPT','protect_CBH_VILLIFORM','protect_CBH_BPLASMAC','protect_CBH_BLYMPHOID','protect_CBH_PANETH','protect_EOSGRADE_C2','protect_ALB_C2')
# YN_prms = c('protect_CBH_VILLIFORM','protect_EOSGRADE_C2','protect_ALB_C2')
YN_prms = c('protect_ALB_C2')
for ( prm in YN_prms )
{
  datTraits[[prm]] = ifelse(datTraits[[prm]]=="Y", 1, 0)
}


datTraits =datTraits[  ]
# allTraits = allTraits

for ( i in 1:dim(datTraits)[2] )
{
  print(i)
  print( class(datTraits[,i]))
}

colnames(datTraits)[colnames(datTraits) == 'Age_years'] = 'Age'
colnames(datTraits)[colnames(datTraits) == 'protect_BL_BMIZ'] = 'BMI'
colnames(datTraits)[colnames(datTraits) == 'healthy'] = 'protect_Ulcerative colitis'
colnames(datTraits)[colnames(datTraits) == 'protect_BASELINE_Pucai'] = 'protect_PUCAI'
colnames(datTraits)[colnames(datTraits) == 'protect_BASELINE_Pucai_Category'] = 'protect_PUCAI category'
colnames(datTraits)[colnames(datTraits) == 'protect_ENDO_MAYO'] = 'protect_Endoscopic Mayo'
colnames(datTraits)[colnames(datTraits) == 'protect_TOTAL_MAYO_EEF'] = 'protect_Total Mayo'
colnames(datTraits)[colnames(datTraits) == 'protect_TOTAL_MAYO_C3'] = 'protect_Total Mayo category'
colnames(datTraits)[colnames(datTraits) == 'protect_PARIS_CLASS'] = 'protect_Paris classification'
## remove # colnames(datTraits)[colnames(datTraits) == 'protect_Nancy_total_histology'] = 'protect_Baseline Nancy histology score'
# colnames(datTraits)[colnames(datTraits) == 'protect_CBH_VILLIFORM'] = 'protect_Baseline histology villiform changes'
# colnames(datTraits)[colnames(datTraits) == 'protect_CBH_INFLAM'] = 'protect_Baseline histology inflammation'
# colnames(datTraits)[colnames(datTraits) == 'protect_CBH_EOSCOUNT'] = 'protect_Baseline histology eosinophils count'
# colnames(datTraits)[colnames(datTraits) == 'protect_EOSGRADE_C2'] = 'protect_Baseline histology eosinophils grade'
colnames(datTraits)[colnames(datTraits) == 'protect_ALB_C2'] = 'protect_Albumin'
colnames(datTraits)[colnames(datTraits) == 'protect_W4R'] = 'protect_W4 no response'
# colnames(datTraits)[colnames(datTraits) == 'protect_W12_CSfreeR'] = 'protect_W12 CS free remsion '
colnames(datTraits)[colnames(datTraits) == 'protect_CSFREE_REMISSION_WK52_KM'] = 'protect_W52 no CS free remssion'
colnames(datTraits)[colnames(datTraits) == 'protect_3years_colectomy'] = 'protect_Colectomy within 3Y'

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
# pdf(sprintf('res/WGCNA/%s_sample_dendro_%s.pdf',c,type),width = 13,height = 12 )
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
# dev.off()
# save(datExpr, datTraits, file = "try-dataInput.RData")


#### netweork construction

# lnames = load(file = 'try-dataInput.RData');
#The variable lnames contains the names of loaded variables.
# lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# pwr = 12 # 6 # protect
# # pwr = 10 # SEEM
# # pwr = 7 # SOURCE
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = pwr, maxBlockSize = 20000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       # saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)

table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
# pdf(sprintf('res/WGCNA/%s_cluster_dendrogram_%s.pdf',c,type),width = 10,height = 5 )
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
# save(MEs, moduleLabels, moduleColors, geneTree, file = "try-networkConstruction-auto.RData")


# #### relate models to metadata
# lnames = load(file = "Protect_try-dataInput.RData");
# #The variable lnames contains the names of loaded variables.
# lnames
# # Load network data saved in the second part.
# lnames = load(file = "protect_try-networkConstruction-auto.RData");
# lnames

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3.1, 3));
# Display the correlation values within a heatmap plot

rm_gray_flag = T
if (rm_gray_flag)
{
  textMatrix_old = textMatrix
  textMatrix = textMatrix[row.names(moduleTraitCor)!='MEgrey',]
  moduleTraitCor_old = moduleTraitCor
  moduleTraitCor = moduleTraitCor[row.names(moduleTraitCor)!='MEgrey',]
  moduleTraitPvalue_old = moduleTraitPvalue
  moduleTraitPvalue = moduleTraitPvalue[row.names(moduleTraitPvalue)!='MEgrey',]
  MEs_old = MEs
  MEs = MEs[,names(MEs)!='MEgrey']
}

# jpeg(sprintf('res/WGCNA/%s_module_trait_%s_filtered01.jpeg',c,type),width = 1000,height = 600)
# pdf(sprintf('res/WGCNA/%s_module_trait_%s.pdf',c,type),width = 8,height = 6 )
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = gsub(pattern = 'protect_',replacement = '',x = names(datTraits)),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               xLabelsAdj = 0.9,
               cex.lab=0.9,
               xLabelsAngle = 30,
               main = paste("Module-trait relationships"))
# dev.off()

## filter heatmap to significant in healthy/sick only
good_cols = row.names(moduleTraitPvalue)[ moduleTraitPvalue[,'protect_Ulcerative colitis'] <= 0.001 ]
textMatrix_f =  paste(signif(moduleTraitCor[good_cols,], 2), "\n(",
                      signif(moduleTraitPvalue[good_cols,], 1), ")", sep = "");

h=7
symb = gsub('ME','',good_cols)
if ( type == 'lncRNA' )
{
  symb = c('M1','M2','M3','M4','M5')
  h=6
}

pdf(sprintf('code_for_paper/rnaSeq/figs/fig1/%s_module_trait_%s_filtered001.pdf',c,type),width = 9,height = h )
labeledHeatmap(Matrix = moduleTraitCor[good_cols,],
               xLabels = gsub(pattern = 'protect_',replacement = '',x = names(datTraits)),
               yLabels = good_cols, 
               ySymbols =  symb,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_f,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               xLabelsAdj = 1,
               cex.lab=1.5,
               xLabelsAngle = 30,
               main = paste("Module-trait relationships"))
dev.off()

# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$Gender);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

ADJ1=abs(cor(datExpr,use="p"))^pwr
Alldegrees=intramodularConnectivity(ADJ1, moduleColors)

# head(Alldegrees_protect)
# module = "blue"
# column = match(module, modNames);
# moduleGenes = moduleColors==module;
# 
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = "Gene significance for body weight",
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# 
# names(datExpr)
# names(datExpr)[moduleColors=="blue"]

df = data.frame(row.names = names(datExpr), module_color = moduleColors)

df = merge(x = df, y = Alldegrees, by = 0); rownames(df)=df$Row.names; df=df[2:dim(df)[2]];
df = merge(x = df, y = geneTraitSignificance, by = 0); rownames(df)=df$Row.names; df=df[2:dim(df)[2]];
df = merge(x = df, y = GSPvalue, by = 0); rownames(df)=df$Row.names; df=df[2:dim(df)[2]];
df = merge(x = df, y = geneModuleMembership, by = 0); rownames(df)=df$Row.names; df=df[2:dim(df)[2]];
df = merge(x = df, y = MMPvalue, by = 0); rownames(df)=df$Row.names; df=df[2:dim(df)[2]];

if (type == 'lncPrtCd')
{
  # add gene type
  prt_file = '/pita/users/tzipi/projects/rnaSeq/lncRNA/tables/from_yael/19815_Protein_coding_genes_GencodeV24_fixed.txt'
  prt_df = read.table(prt_file, header = T)
  
  lnc_file = '/pita/users/tzipi/projects/rnaSeq/lncRNA/tables/from_yael/lncRNA_gencode24_forR.txt'
  lnc_df = read.table(lnc_file, header = T)
  
  gene_type = vector(mode = 'character',length = dim(df)[1])
  gene_type[row.names(df) %in% prt_df$Gene] = 'Protein_coding'
  gene_type[row.names(df) %in% lnc_df$Gene] = 'lncRNA'
  df$Gene_type = gene_type
}
df$Gene = row.names(df)
l = dim(df)[2]
df = df[,c(l,1:(l-1))]

# write.table(x = df, file = sprintf('res/WGCNA/%s_gene_module_table_%s.txt',c, type),quote = F,row.names = F,sep = '\t')
