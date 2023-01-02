
path = '../res_v2//rnaSeq_merged/'

type = 'stool'
# type = 'biopsy'

tpm_file = sprintf('%s/TPM_%s.txt', path, type)
taxa_file = sprintf('%s/taxa_%s.txt', path, type)

tpm = read.table(file = tpm_file, header = T, sep = '\t', row.names = 1)
taxa = read.table(file = taxa_file, header = T, sep = '\t', row.names = 1)
pos = which(names(taxa) == 'Description')
taxa = taxa[, (pos+1):dim(taxa)[2] ]

tpm = as.data.frame(t(tpm))
taxa = as.data.frame(t(taxa))

taxa_de_file = sprintf('../python/res/diff-PUCAI_13_a1_%s.txt_feature.txt', type)
taxa_de = read.table(taxa_de_file, header = T, sep = '\t', row.names = 1)

# set taxonomy instead of full sv sequence
source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
taxa = taxa[taxa_de$X_feature_id, ]
sv = row.names(taxa)
tx = sv_2_tax_v2(sv)
tx = tx[order(tx$Feature_ID),]
taxa = taxa[sort(row.names(taxa)),]
row.names(taxa) = tx$clean_name

# tpm_de_file = '/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/Protect/res_PUCAI/Protect_lncRNA_protect_BASELINE_Pucai_Category_1_vs_3_deseq2_res_filtered.tsv'
tpm_de_file = '/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/Protect/res_PUCAI/Protect_ProteinCoding_protect_BASELINE_Pucai_Category_1_vs_3_deseq2_res_filtered.tsv'
tpm_de = read.table(tpm_de_file, header = T, sep = '\t')

tpm = tpm[make.names(tpm_de$Gene[1:30]),]

out_path = sprintf('halla/res/%s/',type)
dir.create(out_path)

taxa_file = sprintf('%s/taxa_pucai_DE.txt', out_path)
# write.table(x = taxa, file = taxa_file, quote = F, sep = '\t', row.names = T)

# tpm_file = sprintf('%s/tpm_pucai_DE.txt', out_path)
tpm_file = sprintf('%s/tpmPC_pucai_DE.txt', out_path)
write.table(x = tpm, file = tpm_file, quote = F, sep = '\t', row.names = T)

