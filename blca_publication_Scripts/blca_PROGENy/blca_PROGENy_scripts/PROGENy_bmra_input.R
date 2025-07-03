library(progeny)
library(tidyr)
library(tibble)
library(openxlsx)

wkdir <- '/home/jing/Phd_project/project_UCD_blca/blca_publication/blca_publication_Scripts/blca_PROGENy/blca_PROGENy_DATA'
setwd(wkdir)

bc3c_sc <-read.xlsx('PROGENy_Data_mean_2020_BC3C.xlsx',rowNames = TRUE)
head(colnames(bc3c_sc))

intersection_non_inf_df <- read.xlsx('non_inferred_gene_df.xlsx',rowNames = TRUE)
length(rownames(progeny_gene_list))

intersect(colnames(bc3c_sc),rownames(intersection_non_inf_df))

sig_info <- read.xlsx('BC3C_sig_info_df_drugs_full.xlsx',rowNames = TRUE)

#Add three drugs cabozatinib, EMD-1214063 and crizotinib
all_sig <- read.xlsx('/home/jing/Phd_project/project_UCD_blca/blca_DATA/blca_DATA_LINCS_output/00_outputs_2020_BC3C/sig_info_2020_BC3C.xlsx',
                     rowNames = TRUE)
selected_drugs <- c('cabozantinib', 'EMD-1214063', 'crizotinib')
all_sig_met <- all_sig[all_sig$pert_drug %in% selected_drugs, ]


sig_info <- rbind(sig_info, all_sig_met)
library(progeny)

bc3c_sc_non_inf <- bc3c_sc[,rownames(intersection_non_inf_df)]

bc3c_sc_non_inf_t <- t(bc3c_sc_non_inf)
bc3c_sc_non_inf_t_selected <- bc3c_sc_non_inf_t[,rownames(sig_info)]

progeny_scores <-progeny(bc3c_sc_non_inf_t_selected,scale=TRUE)
head(progeny_scores)


#Full 
full_progeny_scores <-progeny(bc3c_sc_non_inf_t,scale=TRUE)

hist(full_progeny_scores, main = ' ')

drugs_check <- names(head(sort(full_progeny_scores[, "TGFb"],decreasing = F),100))
write.csv(all_sig[rownames(all_sig) %in% drugs_check,],'/home/jing/Phd_project/project_UCD_blca/blca_publication/blca_publication_Scripts/blca_PROGENy/blca_PROGENy_OUTPUT/top_tgfb_inh.csv')

top_tgfb <- sort(full_progeny_scores[, "TGFb"], decreasing = FALSE)
top_tgfb_df <- data.frame(sample_id = names(top_tgfb), score = top_tgfb)

top_tgfb_df <- head(top_tgfb_df, 100)
all_sig$sample_id <- rownames(all_sig)

merged_df <- merge(top_tgfb_df, all_sig, by = "sample_id")
write.csv(merged_df, '/home/jing/Phd_project/project_UCD_blca/blca_publication/blca_publication_Scripts/blca_PROGENy/blca_PROGENy_OUTPUT/top_tgfb_inh.csv', row.names = FALSE)


pheatmap(full_progeny_scores,cluster_rows= TRUE, cluster_columns=TRUE,show_rownames = ,
         show_colnames = TRUE)

progeny_scores[rownames(all_sig_met),][,'TGFb']

library(pheatmap)
pheatmap(progeny_scores[rownames(all_sig_met),],cluster_rows= TRUE, cluster_columns=TRUE,show_rownames = ,
         show_colnames = TRUE)

pheatmap(progeny_scores,cluster_rows= TRUE, cluster_columns=TRUE,show_rownames = ,
         show_colnames = TRUE)


sig_info$first_target <-first_target
sig_info<-sig_info[order(sig_info$first_target),]
first_target <- str_sub(sig_info$targets,end=4)

progeny_scores<-progeny_scores[row.names(sig_info),]

reference_names <- paste(rownames(sig_info),first_target,sep='_')
rownames(progeny_scores) <-reference_names
pheatmap(progeny_scores,cluster_rows= FALSE, cluster_columns=TRUE,show_rownames =TRUE ,
         show_colnames = TRUE)
