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

sig_info <- read.xlsx('sig_info_2020_BC3C.xlsx',rowNames = TRUE)

library(progeny)

bc3c_sc_non_inf <- bc3c_sc[,rownames(intersection_non_inf_df)]
bc3c_sc_non_inf_t <- t(bc3c_sc_non_inf)


#Full 
full_progeny_scores <-progeny(bc3c_sc_non_inf_t,scale=TRUE)

hist(full_progeny_scores, main = ' ')

drugs_check <- names(head(sort(full_progeny_scores[, "TGFb"],decreasing = F),100))
#write.csv(all_sig[rownames(all_sig) %in% drugs_check,],'/home/jing/Phd_project/project_UCD_blca/blca_publication/blca_publication_Scripts/blca_PROGENy/blca_PROGENy_OUTPUT/top_tgfb_inh.csv')

top_tgfb <- sort(full_progeny_scores[, "TGFb"], decreasing = FALSE)
top_tgfb_df <- data.frame(sample_id = names(top_tgfb), score = top_tgfb)

top_tgfb_df <- head(top_tgfb_df, 100)
sig_info$sample_id <- rownames(sig_info)

merged_df <- merge(top_tgfb_df, sig_info, by = "sample_id")

#write.csv(merged_df, '/home/jing/Phd_project/project_UCD_blca/blca_publication/blca_publication_Scripts/blca_PROGENy/blca_PROGENy_OUTPUT/top_tgfb_inh.csv', row.names = FALSE)

library(pheatmap)
pheatmap(full_progeny_scores,cluster_rows= TRUE, cluster_columns=TRUE,show_rownames = ,
         show_colnames = TRUE)

progeny_scores[rownames(all_sig_met),][,'TGFb']
