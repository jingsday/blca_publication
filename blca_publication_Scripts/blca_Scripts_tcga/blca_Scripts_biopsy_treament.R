library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)

setwd("/home/jing/Phd_project/project_UCD_blca/blca_DATA/blca_DATA_tcga_pan_can_atlas_2018/blca_DATA_tcga_pan_can_atlas_2018")

clin_raw <- read.delim("data_clinical_patient.txt", sep = '\t',skip = 4)
table(clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE)

#Filtered by stages of interest
clin_raw <- clin_raw[clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE II","STAGE III", "STAGE IV"), ]

#double check
table(clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE)
rownames(clin_raw) <- clin_raw$PATIENT_ID

RNA_raw <- read.delim("data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",check.names = FALSE)
RNA_raw[is.na(RNA_raw)] <- 0
RNA_raw <- RNA_raw[RNA_raw$Hugo_Symbol!='',]
RNA_raw <- RNA_raw[!duplicated(RNA_raw$Hugo_Symbol),]
rownames(RNA_raw) <- RNA_raw$Hugo_Symbol
RNA <- as.data.frame(t(RNA_raw[-1:-2]))

#retrieve RNAs of interest

clin <- clin_raw[str_sub(row.names(RNA), end = -4),]
clin <- clin[!is.na(clin$OS_MONTHS) & !is.na(clin$OS_STATUS), ]

RNA <- RNA[str_sub(row.names(RNA), end = -4) %in% row.names(clin), ]

# create a survival object consisting of times & censoring
surv_obj <- Surv(time = clin$OS_MONTHS, 
                 event = clin$OS_STATUS=="1:DECEASED")

# now creating object without zero times
clin_filt <- clin[clin$OS_MONTHS > 0,]
#RNA_filt<- RNA[str_sub(row.names(RNA), end = -4) %in% row.names(clin_filt), ]

RNA_filt <- RNA[clin$OS_MONTHS > 0,]

#Treatment and biopsy dataframes
treatment_time <- read.delim2(('data_timeline_treatment.txt'))
rownames(treatment_time) <- treatment_time$PATIENT_ID

biopsy_time <- read.delim2('data_timeline_sample_acquisition.txt')
length(intersect(clin_filt$PATIENT_ID,treatment_time$PATIENT_ID)) #127 overlapping

#Removing patients that received immuno 
#and those who were treated before biopsies were taken

clin_filt$biopsy_date <- biopsy_time$START_DATE[match(rownames(clin_filt), biopsy_time$PATIENT_ID)]
clin_filt$treatment_date <- treatment_time$START_DATE[match(rownames(clin_filt), treatment_time$PATIENT_ID)]
clin_filt$treated <- clin_filt$treatment_date-clin_filt$biopsy_date

rownames(clin_filt[rownames(clin_filt)%in% treatment_time[which(treatment_time$TREATMENT_TYPE == 'Immunotherapy'),'PATIENT_ID'],])
p_id <- append(treatment_time[which(treatment_time$TREATMENT_TYPE == 'Immunotherapy'),'PATIENT_ID'],c('TCGA-XF-AAN1','TCGA-K4-AAQO','TCGA-XF-A9SV','TCGA-XF-AAMZ'))


excluded_df <- clin_filt[p_id, ]
excluded_df <- merge(excluded_df,treatment_time,how='left',on='PATIENT_ID')
#write.csv(excluded_df,file='/home/jing/Phd_project/project_UCD_blca/blca_publication_OUTPUT/blca_publication_OUTPUT_survival/excluded_df.csv')
