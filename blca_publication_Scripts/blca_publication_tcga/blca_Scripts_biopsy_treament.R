
setwd("/home/jing/Phd_project/project_UCD_blca/blca_DATA/blca_DATA_tcga_pan_can_atlas_2018/blca_DATA_tcga_pan_can_atlas_2018")


clin_raw <- read.delim("data_clinical_patient.txt", sep = '\t',skip = 4)

table(clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE)
#Filtered by stages of interest

clin_raw <- clin_raw[clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE II","STAGE III", "STAGE IV"), ]
#double check
table(clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE)

rownames(clin_raw) <- clin_raw$PATIENT_ID
length(clin_raw$PATIENT_ID)

RNA_raw <- read.delim("data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",check.names = FALSE)
RNA_raw[is.na(RNA_raw)] <- 0
RNA_raw <- RNA_raw[RNA_raw$Hugo_Symbol!='',]
RNA_raw <- RNA_raw[!duplicated(RNA_raw$Hugo_Symbol),]
rownames(RNA_raw) <- RNA_raw$Hugo_Symbol
RNA <- as.data.frame(t(RNA_raw[-1:-2]))

#retrieve RNAs of interest
RNA <- RNA[str_sub(row.names(RNA), end = -4) %in% row.names(clin), ]

clin <- clin_raw[str_sub(row.names(RNA), end = -4),]
clin <- clin[!is.na(clin$OS_MONTHS) & !is.na(clin$OS_STATUS), ]


# create a survival object consisting of times & censoring
surv_obj <- Surv(time = clin$OS_MONTHS, 
                 event = clin$OS_STATUS=="1:DECEASED")

# now creating object without zero times
clin_filt <- clin[clin$OS_MONTHS > 0,]
#RNA_filt<- RNA[str_sub(row.names(RNA), end = -4) %in% row.names(clin_filt), ]

RNA_filt <- RNA[clin$OS_MONTHS > 0,]


# create a survival object consisting of times & censoring
surv_filt <- Surv(time = clin_filt$OS_MONTHS, 
                  event = clin_filt$OS_STATUS=="1:DECEASED")

fit <- survfit(surv_filt ~ 1, data = clin_filt)
ggsurvplot(fit, data = clin_filt, xlab = "Month", ylab = "Overall survival")

#! surv_filt might contain null values because clin_flit contains
which(is.na(surv_filt))


treatment_time <- read.delim2(('data_timeline_treatment.txt'))
biopsy_time <- read.delim2('data_timeline_sample_acquisition.txt')

clin_filt$biopsy_date <- biopsy_time$START_DATE[match(rownames(clin_filt), biopsy_time$PATIENT_ID)]
clin_filt$treatment_date <- treatment_time$START_DATE[match(rownames(clin_filt), treatment_time$PATIENT_ID)]

clin_filt$treated <- clin_filt$treatment_date-clin_filt$biopsy_date
