library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)
library(maftools)
library(ggsurvplot)

setwd('/home/jing/Phd_project/project_UCD_blca/blca_DATA/blca_DATA_tcga_pan_can_atlas_2018/blca_DATA_tcga_pan_can_atlas_2018')
blca.maf <-read.delim("data_mutations.txt", sep = '\t')

head(blca.maf$Tumor_Sample_Barcode)
head(blca.maf$Matched_Norm_Sample_Barcode)


clin_raw<- read.delim("data_clinical_patient.txt", sep = '\t',skip = 4,row.names = 'PATIENT_ID')

table(clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE)
#Filtered by stages of interest

clin_raw <- clin_raw[clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE II","STAGE III","STAGE IV"), ]
#double check
table(clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE)



RNA_raw <- read.delim("data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",check.names = FALSE)
RNA_raw[is.na(RNA_raw)] <- 0
RNA_raw <- RNA_raw[RNA_raw$Hugo_Symbol!='',]
RNA_raw <- RNA_raw[!duplicated(RNA_raw$Hugo_Symbol),]
rownames(RNA_raw) <- RNA_raw$Hugo_Symbol
RNA <- as.data.frame(t(RNA_raw[-1:-2]))

# Align clinical data:
clin <- clin_raw[row.names(clin_raw)%in% str_sub(row.names(RNA), end = -4),]
RNA <- RNA[str_sub(row.names(RNA), end = -4) %in% row.names(clin_raw), ]

clin= clin[str_sub(row.names(RNA), end = -4) %in% row.names(clin_raw),]
clin$Tumor_Sample_Barcode =paste0(rownames(clin),'-01')
clin$Tumor_Sample_Barcode

clin_filt <- clin[clin$OS_MONTHS > 0,]
RNA_filt <- RNA[clin$OS_MONTHS > 0,]

surv_filt <- Surv(time = clin_filt$OS_MONTHS, 
                  event = clin_filt$OS_STATUS=="1:DECEASED")
which(is.na(surv_filt))
surv_filt<- surv_filt[-261,]
clin_filt<- clin_filt[-261,]
RNA_filt<- RNA_filt[-261,]

#maf_data_view 

blca <-read.maf(maf=blca.maf,clinicalData= clin_filt)

# Filter MAF data to include only matching samples
matching_barcodes <- intersect(blca.maf$Tumor_Sample_Barcode, clin_filt$Tumor_Sample_Barcode)

blca<- subsetMaf(maf = blca, tsb = matching_barcodes)


maf_data_view <- blca@data
maf_gene.summary <-blca@gene.summary

getSampleSummary(blca)


plotmafSummary(maf=blca,rmOutlier = TRUE,addStat = 'median',dashboard = TRUE,
               titvRaw=FALSE)

mafbarplot(maf=blca)

oncoplot(maf=blca,top=50)

blca.titv = titv(maf=blca,plot=FALSE,useSyn=TRUE)
plotTiTv(res=blca.titv)  

#individual genes TP53, ERCC1,FGFR3
lollipopPlot(maf=blca,
             gene='FGFR3',
             AACol = 'HGVSp_Short',
             showMutationRate = TRUE)

lollipopPlot(maf = blca,
             gene = 'FGFR3',
             AACol = 'HGVSp_Short',
             showMutationRate = TRUE,
             refSeqID = 'NM_022965') # Specify the transcript you want to use

plotProtein(gene="ERCC1")

rainfallPlot(maf = blca, detectChangePoints = TRUE, pointSize = 0.4)

plotVaf(maf=blca)

getClinicalData(blca)

# prepare for survival analysis

blca@clinical.data$deceased <- ifelse(blca@clinical.data$OS_STATUS == "0:LIVING", 1, 0)

blca@clinical.data[,OS_MONTHS :=  as.numeric(OS_MONTHS)]
blca@clinical.data[,OS_MONTHS :=  as.numeric(OS_MONTHS)]

par(mar = c(1, 1, 1, 1))
plot(1:30)
#Survival analysis
mafSurvival(maf=blca, genes='KMT2C',time ='OS_MONTHS',Status = 'deceased')


#
fit<- survfit(Surv(clin$OS_MONTHS, clin$OS_STATUS == '1:DECEASED') ~ SEX, data = clin)
ggsurvplot(fit, data = clin,
           risk.table = TRUE,                  # Add No at risk table
           cumevents = TRUE,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE    ,pval = TRUE)

km_fit <-survfit(Surv(clin$OS_MONTHS, clin$OS_STATUS == '1:DECEASED') ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = clin)
ggsurvplot(km_fit,data=clin, risk.table = TRUE,pval = TRUE)
km_fit
fit <- survfit(surv_obj ~ 1, data = clin)
ggsurvplot(fit, data = clin, xlab = "Month", ylab = "Overall survival")
fit

##Predict genesets with survivals
prog_geneset = survGroup(maf = blca, top = 20, geneSetSize = 2, time = 'OS_MONTHS', Status = 'deceased', verbose = FALSE)

print(head(prog_geneset))


mafSurvGroup(maf = blca, geneSet = c("", ""), time = "days_to_last_followup", Status = "Overall_Survival_Status")

