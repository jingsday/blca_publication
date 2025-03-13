library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)

setwd("/Users/lidiayung/PhD_project/project_UCD_blca/blca_DATA/blca_DATA_tcga_pan_can_atlas_2018")


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
                 event = clin$OS_STATUS=="1:DECEASED",)
#surv_obj 

fit <- survfit(surv_obj ~ 1, data = clin)
ggsurvplot(fit, data = clin, xlab = "Month", ylab = "Overall survival",surv.median.line = "hv")


fit_stage <- survfit(surv_obj ~ clin$AJCC_PATHOLOGIC_TUMOR_STAGE, data = clin)
ggsurvplot(fit_stage,data = clin,pval=T)
t <- ggsurvplot(
  fit_stage,
  data = clin,
  pval = TRUE,
  xlab = 'Month',
  legend.title = "STAGES",
  legend.labs = c("STAGE=II", "STAGE=III", "STAGE IV"),
  legend.size = 1.5,
  legend = "top",
  risk.table = TRUE
)

# Saving the plot
ggsave("~/Downloads/survival_plot.png", dpi = 300, plot = t, width = 4, height = 3, device = "png")


# fit multivariate model (COX proportional hazard) 
fit.coxph <- coxph(surv_obj ~ EGFR + ERBB2 + ERBB3 + ERBB4 + NRG1 + NRG2 + NRG3 + NRG4, 
                   data = RNA)
ggforest(fit.coxph, data = RNA)
#new.coxph <- predict(fit.coxph, newdata = RNA[,c("KRAS", "ERBB2")])
#new.coxph
fit.coxph <- coxph(surv_obj ~ KSR1 + KSR2 + IQGAP1 + GAB1 + GAB2, 
                   data = RNA)
ggforest(fit.coxph, data = RNA)




library(progeny)
zscores = as.matrix(t(RNA))
pathways <- progeny(zscores, scale=TRUE, organism="Human")  #, top = 100, perm = 1)
path_df = as.data.frame(pathways)
# fit multivariate model (COX proportional hazard) 
fit.coxph <- coxph(surv_obj ~ MAPK + PI3K  +  p53 + TGFb  + `JAK-STAT` + TNFa, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)

fit.coxph <- coxph(surv_obj ~ NFkB + Androgen +  VEGF + `JAK-STAT` + TGFb, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)

fit.coxph <- coxph(surv_obj ~ MAPK + PI3K  +  p53 + TGFb  + `JAK-STAT` + TNFa + 
                     WNT + Androgen +  Estrogen + Hypoxia +  Trail + VEGF, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)


library(corrplot)
corrplot(cor(path_df), type = "upper", order = "hclust",tl.col = "black", tl.srt = 45)

path_df

TGFb_df = cbind(TGFb = path_df$TGFb, SMAD2 = RNA$SMAD2, 
                        SMAD3 = RNA$SMAD3, SMAD4 = RNA$SMAD4, 
                        SMAD7 = RNA$SMAD7, TGFB1 = RNA$TGFB1, 
                        TGFB2 = RNA$TGFB2, TGFB3 = RNA$TGFB3, 
                        TGFBR1 = RNA$TGFBR1, TGFBR2 = RNA$TGFBR2, 
                        TGFBR3 = RNA$TGFBR3, ACVR1B = RNA$ACVR1B, 
                        ACVR1C = RNA$ACVR1C, ACVR2A = RNA$ACVR2A, 
                        ACVR2B = RNA$ACVR2B)

corrplot(cor(TGFb_df), type = "upper", order = "hclust",tl.col = "black", tl.srt = 45)

# analysing on what MAPK activity depend
MAPK_df = cbind(MAPK = path_df$MAPK, PI3K = path_df$PI3K,
                KSR1 = RNA$KSR1, KSR2 = RNA$KSR2, IQGAP1 = RNA$IQGAP1, 
                IQGAP2 = RNA$IQGAP2, IQGAP3 = RNA$IQGAP3, GAB1 = RNA$GAB1, GAB2 = RNA$GAB2,
                KRAS = RNA$KRAS, NRAS = RNA$NRAS, HRAS = RNA$HRAS, BRAF=RNA$BRAF, 
                CRAF=RNA$RAF1, ARAF=RNA$ARAF)
corrplot(cor(MAPK_df), type = "upper", order = "hclust",tl.col = "black", tl.srt = 45)

cor(MAPK_df)




# now creating object without zero times
clin_filt <- clin[clin$OS_MONTHS > 0,]
#RNA_filt<- RNA[str_sub(row.names(RNA), end = -4) %in% row.names(clin_filt), ]

RNA_filt <- RNA[clin$OS_MONTHS > 0,]

gene_lincs <- read.csv("/Users/lidiayung/PhD_project/project_UCD_blca/blca_DATA/blca_DATA_tcga_pan_can_atlas_2018/gl_lincs.csv")

RNA_lincs <- RNA_filt[, colnames(RNA_filt) %in% gene_lincs$gene]


path_filt <- path_df[clin$OS_MONTHS > 0,]
# create a survival object consisting of times & censoring
surv_filt <- Surv(time = clin_filt$OS_MONTHS, 
                 event = clin_filt$OS_STATUS=="1:DECEASED")

fit <- survfit(surv_filt ~ 1, data = clin_filt)
ggsurvplot(fit, data = clin_filt, xlab = "Month", ylab = "Overall survival")



#! surv_filt might contain null values because clin_flit contains
which(is.na(surv_filt))

#Remove rows that have null values
surv_filt<- surv_filt[-261,]
RNA_lincs<- RNA_lincs[-261,]


# glmnet
library("glmpath")
library("glmnet")
library("penalized")


fit_glm <- glmnet(RNA_lincs,surv_filt,family="cox")#, alpha = 1) #, alpha = 1, standardize = TRUE, maxit = 1000
print(fit_glm)
# analysing results

cfs = coef(fit_glm,s=0.001717)
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]

coef_data <-data.frame(variable=meaning_coefs,coefficient=meaning_vals)
rownames(coef_data)

sorted_coef_abs <- coef_data[order(-abs(coef_data$coefficient)), ]

top_100 <- head(sorted_coef_abs,100)
list(top_100$variable)

# Assuming top_100 is a data frame containing the top 100 coefficients
paste(top_100$variable, collapse = ",")

# Create bar plot
ggplot(coef_data, aes(x = reorder(variable, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "Genes", y = "Coefficient") +
  theme_minimal()

#write.csv(meaning_vals,file='97th_survival_coeffs_glmnet_blca.csv')
#library(openxlsx)
#write.xlsx(coef_data, "97th_survival_coeffs_glmnet_blca.xlsx")

print(paste(meaning_coefs,collapse=" + "))
# cutting to find most important
ncut = 5
vals_surv = sort(abs(meaning_vals),decreasing = TRUE)[1:ncut]
print(paste(names(vals_surv),collapse=" + "))
# if we want to run coxph we have to copy-paste the string
fit.coxph <- coxph(surv_filt ~ RAB27A + STAP2 + RRP12 + OXCT1 + SATB1 , data = RNA_lincs)
ggforest(fit.coxph, data = RNA_filt)

ggplot(coef_data, aes(x = reorder(variable, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "Genes", y = "Coefficient") +
  theme_minimal()


#Cross validation
cvfit <- cv.glmnet(data.matrix(RNA_lincs),surv_filt,family="cox",type.measure = "C")
plot(cvfit)



#DPD hazard ratio 

#rna data log transformed

rna_raw <- read.delim("data_mrna_seq_v2_rsem.txt",check.names = FALSE)

rna_raw[is.na(rna_raw)] <- 0
rna_raw <- rna_raw[rna_raw$Hugo_Symbol!='',]
rna_raw <- rna_raw[!duplicated(rna_raw$Hugo_Symbol),]
rownames(rna_raw) <- rna_raw$Hugo_Symbol
rna <- as.data.frame(t(rna_raw[-1:-2]))

#retrieve RNAs of interest
rna <- rna[str_sub(row.names(rna), end = -4) %in% row.names(clin_raw), ]
#rna <- rna[-261,]
rna_log2 <- apply(rna, c(1, 2), function(value) log2(value + 1))



max(rna_log2)
min(rna_log2)
library(openxlsx)
stvs <- read.xlsx('/Users/lidiayung/PhD_project/project_UCD_blca/blca_OUTPUT/blca_OUTPUT_depmap/blca_OUTPUT_depmap_stvs.xlsx')

head(colnames(rna_log2))
head(stvs$Gene)

missing_genes <- setdiff(stvs$Gene, colnames(rna_log2))
if(length(missing_genes) > 0) {
  warning("The following genes are not found in rna_log2 columns: ", paste(missing_genes, collapse = ", "))
}

cols<-intersect(stvs$Gene, colnames(rna_log2))
# Subset the matrix
rna_log2_subset <- rna_log2[, cols, drop = FALSE]

#dot product
stvs <- stvs[!stvs$Gene %in% missing_genes,]


rna_log2_subset <- rna_log2[, stvs$cell_lines]
coef_data_subset <-coef_data[,-1]



# Perform the dot product
dot_product_result <- rna_log2_subset %*% stvs$blca_onc
head(stvs$cell_lines)

result <- rna_log2_subset %*% stvs$blca_inv
head(dot_product_result)

clin$DPD_prognosis <- dot_product_result[str_sub(row.names(dot_product_result), end = -4) %in% row.names(clin), ]


table(clin$AJCC_PATHOLOGIC_TUMOR_STAGE)

#write.csv(clin$DPD_prognosis, file = "DPD_prognosis.csv", row.names = FALSE)

clin$outcome <- ifelse(clin$DPD_prognosis > 10, "positive",
                       ifelse(clin$DPD_prognosis == 10, "0", "negative"))

### hazard ratio
clin_nice <- clin
fit2.coxph <- coxph(Surv(time = clin_nice$OS_MONTHS, 
                        event = clin_nice$OS_STATUS=="1:DECEASED") ~ dot_product_result+SEX+AJCC_PATHOLOGIC_TUMOR_STAGE, data = clin_nice)

ggforest(fit2.coxph, data = clin_nice,)
ggsave('hr_dotproduct.svg',dpi=72)
ggsave('hr_dotproduct2.svg',dpi=300)


##kaplan M curve
fit <- survfit(Surv(time = clin$OS_MONTHS, 
                    event = clin$OS_STATUS=="1:DECEASED") ~ outcome, data = clin)
ggsurvplot(fit, data = clin, pval = TRUE,xlab = "Month", ylab = "Overall survival",risk.table=TRUE)



# Check
table(clin$outcome,clin$AJCC_PATHOLOGIC_TUMOR_STAGE)

#Overall expression 

#keep top and bottom genes 
library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
library(ggplotify) ## to convert pheatmap to ggplot2
library(heatmaply) ## for constructing interactive heatmap

pheatmap(RNA)

# Sort the genes based on absolute z-scores
DE_results_sorted <- RNA[order(abs(RNA), decreasing = TRUE), ]

# Select the top 20 genes
top_20_genes <- DE_results_sorted[1:20, ]


# Assuming `exprData` is a matrix where rows represent genes and columns represent patients/groups
heatmap(exprData, Rowv = , Colv =exprData , scale = "none", main = "Gene Expression Heatmap")

  
boxplot(exprData ~ clin$outcome, data = exprData, main = "Gene Expression Boxplot", xlab = "Patient Group", ylab = "Expression Level")


fit.coxph <- coxph(surv_obj ~ MAPK + PI3K  +  p53 + TGFb  + `JAK-STAT` + TNFa + 
                     WNT + Androgen +  Estrogen + Hypoxia +  Trail + VEGF, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)






path_filt <- path_filt[-261,]

# glmnet over progeny
fit_glm <- glmnet(path_filt,surv_filt,family="cox") # , alpha = 1, standardize = TRUE, maxit = 1000
print(fit_glm)
# analysing results
cfs = coef(fit_glm,s=0.02)
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]


coef_data <-data.frame(variable=meaning_coefs,coefficient=meaning_vals)
# Create bar plot
ggplot(coef_data, aes(x = reorder(variable, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "Pathways", y = "Coefficient") +
  theme_minimal()

write.csv(meaning_vals,file='skcm_survival_coeffs_progeny_glmnet.csv')
print(paste(meaning_coefs,collapse=" + "))
# cutting to find most important
ncut = 10
vals_surv = sort(abs(meaning_vals),decreasing = TRUE)[1:ncut]
print(paste(names(vals_surv),collapse=" + "))
# if we want to run coxph we have to copy-paste the string
fit.coxph <- coxph(surv_filt ~ TGFb+Estrogen+VEGF+PI3K+Androgen+EGFR+WNT+p53 +`JAK-STAT`+ NFkB, data = path_filt)
ggforest(fit.coxph, data = RNA_lincs)

#`JAK-STAT`

# plotting Kaplan-Mayer curves
pathway = 'MAPK'
pathway_data = path_filt$MAPK

MDM2_df <- dpd_bc3c[dpd_bc3c$targets == 'MDM2',]
pathway_data <-MDM2_df$x

# sort age 
uni_path = sort(unique(pathway_data))

# store results
results_path = matrix(1,length(uni_path))
# do a for loop for every unique value of age mat
for (i in 2:(length(uni_path)-1)){ # Starting from 2 because the first element would yield no elements higher than.
  path_i = 1*(pathway_data>uni_path[i])
  # survdiff is the function from the survival package 
  logrank = survdiff(surv_filt ~ path_i)
  # store in results_age
  results_path[i] = logrank$pvalue
}
# Plot unique elements of age against p-value
plot(uni_path, results_path, log = 'y')
# Select minimum P-value
min_p_path = which.min(results_path)
# here are 2 good thresholds, -1 and 1
opt_thr = uni_path[min_p_path]
opt_thr
#opt_JAK = opt_thr
#opt_thr = 0.0
# I recalculated the P-value as a sanity check
pval = survdiff(surv_filt ~ pathway_data>opt_thr)$pvalue
nplus = sum(pathway_data>opt_thr)   # how many patients we have in high group
nminus = sum(pathway_data<opt_thr)   # how many patients we have in low group
# fit Kaplan Meier model
path_filt <- path_filt %>% mutate(MAPK_group = ifelse(MAPK >= 1, "high", ifelse(MAPK < 1.001*opt_thr,"low","intermediate")))
nhigh = sum(pathway_data>1)
ninter = sum((pathway_data<1) & (pathway_data > 1.001*opt_thr))
nlow = sum(pathway_data<1.001*opt_thr)
#KM = survfit(surv_filt ~ pathway_data>opt_thr,data = path_filt)
KM = survfit(surv_filt ~ MAPK_group,data = path_filt)
# Plot Kaplan Meier 
#plot(KM, lwd = 3, col = c(1,2), cex.axis = 1.5, xlab = 'Months', ylab = 'Survival Probability' , cex.lab = 1.5)
#ggsurvplot(KM, data = path_filt,pval = TRUE,xlab = 'Overall survival time, months',
#           legend.labs=c(paste('Low ',pathway,' activity, ',nminus,' patient(s)',sep = ''),paste('High ',pathway,' activity, ',nplus,' patient(s)',sep = '')),
#           palette = c('blue','red'),legend.title="")
p <- ggsurvplot(KM, data = path_filt,pval = TRUE,xlab = 'Overall survival time, months',
           legend.labs=c(paste("High MAPK activity,\n",nhigh," patients",sep=""),paste("Intermediate MAPK activity,\n",ninter," patients",sep=""),paste("Low MAPK activity,\n",nlow," patients",sep="")),
           legend.title=""
           )

ggpar(p, 
      font.main = c(13, "bold"),
      font.x = c(16, "bold"),
      font.y = c(16, "bold"),
      font.caption = c(16, "bold"), 
      font.legend = c(13, "bold"), 
      font.tickslab = c(14, "bold"))


# RF
library("randomForestSRC")
# building RF model
B <- 1000
# Building a RSF
status_surv <- clin_filt$OS_STATUS=="1:DECEASED"
time_surv <- clin_filt$OS_MONTHS
RNA_filt <- RNA[clin$OS_MONTHS > 0,]

dataSetRF <- cbind(time_surv,status_surv,RNA_filt)

names(dataSetRF)<-make.names(names(dataSetRF))

RF_obj <- rfsrc(Surv(time_surv,status_surv)~., dataSetRF,  ntree = B,  membership = TRUE, importance=TRUE)
# Printing the RF object  
print(RF_obj)

# Vadiable importance
jk.obj <- subsample(RF_obj)
#pdf("VIMPsur.pdf", width = 15, height = 20)
#par(oma = c(0.5, 10, 0.5, 0.5))
#par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,17,1,1), mgp = c(4, 1, 0))
#plot(jk.obj, xlab = "Variable Importance (x 100)", cex = 1.2)
plot(jk.obj)
#dev.off()

fit.coxph <- coxph(surv_filt ~ TRIM67 + RASA1 + RHOF + NCAPD3 + NEU2 + ARNTL2 + CDK6, data = RNA_filt)
ggforest(fit.coxph, data = RNA_filt)
