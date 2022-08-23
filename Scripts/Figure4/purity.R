### this script is used to organize CPTAC purity data
### two kinds of purity data
### first: estimated from RNA
### second: tumor nuclei from HE stain

setwd("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis")

library(ggplot2)
library(dplyr)

rm(list=ls())

### first organize estimate purity data
### colon2
path_colon <- "/Volumes/davolt01lab/davolt01labspace/Data/CPTAC_PanCan/CPTAC_PanCan_Clinical_Data/CO/AWG_data_freeze/Human__CPTAC_COAD__MS__Clinical__Clinical__03_01_2017__CPTAC__Clinical__BCM.tsi"
sample <- read.delim(path_colon)
sample <- data.frame(t(sample))
colnames(sample) <- sample[1,]
sample <- sample[-1,]
purity_colon <- data.frame(samples=rownames(sample),
                           estiPurity=sample$TumorPurity)
rm(sample, path_colon)

### breast2
path_breast <- "/Volumes/davolt01lab/davolt01labspace/Data/CPTAC_PanCan/CPTAC_PanCan_Clinical_Data/BR/AWG_data_freeze/prosp-brca-v5.3-sample-annotation.csv"
sample <- read.csv(path_breast)
purity_breast <- data.frame(samples=sample$Sample.ID,
                            estiPurity=sample$ESTIMATE.TumorPurity)
rm(sample, path_breast)

### ovarian2
# path_ovarian <- "/Volumes/davolt01lab/davolt01labspace/Data/CPTAC_PanCan/CPTAC_PanCan_Clinical_Data/OV/AWG_data_freeze/Ovary One Year Clinical Data_20160927.txt"
# sample <- read.delim(path_ovarian)
# rm(sample, path_ovarian)

### ccrcc
# path_ccrcc <- "/Volumes/davolt01lab/davolt01labspace/Data/CPTAC_PanCan/CPTAC_PanCan_Clinical_Data/ccRCC/AWG_data_freeze/CCRCC_September2018_case.tsv"
# sample <- read.delim(path_ccrcc)
# rm(sample, path_ccrcc)

### endometrial
path_endometrial <- "/Volumes/davolt01lab/davolt01labspace/Data/CPTAC_PanCan/CPTAC_PanCan_Clinical_Data/EC/AWG_data_freeze/HS_CPTAC_UCEC_CLI.txt"
sample <- read.delim(path_endometrial)
purity_endometrial <- data.frame(samples=gsub("-",".",sample$Proteomics_Participant_ID),
                            estiPurity=sample$Purity_Cancer)
rm(sample, path_endometrial)

### hnscc
path_hnscc <- "/Volumes/davolt01lab/davolt01labspace/Data/CPTAC_PanCan/CPTAC_PanCan_Clinical_Data/HNSCC/AWG_data_freeze/Meta_table.tsv"
sample <- read.delim(path_hnscc)
purity_hnscc <- data.frame(samples=gsub("-",".",sample$case_id),
                           estiPurity=sample$tumor_proportion)

rm(sample, path_hnscc)

### luad
path_luad <- "/Volumes/davolt01lab/davolt01labspace/Data/CPTAC_PanCan/CPTAC_PanCan_Clinical_Data/LUAD/AWG_data_freeze/luad-v3.2-sample-annotation.csv"
sample <- read.csv(path_luad)
purity_luad <- data.frame(samples=sample$Sample.ID,
                          estiPurity=sample$Tumor.Purity.byESTIMATE.RNAseq)
rm(sample, path_luad)

purity_estimate <- rbind(purity_colon, purity_breast, purity_endometrial, purity_hnscc, purity_luad)
purity_estimate <- na.omit(purity_estimate)

### read aneuploidy score
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis/pan_aneuploidy.RData")

indep_pan2 <- merge(indep_pan, purity_estimate, by.x="patients", by.y="samples", all.x=T, sort=F)

### read purity data from image
path_cohort <- "/Volumes/davolt01lab/davolt01labspace/Data/Proteogenomic_Fenyo/CPTAC3/CPTAC3_cohort.csv"
sample <- read.csv(path_cohort)
sample <- sample[sample$Specimen_Type=="tumor_tissue",]
### there are some tissues from the same patients, calculate the average of them
patients <- unique(sample$Patient_ID)
purity_nuclei <- data.frame(samples=patients,
                            nucleiPurity=rep(NA,length(patients)))
rownames(purity_nuclei) <- purity_nuclei$samples
# i=patients[1]
for (i in patients) {
  index <- sample$Patient_ID==i
  purity_nuclei[i,"nucleiPurity"] <- mean(sample$Percent_Tumor_Nuclei[index])
  
}
purity_nuclei$samples <- gsub("-", "." , purity_nuclei$samples)
rm(sample, path_cohort)

indep_pan2 <- merge(indep_pan2, purity_nuclei, by.x="patients", by.y="samples", all.x=T, sort=F)
indep_pan2$estiPurity <- as.numeric(indep_pan2$estiPurity)
indep_pan2$nucleiPurity <- as.numeric(indep_pan2$nucleiPurity)
save(indep_pan2, file = "CPTAC_purity.RData")

indep_pan2_plot <- na.omit(indep_pan2)
ggplot(indep_pan2_plot, aes(x=estiPurity, y=nucleiPurity, color=cancer)) + 
  geom_point() +
  labs(x="purity from genomic data", y="purity from nuclei")

ggplot(indep_pan2, aes(x=cancer, y=estiPurity, color=cancer)) + geom_jitter()

ggplot(indep_pan2, aes(x=cancer, y=nucleiPurity, color=cancer)) + geom_jitter()
