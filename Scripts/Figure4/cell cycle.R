### this script is used to organize CPTAC cell cycle data
### organize the cell cycle scores calculated by raw RNA data
### validate by raw data

setwd("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis")

library(ggplot2)
library(dplyr)

rm(list=ls())

### colon2
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/colon2_addMutation.RData")
colon2 <- data.frame(samples=indepVar$Patient,
                     CellCycle=indepVar$CellCycle)
rm(cnv1, rna_partial_filter, rna_all_filter, protein1, indepVar)

### breast2
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/breast2_addMutation.RData")
breast2 <- data.frame(samples=indepVar$Patient,
                      CellCycle=indepVar$CellCycle)
rm(cnv1, rna_partial_filter, rna_all_filter, protein1, indepVar)

### ovarian2
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data/ovarian2.RData")
cc.signature <- c("CENPE","CCNA2","CCNB2","MCM6","CCNF","BUB1","CDC20","CDC6","CDK1","PLK1")
rna.cc <- rna1[cc.signature,]
rna.cc <- na.omit(rna.cc)
rna.cc <- data.frame(t(rna.cc))
rank.cc <- data.frame(apply(rna.cc, 2, rank))
ranked_colsum <- rank(apply(rank.cc, 1, sum))
rank.cc$CellCycle <- ranked_colsum
rank.cc$patient <- rownames(rank.cc)
ovarian2 <- data.frame(samples=rank.cc$patient,
                       CellCycle=rank.cc$CellCycle)
rm(cnv1, rna1, protein1, AS1, rank.cc, rna.cc)

### ccrcc
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/ccrcc_addMutation.RData")
ccrcc <- data.frame(samples=indepVar$Patient,
                    CellCycle=indepVar$CellCycle)
rm(cnv1, rna_partial_filter, rna_all_filter, protein1, indepVar)

### endometrial
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/endometrail_addMutation.RData")
endometrial <- data.frame(samples=indepVar$Patient,
                          CellCycle=indepVar$CellCycle)
rm(cnv1, rna_partial_filter, rna_all_filter, protein1, indepVar)

### hnscc
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/hnscc_addMutation.RData")
hnscc <- data.frame(samples=indepVar$Patient,
                    CellCycle=indepVar$CellCycle)
rm(cnv1, rna_partial_filter, rna_all_filter, protein1, indepVar)

### luad
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/luad_addMutation.RData")
luad <- data.frame(samples=indepVar$Patient,
                   CellCycle=indepVar$CellCycle)
rm(cnv1, rna_partial_filter, rna_all_filter, protein1, indepVar)

### put them together
CellCycle <- rbind(colon2, breast2, ovarian2, ccrcc, endometrial, hnscc, luad)

### load meta data
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis/CPTAC_purity.RData")
indep_pan3 <- merge(indep_pan2, CellCycle, by.x="patients", by.y="samples", all.x=T, sort=F)
indep_pan3$CellCycle <- as.numeric(indep_pan3$CellCycle)
save(indep_pan3, file = "CPTAC_purity_CellCycle.RData")

ggplot(indep_pan3, aes(x=cancer, y=CellCycle, color=cancer)) + geom_jitter()




