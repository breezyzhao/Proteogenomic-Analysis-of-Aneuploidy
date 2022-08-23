rm(list=ls())

AS <- read.table("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC/CPTAC_pan_lm_Protein_AS.rnk", header=T)
colnames(AS)[2] <- paste0(colnames(AS)[2], "_AS")

AS_cancer <- read.table("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC/CPTAC_pan_lm_Protein_AS_cancer.rnk", header=T)
colnames(AS_cancer)[2] <- paste0(colnames(AS_cancer)[2], "_AS+Cancer")

AS_cnv <- read.table("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC/CPTAC_pan_lm_Protein_AS_cnv.rnk", header=T)
colnames(AS_cnv)[2] <- paste0(colnames(AS_cnv)[2], "_AS+CNV")

AS_estiPurity <- read.table("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC/CPTAC_pan_lm_Protein_AS_estiPurity.rnk", header=T)
colnames(AS_estiPurity)[2] <- paste0(colnames(AS_estiPurity)[2], "_AS+estiPurity")

AS_nucleiPurity <- read.table("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC/CPTAC_pan_lm_Protein_AS_nucleiPurity.rnk", header=T)
colnames(AS_nucleiPurity)[2] <- paste0(colnames(AS_nucleiPurity)[2], "_AS+nucleiPurity")

AS_CellCycle <- read.table("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC/CPTAC_pan_lm_Protein_AS_CellCycle.rnk", header=T)
colnames(AS_CellCycle)[2] <- paste0(colnames(AS_CellCycle)[2], "_AS+CellCycle")

AS_removeMito <- read.table("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC/CPTAC_pan_lm_Protein_AS_removeMito.rnk", header=T)
colnames(AS_removeMito)[2] <- paste0(colnames(AS_removeMito)[2], "_AS+removeMito")

finalTable <- merge(AS, AS_cancer, by.x="gene", by.y="gene", all=T)
finalTable <- merge(finalTable, AS_cnv, by.x="gene", by.y="gene", all=T)
finalTable <- merge(finalTable, AS_estiPurity, by.x="gene", by.y="gene", all=T)
finalTable <- merge(finalTable, AS_nucleiPurity, by.x="gene", by.y="gene", all=T)
finalTable <- merge(finalTable, AS_CellCycle, by.x="gene", by.y="gene", all=T)
finalTable <- merge(finalTable, AS_removeMito, by.x="gene", by.y="gene", all=T)

write.table(finalTable, file="tValue_differntModel_CPTAC_pan_cancer.txt", quote=F, sep="\t", row.names=F, col.names=T)
