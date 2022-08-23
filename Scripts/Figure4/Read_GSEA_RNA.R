setwd("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC")

rm(list=ls())

library(ComplexHeatmap)
library(ggplot2)
library(qusage)
library(patchwork)

# ### GSEA results (same sets of genes)
CPTAC_AS_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Compensation/25_CPTAC_pan_commonGenes_lm_RNA_AS_C5BP/edb/results.edb"
CPTAC_AS_cancer_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Compensation/25_CPTAC_pan_commonGenes_lm_RNA_AS_cancer_C5BP/edb/results.edb"
CPTAC_AS_cnv_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Compensation/25_CPTAC_pan_commonGenes_lm_RNA_AS_cnv_C5BP/edb/results.edb"
CPTAC_AS_estiPurity_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Compensation/25_CPTAC_pan_commonGenes_lm_RNA_AS_estiPurity_C5BP/edb/results.edb"
CPTAC_AS_nucleiPurity_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Compensation/25_CPTAC_pan_commonGenes_lm_RNA_AS_nucleiPurity_C5BP/edb/results.edb"
CPTAC_AS_cellcycle_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Compensation/25_CPTAC_pan_commonGenes_lm_RNA_AS_CellCycle_C5BP/edb/results.edb"
CPTAC_AS_removeMito_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Compensation/25_CPTAC_pan_commonGenes_lm_RNA_AS_removeMito_C5BP/edb/results.edb"

CCLE_AS_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Compensation/26_CCLE_pan_commonGenes_lm_RNA_AS_C5BP/edb/results.edb"
CCLE_AS_cancer_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Compensation/26_CCLE_pan_commonGenes_lm_RNA_AS_cancer_C5BP/edb/results.edb"
CCLE_AS_cnv_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Compensation/26_CCLE_pan_commonGenes_lm_RNA_AS_cnv_C5BP/edb/results.edb"
CCLE_AS_cellcycle_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Compensation/26_CCLE_pan_commonGenes_lm_RNA_AS_CellCycle_C5BP/edb/results.edb"
CCLE_AS_removeMito_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Compensation/26_CCLE_pan_commonGenes_lm_RNA_AS_removeMito_C5BP/edb/results.edb"


KeyElement <- function(CharVec){
  keywords <- c("GENESET=","ES=","NES=","NP=","FDR=","FWER=")
  results <- unlist(lapply(keywords, function(x) CharVec[grep(x,CharVec)[1]]))
  return(results)
}

extractGeneSets <- function(EDBpathway){
  rawtable <- read.csv(EDBpathway, header=TRUE)
  rawtable <- data.frame(rawtable[2:(nrow(rawtable)-1),])
  selectedElement <- lapply(rawtable[,1], function(x) KeyElement(unlist(strsplit(as.character(x)," "))))
  selectedElement <- t(data.frame(selectedElement))
  colnames(selectedElement) <- c("GeneSet","ES","NES","NP","FDR","FWER")
  selectedElement <- data.frame(selectedElement)
  selectedElement$GeneSet <- as.character(sub("GENESET=.*#","",selectedElement[,"GeneSet"]))
  selectedElement$ES <- as.numeric(sub("ES=","",selectedElement[,"ES"]))
  selectedElement$NES <- as.numeric(sub("NES=","",selectedElement[,"NES"]))
  selectedElement$NP <- as.numeric(sub("NP=","",selectedElement[,"NP"]))
  selectedElement$FDR <- as.numeric(sub("FDR=","",selectedElement[,"FDR"]))
  selectedElement$FWER <- as.numeric(sub("FWER=","",selectedElement[,"FWER"]))
  rownames(selectedElement) <- selectedElement$GeneSet
  selectedElement <- selectedElement[order(selectedElement$NES,decreasing=TRUE),]
  return(selectedElement)
}

CPTAC_AS <- extractGeneSets(CPTAC_AS_path)
CPTAC_AS <- CPTAC_AS[,c(1,3,4,5)]
colnames(CPTAC_AS)[2:4] <- paste0(c("NES","NP","FDR"),"_CPTAC_AS")

CPTAC_AS_cancer <- extractGeneSets(CPTAC_AS_cancer_path)
CPTAC_AS_cancer <- CPTAC_AS_cancer[,c(1,3,4,5)]
colnames(CPTAC_AS_cancer)[2:4] <- paste0(c("NES","NP","FDR"),"_CPTAC_AS_cancer")

CPTAC_AS_cnv <- extractGeneSets(CPTAC_AS_cnv_path)
CPTAC_AS_cnv <- CPTAC_AS_cnv[,c(1,3,4,5)]
colnames(CPTAC_AS_cnv)[2:4] <- paste0(c("NES","NP","FDR"),"_CPTAC_AS_cnv")

CPTAC_AS_estiPurity <- extractGeneSets(CPTAC_AS_estiPurity_path)
CPTAC_AS_estiPurity <- CPTAC_AS_estiPurity[,c(1,3,4,5)]
colnames(CPTAC_AS_estiPurity)[2:4] <- paste0(c("NES","NP","FDR"),"_CPTAC_AS_estiPurity")

CPTAC_AS_nucleiPurity <- extractGeneSets(CPTAC_AS_nucleiPurity_path)
CPTAC_AS_nucleiPurity <- CPTAC_AS_nucleiPurity[,c(1,3,4,5)]
colnames(CPTAC_AS_nucleiPurity)[2:4] <- paste0(c("NES","NP","FDR"),"_CPTAC_AS_nucleiPurity")

CPTAC_AS_cellcycle <- extractGeneSets(CPTAC_AS_cellcycle_path)
CPTAC_AS_cellcycle <- CPTAC_AS_cellcycle[,c(1,3,4,5)]
colnames(CPTAC_AS_cellcycle)[2:4] <- paste0(c("NES","NP","FDR"),"_CPTAC_AS_cellcycle")

CPTAC_AS_removeMito <- extractGeneSets(CPTAC_AS_removeMito_path)
CPTAC_AS_removeMito <- CPTAC_AS_removeMito[,c(1,3,4,5)]
colnames(CPTAC_AS_removeMito)[2:4] <- paste0(c("NES","NP","FDR"),"_CPTAC_AS_removeMito")

CCLE_AS <- extractGeneSets(CCLE_AS_path)
CCLE_AS <- CCLE_AS[,c(1,3,4,5)]
colnames(CCLE_AS)[2:4] <- paste0(c("NES","NP","FDR"),"_CCLE_AS")

CCLE_AS_cancer <- extractGeneSets(CCLE_AS_cancer_path)
CCLE_AS_cancer <- CCLE_AS_cancer[,c(1,3,4,5)]
colnames(CCLE_AS_cancer)[2:4] <- paste0(c("NES","NP","FDR"),"_CCLE_AS_cancer")

CCLE_AS_cnv <- extractGeneSets(CCLE_AS_cnv_path)
CCLE_AS_cnv <- CCLE_AS_cnv[,c(1,3,4,5)]
colnames(CCLE_AS_cnv)[2:4] <- paste0(c("NES","NP","FDR"),"_CCLE_AS_cnv")

CCLE_AS_cellcycle <- extractGeneSets(CCLE_AS_cellcycle_path)
CCLE_AS_cellcycle <- CCLE_AS_cellcycle[,c(1,3,4,5)]
colnames(CCLE_AS_cellcycle)[2:4] <- paste0(c("NES","NP","FDR"),"_CCLE_AS_cellcycle")

CCLE_AS_removeMito <- extractGeneSets(CCLE_AS_removeMito_path)
CCLE_AS_removeMito <- CCLE_AS_removeMito[,c(1,3,4,5)]
colnames(CCLE_AS_removeMito)[2:4] <- paste0(c("NES","NP","FDR"),"_CCLE_AS_removeMito")

### combine GSEA results
GSEA_Protein <- merge(CPTAC_AS, CPTAC_AS_cancer, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, CPTAC_AS_cnv, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, CPTAC_AS_estiPurity, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, CPTAC_AS_nucleiPurity, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, CPTAC_AS_cellcycle, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, CPTAC_AS_removeMito, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, CCLE_AS, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, CCLE_AS_cancer, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, CCLE_AS_cnv, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, CCLE_AS_cellcycle, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, CCLE_AS_removeMito, by.x="GeneSet",by.y="GeneSet",all=T)
write.table(GSEA_Protein, file="CCLE_CPTAC_pan_commonGenes_RNA_C5.txt", quote=FALSE, sep="\t", row.names=F, col.names=T, na="")

GSEA_filter1 <- na.omit(GSEA_Protein)
rownames(GSEA_filter1) <- GSEA_filter1$GeneSet
GSEA_filter1 <- GSEA_filter1[,-1]
p_threshold <- 0.01
GSEA_filter1_p <- GSEA_filter1[seq(2, ncol(GSEA_filter1),3)]
GSEA_filter1_NES <- GSEA_filter1[seq(1, ncol(GSEA_filter1),3)]
GSEA_filter1_NES_threshold <- GSEA_filter1_NES
GSEA_filter1_NES_threshold[GSEA_filter1_p>p_threshold] <- "NS"

genesets <- c("GO_DNA_TEMPLATED_TRANSCRIPTION_ELONGATION", "GO_DNA_TEMPLATED_TRANSCRIPTION_INITIATION",
              "GO_DNA_TEMPLATED_TRANSCRIPTION_TERMINATION", "GO_REGULATION_OF_MRNA_SPLICING_VIA_SPLICEOSOME",
              "GO_RNA_3_END_PROCESSING", "GO_RNA_5_END_PROCESSING",
              "GO_MRNA_TRANSPORT", "GO_RNA_LOCALIZATION",
              "GO_MRNA_MODIFICATION", "GO_RNA_CATABOLIC_PROCESS",
              "GO_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS","GO_TRNA_METABOLIC_PROCESS",
              "GO_TRNA_PROCESSING", "GO_RRNA_METABOLIC_PROCESS",
              "GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION", "GO_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_ORGANELLE",
              "GO_PROTEIN_LOCALIZATION_TO_CELL_PERIPHERY","GO_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
              "GO_REGULATION_OF_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS","GO_MEMBRANE_PROTEIN_PROTEOLYSIS")

GSEA_filter1_NES_threshold_pick <- GSEA_filter1_NES_threshold[genesets,]
save(GSEA_filter1_NES_threshold_pick, file="CPTAC_CCLE_GSEA_NES_threshold_LM_RNA_AS.RData")



