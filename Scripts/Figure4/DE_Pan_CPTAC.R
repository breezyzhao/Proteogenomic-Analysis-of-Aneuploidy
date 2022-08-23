### this script is to use lm to analyze log2FC data (pan-CPTAC)

setwd("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC")

library(ggplot2)
library(dplyr)

rm(list=ls())

load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC/pan_omics_log2FC.RData")
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC/pan_aneuploidy.RData")
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC/CPTAC_purity.RData")
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC/CPTAC_purity_CellCycle.RData")
rownames(indep_pan2) <- indep_pan2$patients
rownames(indep_pan3) <- indep_pan2$patients
if (length(intersect(indep_pan$patients, indep_pan2$patients))!=nrow(indep_pan)) {stop("something is wrong!")}
indep_pan2 <- indep_pan2[indep_pan$patients,]
if (length(intersect(indep_pan$patients, indep_pan3$patients))!=nrow(indep_pan)) {stop("something is wrong!")}
indep_pan3 <- indep_pan3[indep_pan$patients,]

### check data structure
if (sum(indep_pan2$patients!=colnames(cnv_pan) | indep_pan2$patients!=colnames(rna_pan) | indep_pan2$patients!=colnames(protein_pan))!=0) {stop("check patient ordering!")}
if (sum(indep_pan3$patients!=colnames(cnv_pan) | indep_pan3$patients!=colnames(rna_pan) | indep_pan3$patients!=colnames(protein_pan))!=0) {stop("check patient ordering!")}
if (sum(indep_pan$patients!=colnames(cnv_pan) | indep_pan$patients!=colnames(rna_pan) | indep_pan$patients!=colnames(protein_pan))!=0) {stop("check patient ordering!")}
if (sum(rownames(cnv_pan)!=rownames(rna_pan) | rownames(cnv_pan)!=rownames(protein_pan))!=0) {stop("check gene ordering!")}

### check aneuploidy score
ggplot(indep_pan, aes(x=cancer, y=aneuploidy, fill=cancer)) + geom_dotplot(binaxis='y', stackdir='center', stackratio=1, binwidth=0.3)

# gene_CPTAC <- rownames(cnv_pan)
# save(gene_CPTAC, file="gene_CPTAC.RData")

# load("gene_CPTAC.RData")
# load("gene_CCLE.RData")
# commonGene <- intersect(gene_CCLE, gene_CPTAC)
# cnv_pan <- cnv_pan[commonGene,]
# rna_pan <- rna_pan[commonGene,]
# protein_pan <- protein_pan[commonGene,]

### we have 4328 genes of 682 patients
### first try: only use aneuploidy score to fit protein log2FC
# indepVar <- indep_pan2
# depVar <- protein_pan[1,]
lm_self <- function(depVar, indepVar) {
  n <- sum(!is.na(depVar))
  if (n<0.2*ncol(protein_pan)) {
    result <- c("n_noNA"=n,
                "AS_coeff"=NA,
                "AS_Tvalue"=NA,
                "AS_Pvalue"=NA)
  } else {
    lm_data <- indepVar
    lm_data$Expression <- as.numeric(t(depVar))
    # lm_model<-lm(Expression~aneuploidy,data=lm_data,na.action=na.omit)
    lm_model<-lm(Expression~aneuploidy+cancer,data=lm_data,na.action=na.omit)
    # lm_model<-lm(Expression~aneuploidy+estiPurity,data=lm_data,na.action=na.omit)
    # lm_model<-lm(Expression~aneuploidy+nucleiPurity,data=lm_data,na.action=na.omit)
    # lm_model<-lm(Expression~aneuploidy+CellCycle,data=lm_data,na.action=na.omit)
    if (sum(is.na(lm_model$coefficients))==0) {
      result <- c("n_noNA"=n,
                  "AS_coeff"=summary(lm_model)$coefficients[2,1],
                  "AS_Tvalue"=summary(lm_model)$coefficients[2,3],
                  "AS_Pvalue"=summary(lm_model)$coefficients[2,4])
    } else {
      result <- c("n_noNA"=n,
                  "AS_coeff"=NA,
                  "AS_Tvalue"=NA,
                  "AS_Pvalue"=NA)
    }
  }
}

lm_result <- apply(protein_pan,1,function(x) lm_self(x,indep_pan3))
lm_result <- as.data.frame(t(lm_result))
lm_result <- lm_result[!is.na(lm_result$AS_Tvalue),]
lm_result$AS_adjP <- p.adjust(lm_result$AS_Pvalue, "BH")
lm_result <- lm_result[lm_result$n_noNA>0.2*ncol(protein_pan),]
lm_result <- lm_result[rownames(lm_result)!="",]
lm_result <- lm_result[order(lm_result$AS_Tvalue,decreasing=TRUE),]
lm_result_export <- lm_result
lm_result_export$gene <- rownames(lm_result_export)
lm_result_export <- lm_result_export[,c(6,3)]

# lm_rna <- lm_result
lm_protein <- lm_result
lm_CPTAC_pan_expression_aneuploidy_cancerType <- merge(lm_rna, lm_protein, by.x="row.names", by.y="row.names", all=T)
write.table(lm_CPTAC_pan_expression_aneuploidy_cancerType, "lm_CPTAC_pan_expression_aneuploidy_cancerType.txt",row.names=F, col.names=T, sep="\t", quote=FALSE)
save(lm_CPTAC_pan_expression_aneuploidy_cancerType, file="lm_CPTAC_pan_expression_aneuploidy_cancerType.RData")

### (optional) remove mito genes
# mito <- read.csv("/Users/pc2644/Documents/DM_Aneuploidy/mitoChromosome/Human.MitoCarta2.0.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# lm_result_export <- lm_result_export[!(lm_result_export$gene %in% mito$Symbol),]
# 
# write.table(lm_result_export, "CPTAC_pan_commonGenes_lm_RNA_AS_removeMito.rnk",row.names=F, sep="\t", quote=FALSE)

# lm_result$gene <- rownames(lm_result)
# write.table(lm_result, "CPTAC_pan_commonGenes_lm_Protein_AS_statistic.rnk",row.names=F, sep="\t", quote=FALSE)

### add cnv into the model
rna_cnv <- cbind(rna_pan,cnv_pan)
protein_cnv <- cbind(protein_pan,cnv_pan)
# depVar <- protein_cnv[1,]
# indepVar <- indep_pan
lm_self2 <- function(depVar, indepVar) {
  n <- sum(!is.na(depVar[1:ncol(protein_pan)]))
  if (n<0.2*ncol(protein_pan)) {
    result <- c("n_noNA"=n,
                "AS_coeff"=NA,
                "AS_Tvalue"=NA,
                "AS_Pvalue"=NA)
  } else {
    lm_data <- indepVar
    lm_data$Expression <- as.numeric(t(depVar[1:ncol(protein_pan)]))
    lm_data$cnv <- as.numeric(t(depVar[(ncol(protein_pan)+1):ncol(protein_cnv)]))
    # lm_model<-lm(Expression~aneuploidy,data=lm_data,na.action=na.omit)
    # lm_model<-lm(Expression~aneuploidy+cancer,data=lm_data,na.action=na.omit)
    # lm_model<-lm(Expression~aneuploidy+cancer+cnv,data=lm_data,na.action=na.omit)
    lm_model<-lm(Expression~aneuploidy+cnv,data=lm_data,na.action=na.omit)
    if (sum(is.na(lm_model$coefficients))==0) {
      result <- c("n_noNA"=n,
                  "AS_coeff"=summary(lm_model)$coefficients[2,1],
                  "AS_Tvalue"=summary(lm_model)$coefficients[2,3],
                  "AS_Pvalue"=summary(lm_model)$coefficients[2,4])
    } else {
      result <- c("n_noNA"=n,
                  "AS_coeff"=NA,
                  "AS_Tvalue"=NA,
                  "AS_Pvalue"=NA)
    }
  }
  return(result)
}

lm_result <- apply(rna_cnv,1,function(x) lm_self2(x,indep_pan))
lm_result <- as.data.frame(t(lm_result))
lm_result <- lm_result[!is.na(lm_result$AS_Tvalue),]
lm_result$AS_adjP <- p.adjust(lm_result$AS_Pvalue, "BH")
lm_result <- lm_result[lm_result$n_noNA>0.2*ncol(protein_pan),]
lm_result <- lm_result[rownames(lm_result)!="",]
lm_result <- lm_result[order(lm_result$AS_Tvalue,decreasing=TRUE),]
lm_result_export <- lm_result
lm_result_export$gene <- rownames(lm_result_export)
lm_result_export <- lm_result_export[,c(6,3)]
write.table(lm_result_export, "CPTAC_pan_commonGenes_lm_RNA_AS_cnv.rnk",row.names=F, sep="\t", quote=FALSE)



