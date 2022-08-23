library(ggplot2)
library(ggfortify)
library(ggrepel)
rm(list=ls())

### colon2 ###
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/colon2_addMutation.RData")
if (class(cnv1[,1])!="numeric") {
  cnv0 <- data.frame(apply(cnv1,2,as.numeric))
  rownames(cnv0) <- rownames(cnv1)
  cnv1 <- cnv0
  rm(cnv0)
}
if (class(rna_all_filter[,1])!="numeric") {
  rna0 <- data.frame(apply(rna_all_filter,2,as.numeric))
  rownames(rna0) <- rownames(rna_all_filter)
  rna1 <- rna0
  rm(rna0)
} else {rna1 <- rna_all_filter}
if (class(protein1[,1])!="numeric") {
  protein0 <- data.frame(apply(protein1,2,as.numeric))
  rownames(protein0) <- rownames(protein1)
  protein1 <- protein0
  rm(protein0)
}
cnv_colon2 <- cnv1
rna_colon2 <- data.frame(scale(rna1, center=T, scale=F))
protein_colon2 <- protein1
indep_colon2 <- indepVar

rm(list=c("cnv1", "rna_all_filter", "rna_partial_filter", "rna1", "protein1", "indepVar"))
### colon2 ###

### breast2 ###
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/breast2_addMutation.RData")
if (class(cnv1[,1])!="numeric") {
  cnv0 <- data.frame(apply(cnv1,2,as.numeric))
  rownames(cnv0) <- rownames(cnv1)
  cnv1 <- cnv0
  rm(cnv0)
}
if (class(rna_all_filter[,1])!="numeric") {
  rna0 <- data.frame(apply(rna_all_filter,2,as.numeric))
  rownames(rna0) <- rownames(rna_all_filter)
  rna1 <- rna0
  rm(rna0)
} else {rna1 <- rna_all_filter}
if (class(protein1[,1])!="numeric") {
  protein0 <- data.frame(apply(protein1,2,as.numeric))
  rownames(protein0) <- rownames(protein1)
  protein1 <- protein0
  rm(protein0)
}
cnv_breast2 <- cnv1
rna_breast2 <- rna1
protein_breast2 <- protein1
indep_breast2 <- indepVar

rm(list=c("cnv1", "rna_all_filter", "rna_partial_filter", "rna1", "protein1", "indepVar"))
### breast2 ###

### ovarian2 ###
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data/ovarian2.RData")
if (class(cnv1[,1])!="numeric") {
  cnv0 <- data.frame(apply(cnv1,2,as.numeric))
  rownames(cnv0) <- rownames(cnv1)
  cnv1 <- cnv0
  rm(cnv0)
}
if (class(rna1[,1])!="numeric") {
  rna0 <- data.frame(apply(rna1,2,as.numeric))
  rownames(rna0) <- rownames(rna1)
  rna1 <- rna0
  rm(rna0)
}
if (class(protein1[,1])!="numeric") {
  protein0 <- data.frame(apply(protein1,2,as.numeric))
  rownames(protein0) <- rownames(protein1)
  protein1 <- protein0
  rm(protein0)
}
cnv_ovarian2 <- cnv1
rna_ovarian2 <- log2(rna1+1)
protein_ovarian2 <- protein1
indep_ovarian2 <- AS1

rm(list=c("cnv1", "rna1", "protein1", "AS1"))
### ovarian2 ###

### ccrcc ###
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/ccrcc_addMutation.RData")
if (class(cnv1[,1])!="numeric") {
  cnv0 <- data.frame(apply(cnv1,2,as.numeric))
  rownames(cnv0) <- rownames(cnv1)
  cnv1 <- cnv0
  rm(cnv0)
}
if (class(rna_all_filter[,1])!="numeric") {
  rna0 <- data.frame(apply(rna_all_filter,2,as.numeric))
  rownames(rna0) <- rownames(rna_all_filter)
  rna1 <- rna0
  rm(rna0)
} else {rna1 <- rna_all_filter}
if (class(protein1[,1])!="numeric") {
  protein0 <- data.frame(apply(protein1,2,as.numeric))
  rownames(protein0) <- rownames(protein1)
  protein1 <- protein0
  rm(protein0)
}
cnv_ccrcc <- cnv1
rna_ccrcc <- rna1
protein_ccrcc <- protein1
indep_ccrcc <- indepVar

rm(list=c("cnv1", "rna_all_filter", "rna_partial_filter", "rna1", "protein1", "indepVar"))
### ccrcc ###

### endometrial ###
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/endometrail_addMutation.RData")
if (class(cnv1[,1])!="numeric") {
  cnv0 <- data.frame(apply(cnv1,2,as.numeric))
  rownames(cnv0) <- rownames(cnv1)
  cnv1 <- cnv0
  rm(cnv0)
}
if (class(rna_all_filter[,1])!="numeric") {
  rna0 <- data.frame(apply(rna_all_filter,2,as.numeric))
  rownames(rna0) <- rownames(rna_all_filter)
  rna1 <- rna0
  rm(rna0)
} else {rna1 <- rna_all_filter}
if (class(protein1[,1])!="numeric") {
  protein0 <- data.frame(apply(protein1,2,as.numeric))
  rownames(protein0) <- rownames(protein1)
  protein1 <- protein0
  rm(protein0)
}
if ((sum(indepVar$Patient!=colnames(cnv1)) + sum(indepVar$Patient!=colnames(rna_all_filter)) + sum(indepVar$Patient!=colnames(protein1)))!=0) {stop("check the ordering of samples")}
index <- !grepl("C3N.01003", colnames(cnv1))
cnv_endometrial <- cnv1[,index]
rna_endometrial <- rna1[,index]
protein_endometrial <- protein1[,index]
indep_endometrial <- indepVar[index,]

rm(list=c("cnv1", "rna_all_filter", "rna_partial_filter", "rna1", "protein1", "indepVar"))
### endometrail ###

### hnscc ###
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/hnscc_addMutation.RData")
if (class(cnv1[,1])!="numeric") {
  cnv0 <- data.frame(apply(cnv1,2,as.numeric))
  rownames(cnv0) <- rownames(cnv1)
  cnv1 <- cnv0
  rm(cnv0)
}
if (class(rna_all_filter[,1])!="numeric") {
  rna0 <- data.frame(apply(rna_all_filter,2,as.numeric))
  rownames(rna0) <- rownames(rna_all_filter)
  rna1 <- rna0
  rm(rna0)
} else {rna1 <- rna_all_filter}
if (class(protein1[,1])!="numeric") {
  # protein0 <- data.frame(apply(protein1,2,as.numeric))
  protein0 <- protein1
  rownames(protein0) <- rownames(protein1)
  protein1 <- protein0
  rm(protein0)
}
cnv_hnscc <- cnv1
rna_hnscc <- rna1
protein_hnscc <- data.frame(scale(protein1, center=T, scale=F))
indep_hnscc <- indepVar

rm(list=c("cnv1", "rna_all_filter", "rna_partial_filter", "rna1", "protein1", "indepVar"))
### hnscc ###

### LUAD ###
load("/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/luad_addMutation.RData")
if (class(cnv1[,1])!="numeric") {
  cnv0 <- data.frame(apply(cnv1,2,as.numeric))
  rownames(cnv0) <- rownames(cnv1)
  cnv1 <- cnv0
  rm(cnv0)
}
if (class(rna_all_filter[,1])!="numeric") {
  rna0 <- data.frame(apply(rna_all_filter,2,as.numeric))
  rownames(rna0) <- rownames(rna_all_filter)
  rna1 <- rna0
  rm(rna0)
} else {rna1 <- rna_all_filter}
if (class(protein1[,1])!="numeric") {
  protein0 <- data.frame(apply(protein1,2,as.numeric))
  rownames(protein0) <- rownames(protein1)
  protein1 <- protein0
  rm(protein0)
}
if ((sum(indepVar$Patient!=colnames(cnv1)) + sum(indepVar$Patient!=colnames(rna_all_filter)) + sum(indepVar$Patient!=colnames(protein1)))!=0) {stop("check the ordering of samples")}
index <- !grepl("C3L.00510", colnames(cnv1))
cnv_luad <- cnv1[,index]
rna_luad <- rna1[,index]
protein_luad <- protein1[,index]
indep_luad <- indepVar[index,]

rm(list=c("cnv1", "rna_all_filter", "rna_partial_filter", "rna1", "protein1", "indepVar"))
### luad ###

### check the dirstribution of CNV
# cnv_plot_colon2 <- data.frame(cnv=as.vector(as.matrix(cnv_colon2)), cancer="colon2")
# cnv_plot_breast2 <- data.frame(cnv=as.vector(as.matrix(cnv_breast2)), cancer="breast2")
# cnv_plot_ovarian2 <- data.frame(cnv=as.vector(as.matrix(cnv_ovarian2)), cancer="ovarian2")
# cnv_plot_ccrcc <- data.frame(cnv=as.vector(as.matrix(cnv_ccrcc)), cancer="ccrcc")
# cnv_plot_endometrial <- data.frame(cnv=as.vector(as.matrix(cnv_endometrial)), cancer="endometrial")
# cnv_plot_hnscc <- data.frame(cnv=as.vector(as.matrix(cnv_hnscc)), cancer="hnscc")
# cnv_plot_luad <- data.frame(cnv=as.vector(as.matrix(cnv_hnscc)), cancer="luad")
# cnv_plot <- rbind(cnv_plot_colon2, cnv_plot_breast2, cnv_plot_ovarian2, cnv_plot_ccrcc, cnv_plot_endometrial, cnv_plot_hnscc, cnv_plot_luad)
# 
# ggplot(cnv_plot, aes(cnv, stat(density))) +
#   geom_density(alpha=0.1, size = 1) +
#   xlab("CNV") +
#   facet_wrap(~cancer) +
#   coord_cartesian(xlim=c(-5,5))

### conbine cnv data
cnv_pan <- merge(cnv_colon2, cnv_breast2, by.x="row.names", by.y="row.names")
cnv_pan <- merge(cnv_pan, cnv_ovarian2, by.x="Row.names", by.y="row.names")
cnv_pan <- merge(cnv_pan, cnv_ccrcc, by.x="Row.names", by.y="row.names")
cnv_pan <- merge(cnv_pan, cnv_endometrial, by.x="Row.names", by.y="row.names")
cnv_pan <- merge(cnv_pan, cnv_hnscc, by.x="Row.names", by.y="row.names")
cnv_pan <- merge(cnv_pan, cnv_luad, by.x="Row.names", by.y="row.names")
rownames(cnv_pan) <- cnv_pan$Row.names
cnv_pan <- cnv_pan[,-1]

# ### check the dirstribution of RNA
# rna_plot_colon2 <- data.frame(rna=as.vector(as.matrix(rna_colon2)), cancer="colon2")
# rna_plot_breast2 <- data.frame(rna=as.vector(as.matrix(rna_breast2)), cancer="breast2")
# rna_plot_ovarian2 <- data.frame(rna=as.vector(as.matrix(rna_ovarian2)), cancer="ovarian2")
# rna_plot_ccrcc <- data.frame(rna=as.vector(as.matrix(rna_ccrcc)), cancer="ccrcc")
# rna_plot_endometrial <- data.frame(rna=as.vector(as.matrix(rna_endometrial)), cancer="endometrial")
# rna_plot_hnscc <- data.frame(rna=as.vector(as.matrix(rna_hnscc)), cancer="hnscc")
# rna_plot_luad <- data.frame(rna=as.vector(as.matrix(rna_hnscc)), cancer="luad")
# rna_plot <- rbind(rna_plot_colon2, rna_plot_breast2, rna_plot_ovarian2, rna_plot_ccrcc, rna_plot_endometrial, rna_plot_hnscc, rna_plot_luad)
# 
# ggplot(rna_plot, aes(rna, stat(density))) +
#   geom_density(alpha=0.1, size = 1) +
#   xlab("RNA") +
#   coord_cartesian(xlim=c(-5,5)) +
#   facet_wrap(~cancer)

### check the dirstribution of protein
# protein_plot_colon2 <- data.frame(protein=as.vector(as.matrix(protein_colon2)), cancer="colon2")
# protein_plot_breast2 <- data.frame(protein=as.vector(as.matrix(protein_breast2)), cancer="breast2")
# protein_plot_ovarian2 <- data.frame(protein=as.vector(as.matrix(protein_ovarian2)), cancer="ovarian2")
# protein_plot_ccrcc <- data.frame(protein=as.vector(as.matrix(protein_ccrcc)), cancer="ccrcc")
# protein_plot_endometrial <- data.frame(protein=as.vector(as.matrix(protein_endometrial)), cancer="endometrial")
# protein_plot_hnscc <- data.frame(protein=as.vector(as.matrix(protein_hnscc)), cancer="hnscc")
# protein_plot_luad <- data.frame(protein=as.vector(as.matrix(protein_hnscc)), cancer="luad")
# protein_plot <- rbind(protein_plot_colon2, protein_plot_breast2, protein_plot_ovarian2, protein_plot_ccrcc, protein_plot_endometrial, protein_plot_hnscc, protein_plot_luad)
# 
# ggplot(protein_plot, aes(protein, stat(density))) +
#   geom_density(alpha=0.1, size = 1) +
#   xlab("Protein") +
#   coord_cartesian(xlim=c(-10,10)) +
#   facet_wrap(~cancer)

### conbine protein data
protein_pan <- merge(protein_colon2, protein_breast2, by.x="row.names", by.y="row.names", all=T)
protein_pan <- merge(protein_pan, protein_ovarian2, by.x="Row.names", by.y="row.names", all=T)
protein_pan <- merge(protein_pan, protein_ccrcc, by.x="Row.names", by.y="row.names", all=T)
protein_pan <- merge(protein_pan, protein_endometrial, by.x="Row.names", by.y="row.names", all=T)
protein_pan <- merge(protein_pan, protein_hnscc, by.x="Row.names", by.y="row.names", all=T)
protein_pan <- merge(protein_pan, protein_luad, by.x="Row.names", by.y="row.names", all=T)
# protein_pan <- merge(protein_colon2, protein_breast2, by.x="row.names", by.y="row.names")
# protein_pan <- merge(protein_pan, protein_ovarian2, by.x="Row.names", by.y="row.names")
# protein_pan <- merge(protein_pan, protein_ccrcc, by.x="Row.names", by.y="row.names")
# protein_pan <- merge(protein_pan, protein_endometrial, by.x="Row.names", by.y="row.names")
# protein_pan <- merge(protein_pan, protein_hnscc, by.x="Row.names", by.y="row.names")
# protein_pan <- merge(protein_pan, protein_luad, by.x="Row.names", by.y="row.names")
rownames(protein_pan) <- protein_pan$Row.names
protein_pan <- protein_pan[,-1]
protein_pan <- data.frame(scale(protein_pan, center=T, scale=F))

### combine aneuploidy score
indep_colon2_selected <- data.frame(cancer="colon2", patients=indep_colon2$Patient, aneuploidy=indep_colon2$AS)
indep_breast2_selected <- data.frame(cancer="breast2", patients=indep_breast2$Patient, aneuploidy=indep_breast2$AS)
indep_ovarian2_selected <- data.frame(cancer="ovarian2", patients=indep_ovarian2$sample, aneuploidy=indep_ovarian2$AS.0.2.noX)
indep_ccrcc_selected <- data.frame(cancer="ccrcc", patients=indep_ccrcc$Patient, aneuploidy=indep_ccrcc$AS)
indep_endometrial_selected <- data.frame(cancer="endometrial", patients=indep_endometrial$Patient, aneuploidy=indep_endometrial$AS)
indep_hnscc_selected <- data.frame(cancer="hnscc", patients=indep_hnscc$Patient, aneuploidy=indep_hnscc$AS)
indep_luad_selected <- data.frame(cancer="luad", patients=indep_luad$Patient, aneuploidy=indep_luad$AS)
indep_pan <- rbind(indep_colon2_selected, indep_breast2_selected, indep_ovarian2_selected, indep_ccrcc_selected, indep_endometrial_selected, indep_hnscc_selected, indep_luad_selected)

### calculate the median and variation of each sample
if ((sum(indep_pan$patients!=colnames(cnv_pan)) + sum(indep_pan$patients!=colnames(protein_pan)))!=0) {stop("check the ordering of samples")}

save(indep_pan, file = "/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/pan_aneuploidy.RData")

summary_pan <- data.frame(patients=indep_pan$patients,
                          cancer=factor(indep_pan$cancer, levels=unique(indep_pan$cancer)),
                          cnv_median=apply(cnv_pan,2,median,na.rm=T),
                          cnv_sd=apply(cnv_pan,2,sd,na.rm=T),
                          protein_median=apply(protein_pan,2,median,na.rm=T),
                          protein_sd=apply(protein_pan,2,sd,na.rm=T)
                          )
ggplot(summary_pan, aes(cnv_median, stat(density), color=cancer)) +
  geom_density(alpha=0.1, size = 1) +
  xlab("Median of CNV of samples") +
  # facet_wrap(~cancer, ncol=3) +
  coord_cartesian(xlim=c(-0.25,0.25))

ggplot(summary_pan, aes(cnv_sd, stat(density), color=cancer)) +
  geom_density(alpha=0.1, size = 1) +
  xlab("SD of CNV of samples") +
  # facet_wrap(~cancer, ncol=3) +
  coord_cartesian(xlim=c(0,1))

PCA <- na.omit(cnv_pan)
PCA2 <- t(PCA)
pcaResults <- prcomp(PCA2)
PCA2 <- as.data.frame(PCA2)
PCA2$cancer <- factor(indep_pan$cancer, levels=unique(indep_pan$cancer))
# autoplot(pcaResults,data=PCA2,colour="cancer")+geom_text_repel(label=rownames(PCA2),aes(color=PCA2$aneuploidy))
autoplot(pcaResults,data=PCA2,colour="cancer")

PCA <- na.omit(protein_pan)
PCA2 <- t(PCA)
pcaResults <- prcomp(PCA2)
PCA2 <- as.data.frame(PCA2)
PCA2$cancer <- factor(indep_pan$cancer, levels=unique(indep_pan$cancer))
# autoplot(pcaResults,data=PCA2,colour="cancer")+geom_text_repel(label=rownames(PCA2),aes(color=PCA2$aneuploidy))
autoplot(pcaResults,data=PCA2,colour="cancer")

ggplot(summary_pan, aes(protein_median, stat(density), color=cancer)) +
  geom_density(alpha=0.1, size = 1) +
  xlab("Median of protein of samples") +
  # facet_wrap(~cancer, ncol=3) +
  coord_cartesian(xlim=c(-0.25,1.25))

ggplot(summary_pan, aes(protein_sd, stat(density), color=cancer)) +
  geom_density(alpha=0.1, size = 1) +
  xlab("SD of protein of samples") +
  # facet_wrap(~cancer, ncol=3) +
  coord_cartesian(xlim=c(0,3))



### cell cycle score? purity?





