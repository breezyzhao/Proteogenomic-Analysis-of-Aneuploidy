library(ggplot2)
library(ggfortify)
library(ggrepel)
rm(list=ls())

### colon2 ###
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/colon_omics_log2FC.RData")
cnv_colon2 <- dna_fc
rna_colon2 <- rna_fc
protein_colon2 <- protein_fc

rm(list=c("dna_fc", "rna_fc", "protein_fc"))
### colon2 ###

### breast2 ###
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/breast_omics_log2FC.RData")
cnv_breast2 <- dna_fc
rna_breast2 <- rna_fc
protein_breast2 <- protein_fc

rm(list=c("dna_fc", "rna_fc", "protein_fc"))
### breast2 ###

### ovarian2 ###
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/ovarian_omics_log2FC.RData")
cnv_ovarian2 <- dna_fc
rna_ovarian2 <- rna_fc
protein_ovarian2 <- protein_fc

rm(list=c("dna_fc", "rna_fc", "protein_fc"))
### ovarian2 ###

### ccrcc ###
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/ccrcc_omics_log2FC.RData")
cnv_ccrcc <- dna_fc
rna_ccrcc <- rna_fc
protein_ccrcc <- protein_fc

rm(list=c("dna_fc", "rna_fc", "protein_fc"))
### ccrcc ###

### endometrial ###
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/endometrail_omics_log2FC.RData")
index <- !grepl("C3N.01003", colnames(dna_fc))
cnv_endometrial <- dna_fc[,index]
rna_endometrial <- rna_fc[,index]
protein_endometrial <- protein_fc[,index]


rm(list=c("dna_fc", "rna_fc", "protein_fc"))
### endometrail ###

### hnscc ###
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/hnscc_omics_log2FC.RData")
cnv_hnscc <- dna_fc
rna_hnscc <- rna_fc
protein_hnscc <- protein_fc

rm(list=c("dna_fc", "rna_fc", "protein_fc"))
### hnscc ###

### LUAD ###
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/luad_omics_log2FC.RData")
index <- !grepl("C3L.00510", colnames(dna_fc))
cnv_luad <- dna_fc[,index]
rna_luad <- rna_fc[,index]
protein_luad <- protein_fc[,index]

rm(list=c("dna_fc", "rna_fc", "protein_fc"))
### luad ###

### check the dirstribution of CNV
cnv_plot_colon2 <- data.frame(cnv=as.vector(as.matrix(cnv_colon2)), cancer="colon2")
cnv_plot_breast2 <- data.frame(cnv=as.vector(as.matrix(cnv_breast2)), cancer="breast2")
cnv_plot_ovarian2 <- data.frame(cnv=as.vector(as.matrix(cnv_ovarian2)), cancer="ovarian2")
cnv_plot_ccrcc <- data.frame(cnv=as.vector(as.matrix(cnv_ccrcc)), cancer="ccrcc")
cnv_plot_endometrial <- data.frame(cnv=as.vector(as.matrix(cnv_endometrial)), cancer="endometrial")
cnv_plot_hnscc <- data.frame(cnv=as.vector(as.matrix(cnv_hnscc)), cancer="hnscc")
cnv_plot_luad <- data.frame(cnv=as.vector(as.matrix(cnv_hnscc)), cancer="luad")
cnv_plot <- rbind(cnv_plot_colon2, cnv_plot_breast2, cnv_plot_ovarian2, cnv_plot_ccrcc, cnv_plot_endometrial, cnv_plot_hnscc, cnv_plot_luad)
cnv_plot$cancer <- factor(cnv_plot$cancer, levels=unique(cnv_plot$cancer))

ggplot(cnv_plot, aes(cnv, stat(density))) +
  geom_density(alpha=0.1, size = 1) +
  xlab("CNV") +
  facet_wrap(~cancer, nrow=1) +
  coord_cartesian(xlim=c(-2,2))

### conbine cnv data
cnv_pan <- merge(cnv_colon2, cnv_breast2, by.x="row.names", by.y="row.names")
cnv_pan <- merge(cnv_pan, cnv_ovarian2, by.x="Row.names", by.y="row.names")
cnv_pan <- merge(cnv_pan, cnv_ccrcc, by.x="Row.names", by.y="row.names")
cnv_pan <- merge(cnv_pan, cnv_endometrial, by.x="Row.names", by.y="row.names")
cnv_pan <- merge(cnv_pan, cnv_hnscc, by.x="Row.names", by.y="row.names")
cnv_pan <- merge(cnv_pan, cnv_luad, by.x="Row.names", by.y="row.names")
rownames(cnv_pan) <- cnv_pan$Row.names
cnv_pan <- cnv_pan[,-1]

### check the dirstribution of RNA
rna_plot_colon2 <- data.frame(rna=as.vector(as.matrix(rna_colon2)), cancer="colon2")
rna_plot_breast2 <- data.frame(rna=as.vector(as.matrix(rna_breast2)), cancer="breast2")
rna_plot_ovarian2 <- data.frame(rna=as.vector(as.matrix(rna_ovarian2)), cancer="ovarian2")
rna_plot_ccrcc <- data.frame(rna=as.vector(as.matrix(rna_ccrcc)), cancer="ccrcc")
rna_plot_endometrial <- data.frame(rna=as.vector(as.matrix(rna_endometrial)), cancer="endometrial")
rna_plot_hnscc <- data.frame(rna=as.vector(as.matrix(rna_hnscc)), cancer="hnscc")
rna_plot_luad <- data.frame(rna=as.vector(as.matrix(rna_hnscc)), cancer="luad")
rna_plot <- rbind(rna_plot_colon2, rna_plot_breast2, rna_plot_ovarian2, rna_plot_ccrcc, rna_plot_endometrial, rna_plot_hnscc, rna_plot_luad)
rna_plot$cancer <- factor(rna_plot$cancer, levels=unique(rna_plot$cancer))

ggplot(rna_plot, aes(rna, stat(density))) +
  geom_density(alpha=0.1, size = 1) +
  xlab("RNA") +
  coord_cartesian(xlim=c(-2,2)) +
  facet_wrap(~cancer, nrow=1)

rna_pan <- merge(rna_colon2, rna_breast2, by.x="row.names", by.y="row.names")
rna_pan <- merge(rna_pan, rna_ovarian2, by.x="Row.names", by.y="row.names")
rna_pan <- merge(rna_pan, rna_ccrcc, by.x="Row.names", by.y="row.names")
rna_pan <- merge(rna_pan, rna_endometrial, by.x="Row.names", by.y="row.names")
rna_pan <- merge(rna_pan, rna_hnscc, by.x="Row.names", by.y="row.names")
rna_pan <- merge(rna_pan, rna_luad, by.x="Row.names", by.y="row.names")
rownames(rna_pan) <- rna_pan$Row.names
rna_pan <- rna_pan[,-1]

### check the dirstribution of protein
protein_plot_colon2 <- data.frame(protein=as.vector(as.matrix(protein_colon2)), cancer="colon2")
protein_plot_breast2 <- data.frame(protein=as.vector(as.matrix(protein_breast2)), cancer="breast2")
protein_plot_ovarian2 <- data.frame(protein=as.vector(as.matrix(protein_ovarian2)), cancer="ovarian2")
protein_plot_ccrcc <- data.frame(protein=as.vector(as.matrix(protein_ccrcc)), cancer="ccrcc")
protein_plot_endometrial <- data.frame(protein=as.vector(as.matrix(protein_endometrial)), cancer="endometrial")
protein_plot_hnscc <- data.frame(protein=as.vector(as.matrix(protein_hnscc)), cancer="hnscc")
protein_plot_luad <- data.frame(protein=as.vector(as.matrix(protein_hnscc)), cancer="luad")
protein_plot <- rbind(protein_plot_colon2, protein_plot_breast2, protein_plot_ovarian2, protein_plot_ccrcc, protein_plot_endometrial, protein_plot_hnscc, protein_plot_luad)
protein_plot$cancer <- factor(protein_plot$cancer, levels=unique(protein_plot$cancer))

ggplot(protein_plot, aes(protein, stat(density))) +
  geom_density(alpha=0.1, size = 1) +
  xlab("Protein") +
  coord_cartesian(xlim=c(-2,2)) +
  facet_wrap(~cancer, nrow=1)

### conbine protein data
# protein_pan <- merge(protein_colon2, protein_breast2, by.x="row.names", by.y="row.names", all=T)
# protein_pan <- merge(protein_pan, protein_ovarian2, by.x="Row.names", by.y="row.names", all=T)
# protein_pan <- merge(protein_pan, protein_ccrcc, by.x="Row.names", by.y="row.names", all=T)
# protein_pan <- merge(protein_pan, protein_endometrial, by.x="Row.names", by.y="row.names", all=T)
# protein_pan <- merge(protein_pan, protein_hnscc, by.x="Row.names", by.y="row.names", all=T)
# protein_pan <- merge(protein_pan, protein_luad, by.x="Row.names", by.y="row.names", all=T)
protein_pan <- merge(protein_colon2, protein_breast2, by.x="row.names", by.y="row.names")
protein_pan <- merge(protein_pan, protein_ovarian2, by.x="Row.names", by.y="row.names")
protein_pan <- merge(protein_pan, protein_ccrcc, by.x="Row.names", by.y="row.names")
protein_pan <- merge(protein_pan, protein_endometrial, by.x="Row.names", by.y="row.names")
protein_pan <- merge(protein_pan, protein_hnscc, by.x="Row.names", by.y="row.names")
protein_pan <- merge(protein_pan, protein_luad, by.x="Row.names", by.y="row.names")
rownames(protein_pan) <- protein_pan$Row.names
protein_pan <- protein_pan[,-1]
protein_pan <- data.frame(scale(protein_pan, center=T, scale=F))

### calculate the median and variation of each sample
if ((sum(colnames(cnv_pan)!=colnames(rna_pan)) + sum(colnames(cnv_pan)!=colnames(protein_pan)))!=0) {stop("check the ordering of samples")}

summary_pan <- data.frame(patients=colnames(cnv_pan),
                          cancer=c(rep("colon",95),rep("breast",87),rep("ovarian",81),rep("ccrcc",110),rep("endometrial",94),rep("hnscc",107),rep("luad",108)),
                          cnv_median=apply(cnv_pan,2,median,na.rm=T),
                          cnv_sd=apply(cnv_pan,2,sd,na.rm=T),
                          rna_median=apply(rna_pan,2,median,na.rm=T),
                          rna_sd=apply(rna_pan,2,sd,na.rm=T),
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
  coord_cartesian(xlim=c(0,2))

ggplot(summary_pan, aes(rna_median, stat(density), color=cancer)) +
  geom_density(alpha=0.1, size = 1) +
  xlab("Median of RNA of samples") +
  # facet_wrap(~cancer, ncol=3) +
  coord_cartesian(xlim=c(-0.25,0.25))

ggplot(summary_pan, aes(rna_sd, stat(density), color=cancer)) +
  geom_density(alpha=0.1, size = 1) +
  xlab("SD of RNA of samples") +
  # facet_wrap(~cancer, ncol=3) +
  coord_cartesian(xlim=c(0,2))

ggplot(summary_pan, aes(protein_median, stat(density), color=cancer)) +
  geom_density(alpha=0.1, size = 1) +
  xlab("Median of protein of samples") +
  # facet_wrap(~cancer, ncol=3) +
  coord_cartesian(xlim=c(-0.25,0.25))

ggplot(summary_pan, aes(protein_sd, stat(density), color=cancer)) +
  geom_density(alpha=0.1, size = 1) +
  xlab("SD of protein of samples") +
  # facet_wrap(~cancer, ncol=3) +
  coord_cartesian(xlim=c(0,2))

### PCA to check
PCA <- na.omit(cnv_pan)
PCA2 <- t(PCA)
pcaResults <- prcomp(PCA2)
PCA2 <- as.data.frame(PCA2)
PCA2$cancer <- summary_pan$cancer
# autoplot(pcaResults,data=PCA2,colour="cancer")+geom_text_repel(label=rownames(PCA2),aes(color=PCA2$aneuploidy))
autoplot(pcaResults,data=PCA2,colour="cancer")

PCA <- na.omit(rna_pan)
PCA2 <- t(PCA)
pcaResults <- prcomp(PCA2)
PCA2 <- as.data.frame(PCA2)
PCA2$cancer <- summary_pan$cancer
# autoplot(pcaResults,data=PCA2,colour="cancer")+geom_text_repel(label=rownames(PCA2),aes(color=PCA2$aneuploidy))
autoplot(pcaResults,data=PCA2,colour="cancer")

PCA <- na.omit(protein_pan)
PCA2 <- t(PCA)
pcaResults <- prcomp(PCA2)
PCA2 <- as.data.frame(PCA2)
PCA2$cancer <- summary_pan$cancer
# autoplot(pcaResults,data=PCA2,colour="cancer")+geom_text_repel(label=rownames(PCA2),aes(color=PCA2$aneuploidy))
autoplot(pcaResults,data=PCA2,colour="cancer")

cnv_pan_plot <- data.frame(cancer=c(rep("colon",95),rep("breast",87),rep("ovarian",81),rep("ccrcc",110),rep("endometrial",94),rep("hnscc",107),rep("luad",108)),
                           cnv_log2FC=as.vector(as.matrix(cnv_pan)))
ggplot(cnv_pan_plot, aes(x=cnv_log2FC)) + 
  geom_histogram(binwidth=0.1, color="black", fill="white") + 
  coord_cartesian(xlim=c(-1,1)) +
  geom_vline(aes(xintercept=0.2), color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=-0.2), color="blue", linetype="dashed", size=1) + 
  geom_vline(aes(xintercept=0.65), color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=-0.65), color="red", linetype="dashed", size=1)

save(cnv_pan, rna_pan, protein_pan, file = "/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/pan_omics_log2FC.RData")

