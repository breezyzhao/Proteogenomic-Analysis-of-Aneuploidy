setwd("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC")

library(ggplot2)
library(dplyr)
library(reshape2)

rm(list=ls())

cancer="luad"

if (cancer=="colon2") {load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/colon_omics_log2FC.RData")}
if (cancer=="ccrcc") {load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/ccrcc_omics_log2FC.RData")}
if (cancer=="hnscc") {load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/hnscc_omics_log2FC.RData")}
if (cancer=="luad") {load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/luad_omics_log2FC.RData")}

### calculate the # of samples in different classes
threshold1 <- 0.2
threshold2 <- 0.65

# x <- dna_fc[1,]
count_class <- function(x, threshold1, threshold2) {
  gain2 <- sum(x>threshold2, na.rm=T)
  gain1 <- sum(x<=threshold2 & x>threshold1, na.rm=T)
  unchange <- sum(x<=threshold1 & x>=threshold1*(-1), na.rm=T)
  loss1 <- sum(x<threshold1*(-1) & x>=threshold2*(-1), na.rm=T)
  loss2 <- sum(x<threshold2*(-1), na.rm=T)
  count <- c("loss+"=loss2, "loss"=loss1, "unchange"=unchange, "gain"=gain1, "gain+"=gain2)
}

count <- data.frame(t(apply(dna_fc, 1, function (x) count_class(x, threshold1, threshold2))))
colnames(count) <- c("loss+", "loss", "unchange", "gain", "gain+")

### add CORUM information
corum <- as.data.frame(read.delim("/Volumes/davolt01lab/davolt01labspace/Raquel/database/coreComplexesv3.0.txt", sep = "\t"))
master_list_names <- strsplit(as.character(corum$subunits.Gene.name.), split = ";")
master_list_names <- unique(as.character(unlist(master_list_names)))
count$CORUM <- rownames(count) %in% master_list_names
count$CORUM <- factor(count$CORUM, levels=c(FALSE, TRUE), labels=c("non-CORUM", "CORUM"))

### optimize for colon2
# count <- count[order(count$`gain+`, decreasing=T),]
# gene_gain_NCORUM <- rownames(count)[which(count$CORUM=="non-CORUM")[2]]
# gene_gain_CORUM <- rownames(count)[which(count$CORUM=="CORUM")[2]]
# count <- count[order(count$`loss`, decreasing=T),]
# gene_loss_NCORUM <- rownames(count)[which(count$CORUM=="non-CORUM")[7]]
# gene_loss_CORUM <- rownames(count)[which(count$CORUM=="CORUM")[3]]

### optimize for ccrcc
# count <- count[order(count$`gain`, decreasing=T),]
# gene_gain_NCORUM <- rownames(count)[which(count$CORUM=="non-CORUM")[6]]
# gene_gain_CORUM <- rownames(count)[which(count$CORUM=="CORUM")[2]]
# count <- count[order(count$`loss`, decreasing=T),]
# gene_loss_NCORUM <- rownames(count)[which(count$CORUM=="non-CORUM")[100]]
# gene_loss_CORUM <- rownames(count)[which(count$CORUM=="CORUM")[19]]

### optimize for hnscc
# count <- count[order(count$`gain`, decreasing=T),]
# gene_gain_NCORUM <- rownames(count)[which(count$CORUM=="non-CORUM")[10]]
# gene_gain_CORUM <- rownames(count)[which(count$CORUM=="CORUM")[2]]
# count <- count[order(count$`loss`, decreasing=T),]
# gene_loss_NCORUM <- rownames(count)[which(count$CORUM=="non-CORUM")[2]]
# gene_loss_CORUM <- rownames(count)[which(count$CORUM=="CORUM")[5]]

### optimize for luad
count <- count[order(count$`gain`, decreasing=T),]
gene_gain_NCORUM <- rownames(count)[which(count$CORUM=="non-CORUM")[5]]
gene_gain_CORUM <- rownames(count)[which(count$CORUM=="CORUM")[18]]
count <- count[order(count$`loss`, decreasing=T),]
gene_loss_NCORUM <- rownames(count)[which(count$CORUM=="non-CORUM")[4]]
gene_loss_CORUM <- rownames(count)[which(count$CORUM=="CORUM")[14]]

### plot RNA and Protein log2FC of selected genes
if (sum(colnames(dna_fc)!=colnames(rna_fc))+sum(colnames(dna_fc)!=colnames(protein_fc))!=0) {stop("check gene ordering!")}
data_gene_NCORUM <- data.frame(DNA_gain_gene=as.numeric(dna_fc[gene_gain_NCORUM,]),
                               RNA_gain_gene=as.numeric(rna_fc[gene_gain_NCORUM,]),
                               Protein_gain_gene=as.numeric(protein_fc[gene_gain_NCORUM,]),
                               DNA_loss_gene=as.numeric(dna_fc[gene_loss_NCORUM,]),
                               RNA_loss_gene=as.numeric(rna_fc[gene_loss_NCORUM,]),
                               Protein_loss_gene=as.numeric(protein_fc[gene_loss_NCORUM,]))
rownames(data_gene_NCORUM) <- colnames(dna_fc)

x <- data_gene_NCORUM$DNA_gain_gene
determine_class <- function(x, threshold1, threshold2) {
  class <- rep(NA, length(x))
  class[x>threshold2] <- "gain+"
  class[x<=threshold2 & x>threshold1] <- "gain"
  class[x<=threshold1 & x>=threshold1*(-1)] <- "unchange"
  class[x<threshold1*(-1) & x>=threshold2*(-1)] <- "loss"
  class[x<threshold2*(-1)] <- "loss+"
  return(class)
}
data_gene_NCORUM$class_gain <- unlist(determine_class(data_gene_NCORUM$DNA_gain_gene, threshold1, threshold2))
data_gene_NCORUM$class_loss <- unlist(determine_class(data_gene_NCORUM$DNA_loss_gene, threshold1, threshold2))
data_gene_NCORUM$class_gain <- factor(data_gene_NCORUM$class_gain, levels=c("loss+","loss","unchange","gain","gain+"))
data_gene_NCORUM$class_loss <- factor(data_gene_NCORUM$class_loss, levels=c("loss+","loss","unchange","gain","gain+"))
data_gene_NCORUM1 <- data_gene_NCORUM[,c(7,1:3)]
data_gene_NCORUM2 <- data_gene_NCORUM[,c(8,4:6)]
data_gene_NCORUM1 <- data_gene_NCORUM1[(data_gene_NCORUM1$class_gain!="loss" & data_gene_NCORUM1$class_gain!="loss+"),]
data_gene_NCORUM2 <- data_gene_NCORUM2[(data_gene_NCORUM2$class_loss!="gain" & data_gene_NCORUM2$class_loss!="gain+"),]
data_gene_NCORUM1 <- reshape2::melt(data_gene_NCORUM1, id="class_gain")
data_gene_NCORUM1$class2 <- "gain"
colnames(data_gene_NCORUM1)[1] <- "class"
data_gene_NCORUM2 <- reshape2::melt(data_gene_NCORUM2, id="class_loss")
data_gene_NCORUM2$class2 <- "loss"
colnames(data_gene_NCORUM2)[1] <- "class"
data_gene_NCORUM_reshape <- rbind(data_gene_NCORUM1, data_gene_NCORUM2)
data_gene_NCORUM_reshape$class3 <- paste0(data_gene_NCORUM_reshape$class2, "_", data_gene_NCORUM_reshape$class)
data_gene_NCORUM_reshape$class3 <- factor(data_gene_NCORUM_reshape$class3, levels=c("loss_loss+", "loss_loss", "loss_unchange", "gain_unchange", "gain_gain", "gain_gain+"))

# pdf(paste0("sample_", cancer, "_", gene_gain_NCORUM, "_", gene_loss_NCORUM, "_log2FC.pdf"), width=8, height=4)
# ggplot(data_gene_NCORUM_reshape, aes(x=class3, y=value, fill=variable)) +
#   geom_boxplot(alpha=0.9, outlier.shape = NA) +
#   coord_cartesian(ylim=c(-2,2)) +
#   labs(x="genes classfied based on CNV", y="log2FC") +
#   # scale_x_discrete(labels=c("loss+","loss","unchange","unchange","gain","gain+")) +
#   # scale_x_discrete(labels=c("loss+","loss","unchange","unchange","gain")) +
#   scale_fill_manual(values=c("DNA_loss_gene"="#488f31","RNA_loss_gene"="#8fbc7b","Protein_loss_gene"="#d2e9c6",
#                              "DNA_gain_gene"="#de425b","RNA_gain_gene"="#f49097","Protein_gain_gene"="#ffd5d7"))
# dev.off()

levels(data_gene_NCORUM_reshape$class3) <- c("Deep loss","Loss", "Neutral_1", "Neutral_2", "Gain", "Profound gain")

pdf(paste0("sample_", cancer, "_", gene_gain_NCORUM, "_", gene_loss_NCORUM, "_log2FC.pdf"), width=9, height=3.2)
ggplot(data_gene_NCORUM_reshape, aes(x=class3, y=value, fill=variable)) +
  geom_boxplot(alpha=0.5, outlier.shape = NA, lwd=0.4) +
  coord_cartesian(ylim=c(-2,2)) +
  labs(x="", y="log2FC") +
  theme_minimal() +
  scale_fill_manual(values=c("DNA_loss_gene"="#488f31","RNA_loss_gene"="#8fbc7b","Protein_loss_gene"="#d2e9c6",
                             "DNA_gain_gene"="#de425b","RNA_gain_gene"="#f49097","Protein_gain_gene"="#ffd5d7")) +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45,size = 12,hjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid=element_line(size = 0.4)
  ) 
dev.off()

### repeat for corum gene
data_gene_CORUM <- data.frame(DNA_gain_gene=as.numeric(dna_fc[gene_gain_CORUM,]),
                              RNA_gain_gene=as.numeric(rna_fc[gene_gain_CORUM,]),
                              Protein_gain_gene=as.numeric(protein_fc[gene_gain_CORUM,]),
                              DNA_loss_gene=as.numeric(dna_fc[gene_loss_CORUM,]),
                              RNA_loss_gene=as.numeric(rna_fc[gene_loss_CORUM,]),
                              Protein_loss_gene=as.numeric(protein_fc[gene_loss_CORUM,]))
rownames(data_gene_CORUM) <- colnames(dna_fc)

data_gene_CORUM$class_gain <- unlist(determine_class(data_gene_CORUM$DNA_gain_gene, threshold1, threshold2))
data_gene_CORUM$class_loss <- unlist(determine_class(data_gene_CORUM$DNA_loss_gene, threshold1, threshold2))
data_gene_CORUM$class_gain <- factor(data_gene_CORUM$class_gain, levels=c("loss+","loss","unchange","gain","gain+"))
data_gene_CORUM$class_loss <- factor(data_gene_CORUM$class_loss, levels=c("loss+","loss","unchange","gain","gain+"))
data_gene_CORUM1 <- data_gene_CORUM[,c(7,1:3)]
data_gene_CORUM2 <- data_gene_CORUM[,c(8,4:6)]
data_gene_CORUM1 <- data_gene_CORUM1[(data_gene_CORUM1$class_gain!="loss" & data_gene_CORUM1$class_gain!="loss+"),]
data_gene_CORUM2 <- data_gene_CORUM2[(data_gene_CORUM2$class_loss!="gain" & data_gene_CORUM2$class_loss!="gain+"),]
data_gene_CORUM1 <- reshape2::melt(data_gene_CORUM1, id="class_gain")
data_gene_CORUM1$class2 <- "gain"
colnames(data_gene_CORUM1)[1] <- "class"
data_gene_CORUM2 <- reshape2::melt(data_gene_CORUM2, id="class_loss")
data_gene_CORUM2$class2 <- "loss"
colnames(data_gene_CORUM2)[1] <- "class"
data_gene_CORUM_reshape <- rbind(data_gene_CORUM1, data_gene_CORUM2)
data_gene_CORUM_reshape$class3 <- paste0(data_gene_CORUM_reshape$class2, "_", data_gene_CORUM_reshape$class)
data_gene_CORUM_reshape$class3 <- factor(data_gene_CORUM_reshape$class3, levels=c("loss_loss+", "loss_loss", "loss_unchange", "gain_unchange", "gain_gain", "gain_gain+"))


# pdf(paste0("sample_", cancer, "_", gene_gain_CORUM, "_", gene_loss_CORUM, "_CORUM_log2FC.pdf"), width=8, height=4)
# ggplot(data_gene_CORUM_reshape, aes(x=class3, y=value, fill=variable)) +
#   geom_boxplot(alpha=0.9, outlier.shape = NA) +
#   coord_cartesian(ylim=c(-2,2)) +
#   labs(x="genes classfied based on CNV", y="log2FC") +
#   scale_x_discrete(labels=c("loss+","loss","unchange","unchange","gain","gain+")) +
#   scale_fill_manual(values=c("DNA_loss_gene"="#488f31","RNA_loss_gene"="#8fbc7b","Protein_loss_gene"="#d2e9c6",
#                              "DNA_gain_gene"="#de425b","RNA_gain_gene"="#f49097","Protein_gain_gene"="#ffd5d7"))
# dev.off()

levels(data_gene_CORUM_reshape$class3) <- c("Deep loss","Loss", "Neutral_1", "Neutral_2", "Gain", "Profound gain")

pdf(paste0("sample_", cancer, "_", gene_gain_CORUM, "_", gene_loss_CORUM, "_CORUM_log2FC.pdf"), width=9, height=3.2)
ggplot(data_gene_CORUM_reshape, aes(x=class3, y=value, fill=variable)) +
  geom_boxplot(alpha=0.5, outlier.shape = NA, lwd=0.4) +
  coord_cartesian(ylim=c(-2,2)) +
  labs(x="", y="log2FC") +
  theme_minimal() +
  scale_fill_manual(values=c("DNA_loss_gene"="#488f31","RNA_loss_gene"="#8fbc7b","Protein_loss_gene"="#d2e9c6",
                             "DNA_gain_gene"="#de425b","RNA_gain_gene"="#f49097","Protein_gain_gene"="#ffd5d7")) +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45,size = 12,hjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid=element_line(size = 0.4)
  ) 
dev.off()