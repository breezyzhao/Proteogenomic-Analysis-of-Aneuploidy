### this script is to quantify the protein compensation in each sample
### test from colon2 data
library(biomaRt)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(boot)
library(boot.pval)
library(simpleboot)

setwd("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC")

rm(list=ls())

cancer="luad"

### read colon2 data
# colon2_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/colon2_addMutation.RData"
# load(colon2_path)

### read breast2 data
# breast2_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/breast2_addMutation.RData"
# load(breast2_path)
# cnv0 <- as.data.frame(apply(cnv1,2,as.numeric))
# rownames(cnv0) <- rownames(cnv1)
# cnv1 <- cnv0
# protein0 <- as.data.frame(apply(protein1,2,as.numeric))
# rownames(protein0) <- rownames(protein1)
# protein1 <- protein0

### read ovarian2 data
# ovarian2_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data/ovarian2.RData"
# load(ovarian2_path)
# rna_all_filter <- log2(rna1+1)
# indepVar <- AS1

### read ccrcc data
# ccrcc_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/ccrcc_addMutation.RData"
# load(ccrcc_path)
# cnv0 <- as.data.frame(apply(cnv1,2,as.numeric))
# rownames(cnv0) <- rownames(cnv1)
# cnv1 <- cnv0
# rna_all_filter <- log2(rna_all_filter+1)
# protein0 <- as.data.frame(apply(protein1,2,as.numeric))
# rownames(protein0) <- rownames(protein1)
# protein1 <- protein0

### read endometrail data
# endometrial_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/endometrail_addMutation.RData"
# load(endometrial_path)

### read hnscc data
# hnscc_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/hnscc_addMutation.RData"
# load(hnscc_path)

### read luad data
luad_path <- "/Volumes/davolt01lab/davolt01labspace/Pan/Aneuploidy/Tumor Dataset/GLM-GSEA/Organized Data_addMutation/luad_addMutation.RData"
load(luad_path)
protein0 <- as.data.frame(apply(protein1,2,as.numeric))
rownames(protein0) <- rownames(protein1)
protein1 <- protein0

### start analysis
id_patient <- intersect(colnames(cnv1), colnames(rna_all_filter))
id_patient <- intersect(id_patient, colnames(protein1))
id_patient <- intersect(id_patient, rownames(indepVar))

### for endometrial and luad, remove weird samples
if (cancer=="endometrial") {id_patient <- id_patient[!grepl("C3N.01003", id_patient)]}
if (cancer=="luad") {id_patient <- id_patient[!grepl("C3L.00510", id_patient)]}

### remove genes of low expression at RNA (bottom 10%)
rna_all_filter_2 <- rna_all_filter[,id_patient]
geneLevel <- data.frame(gene=rownames(rna_all_filter_2),
                        median_RNA=matrixStats::rowMedians(as.matrix(rna_all_filter_2), na.rm=T))
geneLevel <- geneLevel[order(geneLevel$median_RNA, decreasing=T),]
id_gene <- geneLevel$gene[1:round(0.9*nrow(geneLevel))]
id_gene <- intersect(id_gene, rownames(cnv1))
id_gene <- intersect(id_gene, rownames(protein1))

### there are 95 patients and 7386 genes for colon2
### there are 87 patients and 8711 genes for breast2
### there are 81 patients and 7638 genes for breast2
cnv_filter <- cnv1[id_gene, id_patient]
rna_filter <- rna_all_filter[id_gene, id_patient]
protein_filter <- protein1[id_gene, id_patient]
indepVar_filter <- indepVar[id_patient,]

### calculate the basal expression at genewise
threshold_gene <- 0.20
basal <- data.frame(gene=rownames(cnv_filter),
                    cnv_mean=rep(NA, nrow(cnv_filter)),
                    cnv_median=rep(NA, nrow(cnv_filter)),
                    cnv_n=rep(NA, nrow(cnv_filter)),
                    rna_mean=rep(NA, nrow(cnv_filter)),
                    rna_median=rep(NA, nrow(cnv_filter)),
                    rna_n=rep(NA, nrow(cnv_filter)),
                    protein_mean=rep(NA, nrow(cnv_filter)),
                    protein_median=rep(NA, nrow(cnv_filter)),
                    protein_n=rep(NA, nrow(cnv_filter)))

if (sum(colnames(cnv_filter)!=colnames(rna_filter) | colnames(cnv_filter)!=colnames(protein_filter))!=0) {stop(paste0("gene are not corresponding!"))}
for (i in 1:nrow(cnv_filter)) {
  if (rownames(cnv_filter)[i]!=rownames(rna_filter)[i] | rownames(cnv_filter)[i]!=rownames(protein_filter)[i]) {stop(paste0("gene (#",i,") are not corresponding!"))}
  cnv_gene <- cnv_filter[i,]
  rna_gene <- rna_filter[i,]
  protein_gene <- protein_filter[i,]
  index <- (abs(cnv_gene)<=threshold_gene) & (!is.na(cnv_gene))
  if (sum(index)<=3) {next} else {
    basal$cnv_mean[i] <- mean(cnv_gene[index])
    basal$cnv_median[i] <- median(cnv_gene[index])
    basal$cnv_n[i] <- sum(index)
    basal$rna_mean[i] <- mean(rna_gene[index], na.rm=T)
    basal$rna_median[i] <- median(rna_gene[index], na.rm=T)
    basal$rna_n[i] <- sum(!is.na(rna_gene[index]))
    basal$protein_mean[i] <- mean(protein_gene[index], na.rm=T)
    basal$protein_median[i] <- median(protein_gene[index], na.rm=T)
    basal$protein_n[i] <- sum(!is.na(protein_gene[index]))
  }
}

### remove genes whose basal protein consist of very few data points (n<10)
basal <- basal[basal$protein_n>=10 & basal$rna_n>=10 & !is.na(basal$protein_n) & !is.na(basal$rna_n),]
id_gene2 <- basal$gene
cnv_filter <- cnv1[id_gene2, id_patient]
rna_filter <- rna_all_filter[id_gene2, id_patient]
protein_filter <- protein1[id_gene2, id_patient]

if (sum(rownames(cnv_filter)!=basal$gene | rownames(rna_filter)!=basal$gene | rownames(protein_filter)!=basal$gene)!=0) {stop("genes are not corresponding")}
dna_fc <- cnv_filter-basal$cnv_median
rna_fc <- rna_filter-basal$rna_median
protein_fc <- protein_filter-basal$protein_median

# save(dna_fc, rna_fc, protein_fc, file = paste0("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/", cancer, "_omics_log2FC.RData"))

### let see how rna and protein change along with cnv at gene level
fc_gene <- data.frame(gene=rep(rownames(dna_fc),ncol(dna_fc)),
                      DNA_log2FC=as.vector(as.matrix(dna_fc)),
                      RNA_log2FC=as.vector(as.matrix(rna_fc)),
                      Protein_log2FC=as.vector(as.matrix(protein_fc)))
fc_gene <- na.omit(fc_gene)
fc_gene$change <- rep("unchange", nrow(fc_gene))
fc_gene$change[fc_gene$DNA_log2FC>threshold_gene] <- "gain"
fc_gene$change[fc_gene$DNA_log2FC<threshold_gene*-1] <- "loss"

### add protein complex information here
corum <- as.data.frame(read.delim("/Volumes/davolt01lab/davolt01labspace/Raquel/database/coreComplexesv3.0.txt", sep = "\t"))
master_list_names <- strsplit(as.character(corum$subunits.Gene.name.), split = ";")
master_list_names <- unique(as.character(unlist(master_list_names)))
fc_gene$CORUM <- fc_gene$gene %in% master_list_names
fc_gene$CORUM <- factor(fc_gene$CORUM, levels=c(FALSE, TRUE), labels=c("non-CORUM", "CORUM"))


### split gain/loss to more groups
threshold_gene2 <- 0.65
fc_gene$change2 <- rep("unchange", nrow(fc_gene))
fc_gene$change2[fc_gene$DNA_log2FC>threshold_gene & fc_gene$DNA_log2FC<threshold_gene2] <- "gain"
fc_gene$change2[fc_gene$DNA_log2FC<threshold_gene*-1 & fc_gene$DNA_log2FC>threshold_gene2*-1] <- "loss"
fc_gene$change2[fc_gene$DNA_log2FC>=threshold_gene2] <- "gain+"
fc_gene$change2[fc_gene$DNA_log2FC<=threshold_gene2*-1] <- "loss+"
fc_gene_plot2 <- reshape2::melt(fc_gene, id=c("gene","CORUM","change","change2"))
fc_gene_plot2$change2 <- factor(fc_gene_plot2$change2, levels=c("loss+","loss", "unchange", "gain","gain+"))

levels(fc_gene_plot2$variable) <- c("DNA", "RNA", "Protein")
levels(fc_gene_plot2$CORUM) <- c("NoCorum","Corum")
levels(fc_gene_plot2$change2) <- c("Deep loss","Loss", "Neutral", "Gain", "Profound gain")

# pdf(paste0("/Users/pc2644/Desktop/CPTAC_",cancer,"_log2FC.pdf"), width=14, height=5)
ggplot(fc_gene_plot2, aes(x=change2, y=value, fill=CORUM)) +
  geom_boxplot(alpha=0.5, outlier.shape = NA, lwd=0.15) +
  coord_cartesian(ylim=c(-2,2)) +
  facet_wrap(~variable, ncol=3) +
  labs(x="", y="log2FC") +
  theme_minimal() +
  scale_fill_manual(values=c("NoCorum"="#AAA900", "Corum"="#6500AA")) +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45,size = 12,hjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid=element_line(size = 0.15)
  ) 
# dev.off()

### calculate the compensation score (group-level)
fc_gene$compensation_RNA <- rep(NA,nrow(fc_gene))
fc_gene$compensation_Protein <- rep(NA,nrow(fc_gene))
fc_gene$compensation_RNA[fc_gene$change=="gain"] <- fc_gene$RNA_log2FC[fc_gene$change=="gain"]-fc_gene$DNA_log2FC[fc_gene$change=="gain"]
fc_gene$compensation_Protein[fc_gene$change=="gain"] <- fc_gene$Protein_log2FC[fc_gene$change=="gain"]-fc_gene$DNA_log2FC[fc_gene$change=="gain"]
fc_gene$compensation_RNA[fc_gene$change=="loss"] <- (fc_gene$RNA_log2FC[fc_gene$change=="loss"]-fc_gene$DNA_log2FC[fc_gene$change=="loss"])*(-1)
fc_gene$compensation_Protein[fc_gene$change=="loss"] <- (fc_gene$Protein_log2FC[fc_gene$change=="loss"]-fc_gene$DNA_log2FC[fc_gene$change=="loss"])*-1

fc_gene$compensation_RNA[fc_gene$change=="unchange" ] <- (fc_gene$RNA_log2FC[fc_gene$change=="unchange"]-fc_gene$DNA_log2FC[fc_gene$change=="unchange"])
fc_gene$compensation_Protein[fc_gene$change=="unchange"] <- (fc_gene$Protein_log2FC[fc_gene$change=="unchange"]-fc_gene$DNA_log2FC[fc_gene$change=="unchange"])

summary_fc_gene <- fc_gene %>%
  group_by(change2, CORUM) %>%
  summarise(n=n(),
            DNA_log2FC=median(DNA_log2FC, na.rm=T),
            RNA_log2FC=median(RNA_log2FC, na.rm=T),
            Protein_log2FC=median(Protein_log2FC, na.rm=T),
            RNA_compensation=median(compensation_RNA, na.rm=T),
            Protein_compensation=median(compensation_Protein, na.rm=T))

# save(summary_fc_gene, file = paste0("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/", cancer, "_summary.RData"))

### generate color for compensation score
compensation <- summary_fc_gene[1:8,-4:-6]
compensation[,-1:-3] <- -1*compensation[,-1:-3]
compensation <- reshape2::melt(compensation, id=c("change2","CORUM","n"))
compensation <- compensation[order(compensation$variable, compensation$change2),]
compensation$location <- 1:nrow(compensation)

# min(compensation$value)
# max(compensation$value)

# pdf(paste0("/Users/pc2644/Desktop/",cancer,"_CPTAC_log2FC_CompensationScore.pdf"), width=14, height=5)
ggplot() + 
  geom_rect(data = compensation, aes(xmin = location-0.2, xmax = location+0.2, ymin = 1.8, ymax = 2.2, fill=value)) + 
  coord_fixed(ratio = 1) +
  theme_void() +
  scale_fill_gradientn(limits=c(-0.4,0.9),
                       breaks=c(-0.4, 0, 0.1, 0.9),
                       colours=c("#c8c8c8", "#f6f6f6", "#b1eeec", "#00aaa9"))
# dev.off()

### boot first version
# Rsti=10000
# boot_median <- function(data, indices){
#   data2 <- data[indices]
#   return(median(data2,na.rm=T))
# }
# 
# set.seed(626)
# 
# boot_RNA_gain <- boot(fc_gene$compensation_RNA[fc_gene$change2=="gain"], boot_median, R=Rsti)
# boot_RNA_gain2 <- boot(fc_gene$compensation_RNA[fc_gene$change2=="gain+"], boot_median, R=Rsti)
# boot_RNA_loss <- boot(fc_gene$compensation_RNA[fc_gene$change2=="loss"], boot_median, R=Rsti)
# boot_RNA_loss2 <- boot(fc_gene$compensation_RNA[fc_gene$change2=="loss+"], boot_median, R=Rsti)
# 
# boot_Protein_gain <- boot(fc_gene$compensation_Protein[fc_gene$change2=="gain"], boot_median, R=Rsti)
# boot_Protein_gain2 <- boot(fc_gene$compensation_Protein[fc_gene$change2=="gain+"], boot_median, R=Rsti)
# boot_Protein_loss <- boot(fc_gene$compensation_Protein[fc_gene$change2=="loss"], boot_median, R=Rsti)
# boot_Protein_loss2 <- boot(fc_gene$compensation_Protein[fc_gene$change2=="loss+"], boot_median, R=Rsti)
# 
# boot <- data.frame(levels=c(rep("RNA",4),rep("Protein",4)),
#                    conditions=rep(c("gain", "high gain", "loss", "deep loss"),2),
#                    CI=rep(NA,8),
#                    pValue=rep(NA,8),
#                    FDR=rep(NA,8))
# 
# boot[1,3] <- paste0("( ",paste0(-1*round(boot.ci(boot_RNA_gain, type="basic")$basic[5:4],4), collapse=", "), " )")
# boot[1,4] <- boot.pval(boot_RNA_gain, theta_null = 0)
# 
# boot[2,3] <- paste0("( ",paste0(-1*round(boot.ci(boot_RNA_gain2, type="basic")$basic[5:4],4), collapse=", "), " )")
# boot[2,4] <- boot.pval(boot_RNA_gain2, theta_null = 0)
# 
# boot[3,3] <- paste0("( ",paste0(-1*round(boot.ci(boot_RNA_loss, type="basic")$basic[5:4],4), collapse=", "), " )")
# boot[3,4] <- boot.pval(boot_RNA_loss, theta_null = 0)
# 
# boot[4,3] <- paste0("( ",paste0(-1*round(boot.ci(boot_RNA_loss2, type="basic")$basic[5:4],4), collapse=", "), " )")
# boot[4,4] <- boot.pval(boot_RNA_loss2, theta_null = 0)
# 
# boot[5,3] <- paste0("( ",paste0(-1*round(boot.ci(boot_Protein_gain, type="basic")$basic[5:4],4), collapse=", "), " )")
# boot[5,4] <- boot.pval(boot_Protein_gain, theta_null = 0)
# 
# boot[6,3] <- paste0("( ",paste0(-1*round(boot.ci(boot_Protein_gain2, type="basic")$basic[5:4],4), collapse=", "), " )")
# boot[6,4] <- boot.pval(boot_Protein_gain2, theta_null = 0)
# 
# boot[7,3] <- paste0("( ",paste0(-1*round(boot.ci(boot_Protein_loss, type="basic")$basic[5:4],4), collapse=", "), " )")
# boot[7,4] <- boot.pval(boot_Protein_loss, theta_null = 0)
# 
# boot[8,3] <- paste0("( ",paste0(-1*round(boot.ci(boot_Protein_loss2, type="basic")$basic[5:4],4), collapse=", "), " )")
# boot[8,4] <- boot.pval(boot_Protein_loss2, theta_null = 0)
# 
# boot$FDR <-  p.adjust(boot$pValue, "BH")
# 
# write.table(boot, file="/Users/pc2644/Desktop/CPTAC_ccrcc_bootstrape.txt", quote=F, sep="\t", row.names=F)

### boot, second version
### bootstrap to calculate CI
# data <- fc_gene$compensation_RNA
# subgroup <- fc_gene$change2=="gain+"
# Rsti=10000
# boot_median <- function(data, indices){
#   data2 <- data[indices]
#   return(median(data2,na.rm=T))
# }
# 
# set.seed(626)
# boot_RNA_gain_CORUM <- boot(fc_gene$compensation_RNA[fc_gene$change2=="gain" & fc_gene$CORUM=="CORUM"], boot_median, R=Rsti)
# boot_RNA_gain2_CORUM <- boot(fc_gene$compensation_RNA[fc_gene$change2=="gain+" & fc_gene$CORUM=="CORUM"], boot_median, R=Rsti)
# boot_RNA_gain_NoCORUM <- boot(fc_gene$compensation_RNA[fc_gene$change2=="gain" & fc_gene$CORUM=="non-CORUM"], boot_median, R=Rsti)
# boot_RNA_gain2_NoCORUM <- boot(fc_gene$compensation_RNA[fc_gene$change2=="gain+" & fc_gene$CORUM=="non-CORUM"], boot_median, R=Rsti)
# boot_RNA_loss_CORUM <- boot(fc_gene$compensation_RNA[fc_gene$change2=="loss" & fc_gene$CORUM=="CORUM"], boot_median, R=Rsti)
# boot_RNA_loss2_CORUM <- boot(fc_gene$compensation_RNA[fc_gene$change2=="loss+" & fc_gene$CORUM=="CORUM"], boot_median, R=Rsti)
# boot_RNA_loss_NoCORUM <- boot(fc_gene$compensation_RNA[fc_gene$change2=="loss" & fc_gene$CORUM=="non-CORUM"], boot_median, R=Rsti, simple=T)
# boot_RNA_loss2_NoCORUM <- boot(fc_gene$compensation_RNA[fc_gene$change2=="loss+" & fc_gene$CORUM=="non-CORUM"], boot_median, R=Rsti)
# 
# boot_Protein_gain_CORUM <- boot(fc_gene$compensation_Protein[fc_gene$change2=="gain" & fc_gene$CORUM=="CORUM"], boot_median, R=Rsti)
# boot_Protein_gain2_CORUM <- boot(fc_gene$compensation_Protein[fc_gene$change2=="gain+" & fc_gene$CORUM=="CORUM"], boot_median, R=Rsti)
# boot_Protein_gain_NoCORUM <- boot(fc_gene$compensation_Protein[fc_gene$change2=="gain" & fc_gene$CORUM=="non-CORUM"], boot_median, R=Rsti)
# boot_Protein_gain2_NoCORUM <- boot(fc_gene$compensation_Protein[fc_gene$change2=="gain+" & fc_gene$CORUM=="non-CORUM"], boot_median, R=Rsti)
# boot_Protein_loss_CORUM <- boot(fc_gene$compensation_Protein[fc_gene$change2=="loss" & fc_gene$CORUM=="CORUM"], boot_median, R=Rsti)
# boot_Protein_loss2_CORUM <- boot(fc_gene$compensation_Protein[fc_gene$change2=="loss+" & fc_gene$CORUM=="CORUM"], boot_median, R=Rsti)
# boot_Protein_loss_NoCORUM <- boot(fc_gene$compensation_Protein[fc_gene$change2=="loss" & fc_gene$CORUM=="non-CORUM"], boot_median, R=Rsti, simple=T)
# boot_Protein_loss2_NoCORUM <- boot(fc_gene$compensation_Protein[fc_gene$change2=="loss+" & fc_gene$CORUM=="non-CORUM"], boot_median, R=Rsti)
# 
# save.image(paste0("boot_", cancer ,"_CPTAC_log2FC.RData"))
# 
# ## calculate the p value
# boot.ci(CI_RNA_gain2, type = "bca")
# boot.ci(boot_Protein_loss_NoCORUM, type = "norm")
# boot.pval(CI_RNA_gain, theta_null = 0.0533)
# plot(CI_RNA_gain2)
# 
# load("boot_pan_CPTAC_log2FC.RData")
# 
# boot <- data.frame(levels=c(rep("RNA",8),rep("Protein",8)),
#                    conditions=rep(c(rep("gain",2), rep("high gain",2), rep("loss",2), rep("deep loss",2)),2),
#                    CORUM=rep(c("CORUM","NoCORUM"),8),
#                    CI=rep(NA,16),
#                    pValue=rep(NA,16))
# 
# printCI <- function(x) {
#   min <- min(-1*round(boot.ci(x, type="basic")$basic[4:5],4))
#   max <- max(-1*round(boot.ci(x, type="basic")$basic[4:5],4))
#   result <- c(min, max)
#   return(result)
# }
# 
# oneTailp <- function(x) {
#   p <- (sum(((x$t)*(-1)-(x$t0)*(-1))>(x$t0)*(-1))+1)/(x$R+1)
#   return(p)
# }
# 
# boot[1,4] <- paste0("( ",paste0(printCI(boot_RNA_gain_CORUM), collapse=", "), " )")
# boot[1,5] <- oneTailp(boot_RNA_gain_CORUM)
# boot[2,4] <- paste0("( ",paste0(printCI(boot_RNA_gain_NoCORUM), collapse=", "), " )")
# boot[2,5] <- oneTailp(boot_RNA_gain_NoCORUM)
# boot[3,4] <- paste0("( ",paste0(printCI(boot_RNA_gain2_CORUM), collapse=", "), " )")
# boot[3,5] <- oneTailp(boot_RNA_gain2_CORUM)
# boot[4,4] <- paste0("( ",paste0(printCI(boot_RNA_gain2_NoCORUM), collapse=", "), " )")
# boot[4,5] <- oneTailp(boot_RNA_gain2_NoCORUM)
# boot[5,4] <- paste0("( ",paste0(printCI(boot_RNA_loss_CORUM), collapse=", "), " )")
# boot[5,5] <- oneTailp(boot_RNA_loss_CORUM)
# boot[6,4] <- paste0("( ",paste0(printCI(boot_RNA_loss_NoCORUM), collapse=", "), " )")
# boot[6,5] <- oneTailp(boot_RNA_loss_NoCORUM)
# boot[7,4] <- paste0("( ",paste0(printCI(boot_RNA_loss2_CORUM), collapse=", "), " )")
# boot[7,5] <- oneTailp(boot_RNA_loss2_CORUM)
# boot[8,4] <- paste0("( ",paste0(printCI(boot_RNA_loss2_NoCORUM), collapse=", "), " )")
# boot[8,5] <- oneTailp(boot_RNA_loss2_NoCORUM)
# 
# boot[9,4] <- paste0("( ",paste0(printCI(boot_Protein_gain_CORUM), collapse=", "), " )")
# boot[9,5] <- oneTailp(boot_Protein_gain_CORUM)
# boot[10,4] <- paste0("( ",paste0(printCI(boot_Protein_gain_NoCORUM), collapse=", "), " )")
# boot[10,5] <- oneTailp(boot_Protein_gain_NoCORUM)
# boot[11,4] <- paste0("( ",paste0(printCI(boot_Protein_gain2_CORUM), collapse=", "), " )")
# boot[11,5] <- oneTailp(boot_Protein_gain2_CORUM)
# boot[12,4] <- paste0("( ",paste0(printCI(boot_Protein_gain2_NoCORUM), collapse=", "), " )")
# boot[12,5] <- oneTailp(boot_Protein_gain2_NoCORUM)
# boot[13,4] <- paste0("( ",paste0(printCI(boot_Protein_loss_CORUM), collapse=", "), " )")
# boot[13,5] <- oneTailp(boot_Protein_loss_CORUM)
# boot[14,4] <- paste0("( ",paste0(printCI(boot_Protein_loss_NoCORUM), collapse=", "), " )")
# boot[14,5] <- oneTailp(boot_Protein_loss_NoCORUM)
# boot[15,4] <- paste0("( ",paste0(printCI(boot_Protein_loss2_CORUM), collapse=", "), " )")
# boot[15,5] <- oneTailp(boot_Protein_loss2_CORUM)
# boot[16,4] <- paste0("( ",paste0(printCI(boot_Protein_loss2_NoCORUM), collapse=", "), " )")
# boot[16,5] <- oneTailp(boot_Protein_loss2_NoCORUM)
# 
# boot$FDR <- p.adjust(boot$pValue, method="BH")
# 
# write.table(boot, file=paste0("/Users/pc2644/Desktop/CPTAC_", cancer, "_bootstrape_single_tail.txt"), quote=F, sep="\t", row.names=F)

### bootstrap to compare median
A1 <- fc_gene$compensation_RNA[fc_gene$change2=="gain" & fc_gene$CORUM=="CORUM"]
B1 <- fc_gene$compensation_RNA[fc_gene$change2=="gain" & fc_gene$CORUM=="non-CORUM"]
n=length(A1)
m=length(B1)
y=c(A1,B1)
camp=data.frame(group=rep(c(1,2),c(n,m)),y)
dif.median=function(data,i) {
  d=data[i,]
  n1=n+1
  m1=n+m
  median(d$y[1:n])-median(d$y[n1:m1])  
}
# RNA_gain_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group, simple=T)
RNA_gain_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group)

A2 <- fc_gene$compensation_RNA[fc_gene$change2=="gain+" & fc_gene$CORUM=="CORUM"]
B2 <- fc_gene$compensation_RNA[fc_gene$change2=="gain+" & fc_gene$CORUM=="non-CORUM"]
n=length(A2)
m=length(B2)
y=c(A2,B2)
camp=data.frame(group=rep(c(1,2),c(n,m)),y)
dif.median=function(data,i) {
  d=data[i,]
  n1=n+1
  m1=n+m
  median(d$y[1:n])-median(d$y[n1:m1])  
}
# RNA_gain2_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group, simple=T)
RNA_gain2_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group)

A3 <- fc_gene$compensation_RNA[fc_gene$change2=="loss" & fc_gene$CORUM=="CORUM"]
B3 <- fc_gene$compensation_RNA[fc_gene$change2=="loss" & fc_gene$CORUM=="non-CORUM"]
n=length(A3)
m=length(B3)
y=c(A3,B3)
camp=data.frame(group=rep(c(1,2),c(n,m)),y)
dif.median=function(data,i) {
  d=data[i,]
  n1=n+1
  m1=n+m
  median(d$y[1:n])-median(d$y[n1:m1])  
}
# RNA_loss_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group, simple=T)
RNA_loss_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group)

A4 <- fc_gene$compensation_RNA[fc_gene$change2=="loss+" & fc_gene$CORUM=="CORUM"]
B4 <- fc_gene$compensation_RNA[fc_gene$change2=="loss+" & fc_gene$CORUM=="non-CORUM"]
n=length(A4)
m=length(B4)
y=c(A4,B4)
camp=data.frame(group=rep(c(1,2),c(n,m)),y)
dif.median=function(data,i) {
  d=data[i,]
  n1=n+1
  m1=n+m
  median(d$y[1:n])-median(d$y[n1:m1])  
}
RNA_loss2_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group)

A5 <- fc_gene$compensation_RNA[fc_gene$change2=="unchange" & fc_gene$CORUM=="CORUM"]
B5 <- fc_gene$compensation_RNA[fc_gene$change2=="unchange" & fc_gene$CORUM=="non-CORUM"]
n=length(A5)
m=length(B5)
y=c(A5,B5)
camp=data.frame(group=rep(c(1,2),c(n,m)),y)
dif.median=function(data,i) {
  d=data[i,]
  n1=n+1
  m1=n+m
  median(d$y[1:n])-median(d$y[n1:m1])  
}
RNA_unchange_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group, simple=T)

A6 <- fc_gene$compensation_Protein[fc_gene$change2=="gain" & fc_gene$CORUM=="CORUM"]
B6 <- fc_gene$compensation_Protein[fc_gene$change2=="gain" & fc_gene$CORUM=="non-CORUM"]
n=length(A6)
m=length(B6)
y=c(A6,B6)
camp=data.frame(group=rep(c(1,2),c(n,m)),y)
dif.median=function(data,i) {
  d=data[i,]
  n1=n+1
  m1=n+m
  median(d$y[1:n])-median(d$y[n1:m1])  
}
Protein_gain_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group)

A7 <- fc_gene$compensation_Protein[fc_gene$change2=="gain+" & fc_gene$CORUM=="CORUM"]
B7 <- fc_gene$compensation_Protein[fc_gene$change2=="gain+" & fc_gene$CORUM=="non-CORUM"]
n=length(A7)
m=length(B7)
y=c(A7,B7)
camp=data.frame(group=rep(c(1,2),c(n,m)),y)
dif.median=function(data,i) {
  d=data[i,]
  n1=n+1
  m1=n+m
  median(d$y[1:n])-median(d$y[n1:m1])  
}
Protein_gain2_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group)

A8 <- fc_gene$compensation_Protein[fc_gene$change2=="loss" & fc_gene$CORUM=="CORUM"]
B8 <- fc_gene$compensation_Protein[fc_gene$change2=="loss" & fc_gene$CORUM=="non-CORUM"]
n=length(A8)
m=length(B8)
y=c(A8,B8)
camp=data.frame(group=rep(c(1,2),c(n,m)),y)
dif.median=function(data,i) {
  d=data[i,]
  n1=n+1
  m1=n+m
  median(d$y[1:n])-median(d$y[n1:m1])  
}
Protein_loss_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group)

A9 <- fc_gene$compensation_Protein[fc_gene$change2=="loss+" & fc_gene$CORUM=="CORUM"]
B9 <- fc_gene$compensation_Protein[fc_gene$change2=="loss+" & fc_gene$CORUM=="non-CORUM"]
n=length(A9)
m=length(B9)
y=c(A9,B9)
camp=data.frame(group=rep(c(1,2),c(n,m)),y)
dif.median=function(data,i) {
  d=data[i,]
  n1=n+1
  m1=n+m
  median(d$y[1:n])-median(d$y[n1:m1])  
}
Protein_loss2_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group)

A10 <- fc_gene$compensation_Protein[fc_gene$change2=="unchange" & fc_gene$CORUM=="CORUM"]
B10 <- fc_gene$compensation_Protein[fc_gene$change2=="unchange" & fc_gene$CORUM=="non-CORUM"]
n=length(A10)
m=length(B10)
y=c(A10,B10)
camp=data.frame(group=rep(c(1,2),c(n,m)),y)
dif.median=function(data,i) {
  d=data[i,]
  n1=n+1
  m1=n+m
  median(d$y[1:n])-median(d$y[n1:m1])  
}
Protein_unchange_CORUMvsNoCORUM <- boot(camp,dif.median,R=10000, strata=camp$group, simple=T)

save.image(paste0("boot_",cancer,"_CPTAC_log2FC_medianDifference.RData"))

boot <- data.frame(levels=c(rep("RNA",5),rep("Protein",5)),
                   conditions=rep(c("gain","gain+","loss","loss+","unchange"),2),
                   CORUM_NoCORUM=rep(NA,10),
                   CI=rep(NA,10),
                   pValue=rep(NA,10))

printCI <- function(x) {
  min <- min(-1*round(boot.ci(x, type="basic")$basic[4:5],4))
  max <- max(-1*round(boot.ci(x, type="basic")$basic[4:5],4))
  result <- c(min, max)
  return(result)
}

twoTailp <- function(x) {
  p <- (sum(abs(x$t-x$t0)>abs(x$t0))+1)/(x$R+1)
  return(p)
}

boot[1,3] <- RNA_gain_CORUMvsNoCORUM$t0*(-1)
boot[1,4] <- paste0("( ",paste0(printCI(RNA_gain_CORUMvsNoCORUM), collapse=", "), " )")
boot[1,5] <- twoTailp(RNA_gain_CORUMvsNoCORUM)

boot[2,3] <- RNA_gain2_CORUMvsNoCORUM$t0*(-1)
boot[2,4] <- paste0("( ",paste0(printCI(RNA_gain2_CORUMvsNoCORUM), collapse=", "), " )")
boot[2,5] <- twoTailp(RNA_gain2_CORUMvsNoCORUM)

boot[3,3] <- RNA_loss_CORUMvsNoCORUM$t0*(-1)
boot[3,4] <- paste0("( ",paste0(printCI(RNA_loss_CORUMvsNoCORUM), collapse=", "), " )")
boot[3,5] <- twoTailp(RNA_loss_CORUMvsNoCORUM)

boot[4,3] <- RNA_loss2_CORUMvsNoCORUM$t0*(-1)
boot[4,4] <- paste0("( ",paste0(printCI(RNA_loss2_CORUMvsNoCORUM), collapse=", "), " )")
boot[4,5] <- twoTailp(RNA_loss2_CORUMvsNoCORUM)

boot[5,3] <- RNA_unchange_CORUMvsNoCORUM$t0*(-1)
boot[5,4] <- paste0("( ",paste0(printCI(RNA_unchange_CORUMvsNoCORUM), collapse=", "), " )")
boot[5,5] <- twoTailp(RNA_unchange_CORUMvsNoCORUM)

boot[6,3] <- Protein_gain_CORUMvsNoCORUM$t0*(-1)
boot[6,4] <- paste0("( ",paste0(printCI(Protein_gain_CORUMvsNoCORUM), collapse=", "), " )")
boot[6,5] <- twoTailp(Protein_gain_CORUMvsNoCORUM)

boot[7,3] <- Protein_gain2_CORUMvsNoCORUM$t0*(-1)
boot[7,4] <- paste0("( ",paste0(printCI(Protein_gain2_CORUMvsNoCORUM), collapse=", "), " )")
boot[7,5] <- twoTailp(Protein_gain2_CORUMvsNoCORUM)

boot[8,3] <- Protein_loss_CORUMvsNoCORUM$t0*(-1)
boot[8,4] <- paste0("( ",paste0(printCI(Protein_loss_CORUMvsNoCORUM), collapse=", "), " )")
boot[8,5] <- twoTailp(Protein_loss_CORUMvsNoCORUM)

boot[9,3] <- Protein_loss2_CORUMvsNoCORUM$t0*(-1)
boot[9,4] <- paste0("( ",paste0(printCI(Protein_loss2_CORUMvsNoCORUM), collapse=", "), " )")
boot[9,5] <- twoTailp(Protein_loss2_CORUMvsNoCORUM)

boot[10,3] <- Protein_unchange_CORUMvsNoCORUM$t0*(-1)
boot[10,4] <- paste0("( ",paste0(printCI(Protein_unchange_CORUMvsNoCORUM), collapse=", "), " )")
boot[10,5] <- twoTailp(Protein_unchange_CORUMvsNoCORUM)

boot$FDR <- p.adjust(boot$pValue, method="BH")

write.table(boot, file=paste0("/Users/pc2644/Desktop/CPTAC_",cancer,"_bootstrape_CORUMvsNoCORUM.txt"), quote=F, sep="\t", row.names=F)

# ### add statistic analysis to compare CORUM and non-CORUM
# ### choose level from DNA, RNA and Protein
# # change <- "gain+"
# Ttest <- function(fc_gene, change) {
#   subData <- fc_gene[fc_gene$change2==change,]
#   subData_CORUM <- subData[subData$CORUM=="CORUM",]
#   subData_noCORUM <- subData[subData$CORUM=="non-CORUM",]
#   TDNA <- t.test(subData_CORUM$DNA_log2FC, subData_noCORUM$DNA_log2FC)
#   TRNA <- t.test(subData_CORUM$RNA_log2FC, subData_noCORUM$RNA_log2FC)
#   TProtein <- t.test(subData_CORUM$Protein_log2FC, subData_noCORUM$Protein_log2FC)
#   result <- c(change, TDNA$p.value, TRNA$p.value, TProtein$p.value)
#   return(result)
# }
# change <- rownames(Ttest_table)
# result_table <- data.frame(t(matrix(unlist(lapply(change, function(x) Ttest(fc_gene, x))),nrow=4)))

### export log2FC to compare among cancers
# fc_gene$cancer <- "colon2"
# fc_gene_colon2 <- fc_gene
# saveRDS(fc_gene_colon2, file = "/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/fc_gene_colon2.rds")

# fc_gene$cancer <- "breast2"
# fc_gene_breast2 <- fc_gene
# saveRDS(fc_gene_breast2, file = "/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/fc_gene_breast2.rds")

# fc_gene$cancer <- "ovarian2"
# fc_gene_ovarian2 <- fc_gene
# saveRDS(fc_gene_ovarian2, file = "/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/fc_gene_ovarian2.rds")

# fc_gene$cancer <- "ccrcc"
# fc_gene_ccrcc <- fc_gene
# saveRDS(fc_gene_ccrcc, file = "/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/fc_gene_ccrcc.rds")

# fc_gene$cancer <- "endometrail"
# fc_gene_endometrail <- fc_gene
# saveRDS(fc_gene_endometrail, file = "/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/fc_gene_endometrail.rds")

# fc_gene$cancer <- "hnscc"
# fc_gene_hnscc <- fc_gene
# saveRDS(fc_gene_hnscc, file = "/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/fc_gene_hnscc.rds")

# fc_gene$cancer <- "luad"
# fc_gene_luad <- fc_gene
# saveRDS(fc_gene_luad, file = "/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/fc_gene_luad.rds")

# summary <- fc_gene %>% 
#   group_by(change2, CORUM) %>%
#   summarise(n=n(),
#             DNA_log2FC=median(DNA_log2FC, na.rm=T),
#             RNA_log2FC=median(RNA_log2FC, na.rm=T),
#             Protein_log2FC=median(Protein_log2FC, na.rm=T))
# summary$cancer <- "colon2"


