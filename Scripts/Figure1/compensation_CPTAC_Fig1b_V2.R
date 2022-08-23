setwd("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC")

library(ggplot2)
library(dplyr)
library(circlize)
library(ComplexHeatmap)

rm(list=ls())

compensation <- data.frame(cancer=c("Colon", "Breast", "OV", "ccRCC", "Endometrial", "HNSC", "LUAD"),
                          
                           loss2_noCorum_RNA=rep(NA,7),
                           loss2_Corum_RNA=rep(NA,7),
                           loss_noCorum_RNA=rep(NA,7),
                           loss_Corum_RNA=rep(NA,7),
                          
                           gain_noCorum_RNA=rep(NA,7),
                           gain_Corum_RNA=rep(NA,7),
                           gain2_noCorum_RNA=rep(NA,7),
                           gain2_Corum_RNA=rep(NA,7),
                          
                           loss2_noCorum_Protein=rep(NA,7),
                           loss2_Corum_Protein=rep(NA,7),
                           loss_noCorum_Protein=rep(NA,7),
                           loss_Corum_Protein=rep(NA,7),
                          
                           gain_noCorum_Protein=rep(NA,7),
                           gain_Corum_Protein=rep(NA,7),
                           gain2_noCorum_Protein=rep(NA,7),
                           gain2_Corum_Protein=rep(NA,7))

load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/colon2_summary.RData")
summary_colon2 <- summary_fc_gene
compensation[1,-1] <- c(summary_fc_gene$RNA_compensation[c(7,8,5,6,1,2,3,4)], summary_fc_gene$Protein_compensation[c(7,8,5,6,1,2,3,4)]) * (-1)
rm(summary_fc_gene)

load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/breast2_summary.RData")
summary_breast2 <- summary_fc_gene
compensation[2,-1] <- c(summary_fc_gene$RNA_compensation[c(7,8,5,6,1,2,3,4)], summary_fc_gene$Protein_compensation[c(7,8,5,6,1,2,3,4)]) * (-1)
rm(summary_fc_gene)

load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/ovarian2_summary.RData")
summary_ovarian2 <- summary_fc_gene
compensation[3,-1] <- c(summary_fc_gene$RNA_compensation[c(7,8,5,6,1,2,3,4)], summary_fc_gene$Protein_compensation[c(7,8,5,6,1,2,3,4)]) * (-1)
rm(summary_fc_gene)

load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/ccrcc_summary.RData")
summary_ccrcc <- summary_fc_gene
compensation[4,-1] <- c(summary_fc_gene$RNA_compensation[c(7,8,5,6,1,2,3,4)], summary_fc_gene$Protein_compensation[c(7,8,5,6,1,2,3,4)]) * (-1)
rm(summary_fc_gene)

load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/endometrial_summary.RData")
summary_endometrial <- summary_fc_gene
compensation[5,-1] <- c(summary_fc_gene$RNA_compensation[c(7,8,5,6,1,2,3,4)], summary_fc_gene$Protein_compensation[c(7,8,5,6,1,2,3,4)]) * (-1)
rm(summary_fc_gene)

load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/hnscc_summary.RData")
summary_hnscc <- summary_fc_gene
compensation[6,-1] <- c(summary_fc_gene$RNA_compensation[c(7,8,5,6,1,2,3,4)], summary_fc_gene$Protein_compensation[c(7,8,5,6,1,2,3,4)]) * (-1)
rm(summary_fc_gene)

load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC/luad_summary.RData")
summary_luad <- summary_fc_gene
compensation[7,-1] <- c(summary_fc_gene$RNA_compensation[c(7,8,5,6,1,2,3,4)], summary_fc_gene$Protein_compensation[c(7,8,5,6,1,2,3,4)]) * (-1)
rm(summary_fc_gene)

rownames(compensation) <- compensation$cancer
compensation <- compensation[,-1]

### plot the heatmap
column_an<-HeatmapAnnotation(CORUM=rep(c("NoCorum", "Corum"),8),
                             CNV=rep(c(rep("Deep loss",2), rep("loss",2), rep("gain",2), rep("Profound gain",2)),2),
                             which="column",
                             col=list(CORUM=c("NoCorum"="#AAA900", "Corum"="#6500AA"),
                                      CNV=c("Deep loss"="#3300FF", "loss"="#ecdcff", "gain"="#ffd6d0", "Profound gain"="#FF0033")))

mycol <- colorRamp2(c(-0.4, 0, 0.1, 0.9), c("#c8c8c8", "#f6f6f6", "#b1eeec", "#00aaa9"))

pdf("/Users/pc2644/Desktop/CPTAC_choose genes_heatmap_dosage score.pdf", width=6, height=2.6)
Heatmap(as.matrix(compensation),
        col=mycol,
        name="Compensation Score",
        cluster_rows=TRUE,cluster_columns=FALSE,
        row_names_side="left",
        top_annotation=column_an,
        column_names_side="bottom",
        # right_annotation=row_an,
        column_split=factor(c(rep("RNA",8), rep("Protein",8)),levels=c("RNA","Protein")),
        rect_gp=gpar(col = "gray40",lwd=0.3),
        row_names_gp=gpar(fontsize = 8),
        show_column_names = FALSE,
        row_names_max_width=max_text_width(rownames(compensation), gp = gpar(fontsize = 12))
)
dev.off()



