setwd("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/log2FC")

library(ggplot2)
library(dplyr)
library(circlize)
library(ComplexHeatmap)

rm(list=ls())

fc_gene_colon2 <- readRDS("fc_gene_colon2.rds")
summary_colon2 <- fc_gene_colon2 %>% 
  group_by(change2, CORUM) %>%
  summarise(n=n(),
            # DNA_log2FC_25=quantile(DNA_log2FC,1/4,na.rm=T),
            # DNA_log2FC_75=quantile(DNA_log2FC,3/4,na.rm=T),
            # RNA_log2FC_25=quantile(RNA_log2FC,1/4,na.rm=T),
            # RNA_log2FC_75=quantile(RNA_log2FC,3/4,na.rm=T),
            # Protein_log2FC_25=quantile(Protein_log2FC,1/4,na.rm=T),
            # Protein_log2FC_75=quantile(Protein_log2FC,3/4,na.rm=T),
            DNA_log2FC=median(DNA_log2FC, na.rm=T),
            RNA_log2FC=median(RNA_log2FC, na.rm=T),
            Protein_log2FC=median(Protein_log2FC, na.rm=T))
summary_colon2$cancer <- "colon2"
summary_colon2$ratio <- summary_colon2$n/sum(summary_colon2$n)*100

fc_gene_breast2 <- readRDS("fc_gene_breast2.rds")
summary_breast2 <- fc_gene_breast2 %>% 
  group_by(change2, CORUM) %>%
  summarise(n=n(),
            # DNA_log2FC_25=quantile(DNA_log2FC,1/4,na.rm=T),
            # DNA_log2FC_75=quantile(DNA_log2FC,3/4,na.rm=T),
            # RNA_log2FC_25=quantile(RNA_log2FC,1/4,na.rm=T),
            # RNA_log2FC_75=quantile(RNA_log2FC,3/4,na.rm=T),
            # Protein_log2FC_25=quantile(Protein_log2FC,1/4,na.rm=T),
            # Protein_log2FC_75=quantile(Protein_log2FC,3/4,na.rm=T),
            DNA_log2FC=median(DNA_log2FC, na.rm=T),
            RNA_log2FC=median(RNA_log2FC, na.rm=T),
            Protein_log2FC=median(Protein_log2FC, na.rm=T))
summary_breast2$cancer <- "breast2"
summary_breast2$ratio <- summary_breast2$n/sum(summary_breast2$n)*100

fc_gene_ovarian2 <- readRDS("fc_gene_ovarian2.rds")
summary_ovarian2 <- fc_gene_ovarian2 %>% 
  group_by(change2, CORUM) %>%
  summarise(n=n(),
            # DNA_log2FC_25=quantile(DNA_log2FC,1/4,na.rm=T),
            # DNA_log2FC_75=quantile(DNA_log2FC,3/4,na.rm=T),
            # RNA_log2FC_25=quantile(RNA_log2FC,1/4,na.rm=T),
            # RNA_log2FC_75=quantile(RNA_log2FC,3/4,na.rm=T),
            # Protein_log2FC_25=quantile(Protein_log2FC,1/4,na.rm=T),
            # Protein_log2FC_75=quantile(Protein_log2FC,3/4,na.rm=T),
            DNA_log2FC=median(DNA_log2FC, na.rm=T),
            RNA_log2FC=median(RNA_log2FC, na.rm=T),
            Protein_log2FC=median(Protein_log2FC, na.rm=T))
summary_ovarian2$cancer <- "ovarian2"
summary_ovarian2$ratio <- summary_ovarian2$n/sum(summary_ovarian2$n)*100

fc_gene_ccrcc <- readRDS("fc_gene_ccrcc.rds")
summary_ccrcc <- fc_gene_ccrcc %>% 
  group_by(change2, CORUM) %>%
  summarise(n=n(),
            # DNA_log2FC_25=quantile(DNA_log2FC,1/4,na.rm=T),
            # DNA_log2FC_75=quantile(DNA_log2FC,3/4,na.rm=T),
            # RNA_log2FC_25=quantile(RNA_log2FC,1/4,na.rm=T),
            # RNA_log2FC_75=quantile(RNA_log2FC,3/4,na.rm=T),
            # Protein_log2FC_25=quantile(Protein_log2FC,1/4,na.rm=T),
            # Protein_log2FC_75=quantile(Protein_log2FC,3/4,na.rm=T),
            DNA_log2FC=median(DNA_log2FC, na.rm=T),
            RNA_log2FC=median(RNA_log2FC, na.rm=T),
            Protein_log2FC=median(Protein_log2FC, na.rm=T))
summary_ccrcc$cancer <- "ccrcc"
summary_ccrcc$ratio <- summary_ccrcc$n/sum(summary_ccrcc$n)*100

fc_gene_endometrail <- readRDS("fc_gene_endometrail.rds")
summary_endometrail <- fc_gene_endometrail %>% 
  group_by(change2, CORUM) %>%
  summarise(n=n(),
            # DNA_log2FC_25=quantile(DNA_log2FC,1/4,na.rm=T),
            # DNA_log2FC_75=quantile(DNA_log2FC,3/4,na.rm=T),
            # RNA_log2FC_25=quantile(RNA_log2FC,1/4,na.rm=T),
            # RNA_log2FC_75=quantile(RNA_log2FC,3/4,na.rm=T),
            # Protein_log2FC_25=quantile(Protein_log2FC,1/4,na.rm=T),
            # Protein_log2FC_75=quantile(Protein_log2FC,3/4,na.rm=T),
            DNA_log2FC=median(DNA_log2FC, na.rm=T),
            RNA_log2FC=median(RNA_log2FC, na.rm=T),
            Protein_log2FC=median(Protein_log2FC, na.rm=T))
summary_endometrail$cancer <- "endometrail"
summary_endometrail$ratio <- summary_endometrail$n/sum(summary_endometrail$n)*100

fc_gene_hnscc <- readRDS("fc_gene_hnscc.rds")
summary_hnscc <- fc_gene_hnscc %>% 
  group_by(change2, CORUM) %>%
  summarise(n=n(),
            # DNA_log2FC_25=quantile(DNA_log2FC,1/4,na.rm=T),
            # DNA_log2FC_75=quantile(DNA_log2FC,3/4,na.rm=T),
            # RNA_log2FC_25=quantile(RNA_log2FC,1/4,na.rm=T),
            # RNA_log2FC_75=quantile(RNA_log2FC,3/4,na.rm=T),
            # Protein_log2FC_25=quantile(Protein_log2FC,1/4,na.rm=T),
            # Protein_log2FC_75=quantile(Protein_log2FC,3/4,na.rm=T),
            DNA_log2FC=median(DNA_log2FC, na.rm=T),
            RNA_log2FC=median(RNA_log2FC, na.rm=T),
            Protein_log2FC=median(Protein_log2FC, na.rm=T))
summary_hnscc$cancer <- "hnscc"
summary_hnscc$ratio <- summary_hnscc$n/sum(summary_hnscc$n)*100

fc_gene_luad <- readRDS("fc_gene_luad.rds")
summary_luad <- fc_gene_luad %>% 
  group_by(change2, CORUM) %>%
  summarise(n=n(),
            # DNA_log2FC_25=quantile(DNA_log2FC,1/4,na.rm=T),
            # DNA_log2FC_75=quantile(DNA_log2FC,3/4,na.rm=T),
            # RNA_log2FC_25=quantile(RNA_log2FC,1/4,na.rm=T),
            # RNA_log2FC_75=quantile(RNA_log2FC,3/4,na.rm=T),
            # Protein_log2FC_25=quantile(Protein_log2FC,1/4,na.rm=T),
            # Protein_log2FC_75=quantile(Protein_log2FC,3/4,na.rm=T),
            DNA_log2FC=median(DNA_log2FC, na.rm=T),
            RNA_log2FC=median(RNA_log2FC, na.rm=T),
            Protein_log2FC=median(Protein_log2FC, na.rm=T))
summary_luad$cancer <- "luad"
summary_luad$ratio <- summary_luad$n/sum(summary_luad$n)*100

summary_pancan <- rbind(summary_colon2, summary_breast2, summary_ovarian2, summary_ccrcc, summary_endometrail, summary_hnscc, summary_luad)

### the number of genes in each categroies? (fig s1a)
summary_pancan_plot1 <- summary_pancan[summary_pancan$change2!="unchange",]
summary_pancan_plot1$cancer <- factor(summary_pancan_plot1$cancer, levels=unique(summary_pancan_plot1$cancer))
# pdf("/Users/pc2644/Desktop/CPTAC_N_log2FC.pdf", width=17, height=5)
ggplot(summary_pancan_plot1, aes(x=change2, y=ratio, fill=CORUM)) + geom_bar(position="stack", stat="identity") + facet_wrap(~cancer, ncol=7)
# dev.off()

### log2DNA dirstribute
fc_gene <- rbind(fc_gene_colon2, fc_gene_breast2, fc_gene_ovarian2, fc_gene_ccrcc, fc_gene_endometrail, fc_gene_hnscc, fc_gene_luad)

fc_gene_plot1 <- fc_gene[fc_gene$change2!="unchange",]
ggplot(fc_gene_plot1, aes(DNA_log2FC, stat(density))) +
    geom_density(alpha=0.1, size = 1, aes(color=cancer)) +
    xlab("DNA Log2FC") +
    coord_cartesian(xlim=c(-2,2))

### select genes of median at 1 or -0.7 for each cancer
# data <- fc_gene_ccrcc
# median <- -0.6
# n <- 5000
select_gene <- function(data, median, n) {
  data <- data[order(data$DNA_log2FC, decreasing=T),]
  data <- na.omit(data)
  index <- data$DNA_log2FC>median
  index <- grep(FALSE, index)[1]
  if (median>0 & index<n/2) (stop("there is not enough data to meet n"))
  if (median<0 & index>(nrow(data)-n/2)) (stop("there is not enough data to meet n"))
  data_select <- data[(index-round(n/2)):(index+round(n/2)-1),]
}

fc_gene_plot2 <- rbind(select_gene(fc_gene_colon2, 0.7, 5000), 
                       select_gene(fc_gene_colon2, -0.7, 5000),
                       select_gene(fc_gene_breast2, 0.7, 5000), 
                       select_gene(fc_gene_breast2, -0.7, 5000),
                       select_gene(fc_gene_ovarian2, 0.7, 5000), 
                       select_gene(fc_gene_ovarian2, -0.7, 5000),
                       select_gene(fc_gene_ccrcc, 0.7, 1000), 
                       select_gene(fc_gene_ccrcc, -0.7, 1000),
                       select_gene(fc_gene_endometrail, 0.7, 5000), 
                       select_gene(fc_gene_endometrail, -0.7, 1000),
                       select_gene(fc_gene_hnscc, 0.7, 5000), 
                       select_gene(fc_gene_hnscc, -0.7, 2000),
                       select_gene(fc_gene_luad, 0.7, 5000), 
                       select_gene(fc_gene_luad, -0.7, 900))

fc_gene_plot2$cancer <- factor(fc_gene_plot2$cancer, levels=unique(fc_gene_plot2$cancer))

### add to calculate compensation score at gene level
fc_gene_plot2$compensation_RNA <- rep(NA,nrow(fc_gene_plot2))
fc_gene_plot2$compensation_Protein <- rep(NA,nrow(fc_gene_plot2))
fc_gene_plot2$compensation_RNA[fc_gene_plot2$change=="gain"] <- fc_gene_plot2$RNA_log2FC[fc_gene_plot2$change=="gain"]-fc_gene_plot2$DNA_log2FC[fc_gene_plot2$change=="gain"]
fc_gene_plot2$compensation_Protein[fc_gene_plot2$change=="gain"] <- fc_gene_plot2$Protein_log2FC[fc_gene_plot2$change=="gain"]-fc_gene_plot2$DNA_log2FC[fc_gene_plot2$change=="gain"]
fc_gene_plot2$compensation_RNA[fc_gene_plot2$change=="loss"] <- (fc_gene_plot2$RNA_log2FC[fc_gene_plot2$change=="loss"]-fc_gene_plot2$DNA_log2FC[fc_gene_plot2$change=="loss"])*(-1)
fc_gene_plot2$compensation_Protein[fc_gene_plot2$change=="loss"] <- (fc_gene_plot2$Protein_log2FC[fc_gene_plot2$change=="loss"]-fc_gene_plot2$DNA_log2FC[fc_gene_plot2$change=="loss"])*-1
fc_gene_plot_add <- fc_gene_plot2[fc_gene_plot2$change!="unchange",-2:-4]
fc_gene_plot_add <- reshape2::melt(fc_gene_plot_add, id=c("gene","CORUM","change","change2","cancer"))
fc_gene_plot_add$change2 <- factor(fc_gene_plot_add$change2, levels=c("loss+","loss", "unchange", "gain","gain+"))
ggplot(fc_gene_plot_add, aes(x=change, y=value, fill=CORUM)) +
  geom_boxplot(alpha=0.7, outlier.shape = NA) +
  coord_cartesian(ylim=c(-2,2)) +
  facet_wrap(~variable+cancer, ncol=2) +
  geom_hline(yintercept=0, color="blue") +
  labs(x="genes classfied based on CNV", y="Compensation Score")

###To be done, check if there is some gene/patient override

# length(unique(fc_gene_plot2$gene[(fc_gene_plot2$cancer=="colon2" & fc_gene_plot2$change=="gain")]))


# summary_fc_gene_plot2 <- fc_gene_plot2 %>% 
#   group_by(cancer, change) %>%
#   summarise(n=n(),
#             quantile=quantile(DNA_log2FC,0.75))
# summary_fc_gene_plot2$cor <- ifelse(summary_fc_gene_plot2$quantile>0, summary_fc_gene_plot2$quantile+0.3, summary_fc_gene_plot2$quantile-0.3)

### fig. S1b choose genes of certain median of gain/loss
# pdf("/Users/pc2644/Desktop/CPTAC_choose genes.pdf", width=5, height=4)
ggplot(fc_gene_plot2, aes(x=change, y=DNA_log2FC, fill=cancer)) +
  geom_boxplot(alpha=0.7) +
  coord_cartesian(ylim=c(-2,2)) +
  labs(y="log2FC")
# dev.off()

### calculate the RNA and Protein log2FC of select genes
fc_gene_plot3 <- fc_gene_plot2 %>% 
  group_by(cancer, change, CORUM) %>%
  summarise(n=n(),
            DNA_log2FC=median(DNA_log2FC, na.rm=T),
            RNA_log2FC=median(RNA_log2FC, na.rm=T),
            Protein_log2FC=median(Protein_log2FC, na.rm=T))
fc_gene_plot3_2 <- data.frame(DNA_loss_nonCorum=fc_gene_plot3$DNA_log2FC[fc_gene_plot3$change=="loss" & fc_gene_plot3$CORUM=="non-CORUM"],
                              DNA_loss_Corum=fc_gene_plot3$DNA_log2FC[fc_gene_plot3$change=="loss" & fc_gene_plot3$CORUM=="CORUM"],
                              RNA_loss_nonCorum=fc_gene_plot3$RNA_log2FC[fc_gene_plot3$change=="loss" & fc_gene_plot3$CORUM=="non-CORUM"],
                              RNA_loss_Corum=fc_gene_plot3$RNA_log2FC[fc_gene_plot3$change=="loss" & fc_gene_plot3$CORUM=="CORUM"],
                              Protein_loss_nonCorum=fc_gene_plot3$Protein_log2FC[fc_gene_plot3$change=="loss" & fc_gene_plot3$CORUM=="non-CORUM"],
                              Protein_loss_Corum=fc_gene_plot3$Protein_log2FC[fc_gene_plot3$change=="loss" & fc_gene_plot3$CORUM=="CORUM"],
                              DNA_gain_nonCorum=fc_gene_plot3$DNA_log2FC[fc_gene_plot3$change=="gain" & fc_gene_plot3$CORUM=="non-CORUM"],
                              DNA_gain_Corum=fc_gene_plot3$DNA_log2FC[fc_gene_plot3$change=="gain" & fc_gene_plot3$CORUM=="CORUM"],
                              RNA_gain_nonCorum=fc_gene_plot3$RNA_log2FC[fc_gene_plot3$change=="gain" & fc_gene_plot3$CORUM=="non-CORUM"],
                              RNA_gain_Corum=fc_gene_plot3$RNA_log2FC[fc_gene_plot3$change=="gain" & fc_gene_plot3$CORUM=="CORUM"],
                              Protein_gain_nonCorum=fc_gene_plot3$Protein_log2FC[fc_gene_plot3$change=="gain" & fc_gene_plot3$CORUM=="non-CORUM"],
                              Protein_gain_Corum=fc_gene_plot3$Protein_log2FC[fc_gene_plot3$change=="gain" & fc_gene_plot3$CORUM=="CORUM"]
)
rownames(fc_gene_plot3_2) <- unique(fc_gene_plot3$cancer)

### fig S1c: RNA Protein log2FC of selected genes (heatmap)
column_an<-HeatmapAnnotation(CORUM=rep(c("Non-CORUM", "CORUM"),6),
                             # Level=rep(c("RNA","RNA","Protein","Protein"),4),
                             CNV=c(rep("loss",6), rep("gain",6)),
                             which="column",
                             col=list(CORUM=c("Non-CORUM"="#7a5195", "CORUM"="#ffa600"),
                                      # Level=c("RNA"=#003f5c","Protein"="#ef5675"), 
                                      CNV=c("loss"="#a1c292", "gain"="#ef9fa0")))

mycol <- colorRamp2(c(-1.15,0,1.15), c("blue","white","red"))
# pdf("/Users/pc2644/Desktop/CPTAC_choose genes_heatmap.pdf", width=6, height=6)
Heatmap(as.matrix(fc_gene_plot3_2),col=mycol,name="log2FC median",
        cluster_rows=TRUE,cluster_columns=FALSE,
        row_names_side="left",
        top_annotation=column_an,
        column_names_side="bottom",
        # right_annotation=row_an,
        column_split=factor(rep(c("DNA","DNA", "RNA","RNA","Protein","Protein"),2),levels=c("DNA","RNA","Protein")),
        rect_gp=gpar(col = "gray40",lwd=1),
        row_names_gp=gpar(fontsize = 8),
        column_names_gp=gpar(fontsize = 8),
        row_names_max_width=max_text_width(rownames(fc_gene_plot3_2), gp = gpar(fontsize = 12))
)
# dev.off()

### calculate the compensation score
fc_gene_plot3_3 <- fc_gene_plot3_2
fc_gene_plot3_3[,3:4] <- (fc_gene_plot3_3[,3:4]-fc_gene_plot3_3[,1:2])*(-1)
fc_gene_plot3_3[,5:6] <- (fc_gene_plot3_3[,5:6]-fc_gene_plot3_3[,1:2])*(-1)
fc_gene_plot3_3[,1:2] <- (fc_gene_plot3_3[,1:2]-fc_gene_plot3_3[,1:2])*(-1)
fc_gene_plot3_3[,9:10] <- (fc_gene_plot3_3[,9:10]-fc_gene_plot3_3[,7:8])
fc_gene_plot3_3[,11:12] <- (fc_gene_plot3_3[,11:12]-fc_gene_plot3_3[,7:8])
fc_gene_plot3_3[,7:8] <- (fc_gene_plot3_3[,7:8]-fc_gene_plot3_3[,7:8])
fc_gene_plot3_3 <- fc_gene_plot3_3[,c(-1:-2,-7:-8)]

###fig 1c: amp-compensation score
### fig S1c: RNA Protein log2FC of selected genes (heatmap)
column_an<-HeatmapAnnotation(CORUM=rep(c("Non-CORUM", "CORUM"),4),
                             # Level=rep(c("RNA","RNA","Protein","Protein"),4),
                             CNV=c(rep("loss",4), rep("gain",4)),
                             which="column",
                             col=list(CORUM=c("Non-CORUM"="#7a5195", "CORUM"="#ffa600"),
                                      # Level=c("RNA"=#003f5c","Protein"="#ef5675"), 
                                      CNV=c("loss"="#a1c292", "gain"="#ef9fa0")))

mycol <- colorRamp2(c(-0.6, 0, 0.6), c("blue","white","red"))
# pdf("/Users/pc2644/Desktop/CPTAC_choose genes_heatmap_dosage score.pdf", width=6, height=6)
Heatmap(as.matrix(fc_gene_plot3_3),col=mycol,name="Dosage Score",
        cluster_rows=TRUE,cluster_columns=FALSE,
        row_names_side="left",
        top_annotation=column_an,
        column_names_side="bottom",
        # right_annotation=row_an,
        column_split=factor(rep(c("RNA","RNA","Protein","Protein"),2),levels=c("RNA","Protein")),
        rect_gp=gpar(col = "gray40",lwd=1),
        row_names_gp=gpar(fontsize = 12),
        column_names_gp=gpar(fontsize = 8),
        row_names_max_width=max_text_width(rownames(fc_gene_plot3_3), gp = gpar(fontsize = 12))
)
# dev.off()



summary_pancan_plot3 <- data.frame(DNA_loss2_nonCorum=summary_pancan_plot1$DNA_log2FC[summary_pancan_plot1$change2=="loss+" & summary_pancan_plot1$CORUM=="non-CORUM"],
                                   DNA_loss2_Corum=summary_pancan_plot1$DNA_log2FC[summary_pancan_plot1$change2=="loss+" & summary_pancan_plot1$CORUM=="CORUM"],
                                   DNA_loss_nonCorum=summary_pancan_plot1$DNA_log2FC[summary_pancan_plot1$change2=="loss" & summary_pancan_plot1$CORUM=="non-CORUM"],
                                   DNA_loss_Corum=summary_pancan_plot1$DNA_log2FC[summary_pancan_plot1$change2=="loss" & summary_pancan_plot1$CORUM=="CORUM"],
                                   DNA_gain_nonCorum=summary_pancan_plot1$DNA_log2FC[summary_pancan_plot1$change2=="gain" & summary_pancan_plot1$CORUM=="non-CORUM"],
                                   DNA_gain_Corum=summary_pancan_plot1$DNA_log2FC[summary_pancan_plot1$change2=="gain" & summary_pancan_plot1$CORUM=="CORUM"],
                                   DNA_gain2_nonCorum=summary_pancan_plot1$DNA_log2FC[summary_pancan_plot1$change2=="gain+" & summary_pancan_plot1$CORUM=="non-CORUM"],
                                   DNA_gain2_Corum=summary_pancan_plot1$DNA_log2FC[summary_pancan_plot1$change2=="gain+" & summary_pancan_plot1$CORUM=="CORUM"]
)







