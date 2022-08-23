setwd("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC")

rm(list=ls())

library(ComplexHeatmap)
library(ggplot2)
library(qusage)
library(patchwork)
library(VennDiagram)
library(circlize)

### GSEA results (same sets of genes)
COAD_protein_path <- "/Users/pc2644/gsea_home/output/oct13/COAD_lm_Protein_AS_C5BP.GseaPreranked.1634172821100/edb/results.edb"
BRCA_protein_path <- "/Users/pc2644/gsea_home/output/oct13/BRCA_lm_Protein_AS_C5BP.GseaPreranked.1634172834859/edb/results.edb"
OV_protein_path <- "/Users/pc2644/gsea_home/output/oct13/OV_lm_Protein_AS_C5BP.GseaPreranked.1634172850914/edb/results.edb"
ccRCC_protein_path <- "/Users/pc2644/gsea_home/output/oct13/ccRCC_lm_Protein_AS_C5BP.GseaPreranked.1634173279023/edb/results.edb"
UCEC_protein_path <- "/Users/pc2644/gsea_home/output/oct13/UCEC_lm_Protein_AS_C5BP.GseaPreranked.1634173819286/edb/results.edb"
HNSC_protein_path <- "/Users/pc2644/gsea_home/output/oct13/HNSC_lm_Protein_AS_C5BP.GseaPreranked.1634174119369/edb/results.edb"
LUAD_protein_path <- "/Users/pc2644/gsea_home/output/oct13/LUAD_lm_Protein_AS_C5BP.GseaPreranked.1634174355923/edb/results.edb"

COAD_RNA_path <- "/Users/pc2644/gsea_home/output/oct13/COAD_lm_RNA_AS_C5BP.GseaPreranked.1634172812611/edb/results.edb"
BRCA_RNA_path <- "/Users/pc2644/gsea_home/output/oct13/BRCA_lm_RNA_AS_C5BP.GseaPreranked.1634172840859/edb/results.edb"
OV_RNA_path <- "/Users/pc2644/gsea_home/output/oct13/OV_lm_RNA_AS_C5BP.GseaPreranked.1634174831765/edb/results.edb"
ccRCC_RNA_path <- "/Users/pc2644/gsea_home/output/oct13/ccRCC_lm_RNA_AS_C5BP.GseaPreranked.1634174826522/edb/results.edb"
UCEC_RNA_path <- "/Users/pc2644/gsea_home/output/oct13/UCEC_lm_RNA_AS_C5BP.GseaPreranked.1634173806787/edb/results.edb"
HNSC_RNA_path <- "/Users/pc2644/gsea_home/output/oct13/HNSC_lm_RNA_AS_C5BP.GseaPreranked.1634174101379/edb/results.edb"
LUAD_RNA_path <- "/Users/pc2644/gsea_home/output/oct13/LUAD_lm_RNA_AS_C5BP.GseaPreranked.1634174364383/edb/results.edb"

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

COAD_protein <- extractGeneSets(COAD_protein_path)
COAD_protein <- COAD_protein[,c(1,3,4,5)]
colnames(COAD_protein)[2:4] <- paste0(c("NES","NP","FDR"),"_COAD_protein")

BRCA_protein <- extractGeneSets(BRCA_protein_path)
BRCA_protein <- BRCA_protein[,c(1,3,4,5)]
colnames(BRCA_protein)[2:4] <- paste0(c("NES","NP","FDR"),"_BRCA_protein")

OV_protein <- extractGeneSets(OV_protein_path)
OV_protein <- OV_protein[,c(1,3,4,5)]
colnames(OV_protein)[2:4] <- paste0(c("NES","NP","FDR"),"_OV_protein")

ccRCC_protein <- extractGeneSets(ccRCC_protein_path)
ccRCC_protein <- ccRCC_protein[,c(1,3,4,5)]
colnames(ccRCC_protein)[2:4] <- paste0(c("NES","NP","FDR"),"_ccRCC_protein")

UCEC_protein <- extractGeneSets(UCEC_protein_path)
UCEC_protein <- UCEC_protein[,c(1,3,4,5)]
colnames(UCEC_protein)[2:4] <- paste0(c("NES","NP","FDR"),"_UCEC_protein")

HNSC_protein <- extractGeneSets(HNSC_protein_path)
HNSC_protein <- HNSC_protein[,c(1,3,4,5)]
colnames(HNSC_protein)[2:4] <- paste0(c("NES","NP","FDR"),"_HNSC_protein")

LUAD_protein <- extractGeneSets(LUAD_protein_path)
LUAD_protein <- LUAD_protein[,c(1,3,4,5)]
colnames(LUAD_protein)[2:4] <- paste0(c("NES","NP","FDR"),"_LUAD_protein")

COAD_RNA <- extractGeneSets(COAD_RNA_path)
COAD_RNA <- COAD_RNA[,c(1,3,4,5)]
colnames(COAD_RNA)[2:4] <- paste0(c("NES","NP","FDR"),"_COAD_RNA")

BRCA_RNA <- extractGeneSets(BRCA_RNA_path)
BRCA_RNA <- BRCA_RNA[,c(1,3,4,5)]
colnames(BRCA_RNA)[2:4] <- paste0(c("NES","NP","FDR"),"_BRCA_RNA")

OV_RNA <- extractGeneSets(OV_RNA_path)
OV_RNA <- OV_RNA[,c(1,3,4,5)]
colnames(OV_RNA)[2:4] <- paste0(c("NES","NP","FDR"),"_OV_RNA")

ccRCC_RNA <- extractGeneSets(ccRCC_RNA_path)
ccRCC_RNA <- ccRCC_RNA[,c(1,3,4,5)]
colnames(ccRCC_RNA)[2:4] <- paste0(c("NES","NP","FDR"),"_ccRCC_RNA")

UCEC_RNA <- extractGeneSets(UCEC_RNA_path)
UCEC_RNA <- UCEC_RNA[,c(1,3,4,5)]
colnames(UCEC_RNA)[2:4] <- paste0(c("NES","NP","FDR"),"_UCEC_RNA")

HNSC_RNA <- extractGeneSets(HNSC_RNA_path)
HNSC_RNA <- HNSC_RNA[,c(1,3,4,5)]
colnames(HNSC_RNA)[2:4] <- paste0(c("NES","NP","FDR"),"_HNSC_RNA")

LUAD_RNA <- extractGeneSets(LUAD_RNA_path)
LUAD_RNA <- LUAD_RNA[,c(1,3,4,5)]
colnames(LUAD_RNA)[2:4] <- paste0(c("NES","NP","FDR"),"_LUAD_RNA")

### combine GSEA results
GSEA_Protein <- merge(COAD_protein, BRCA_protein, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, OV_protein, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, ccRCC_protein, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, UCEC_protein, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, HNSC_protein, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, LUAD_protein, by.x="GeneSet",by.y="GeneSet",all=T)

GSEA_Protein <- merge(GSEA_Protein, COAD_RNA, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, BRCA_RNA, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, OV_RNA, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, ccRCC_RNA, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, UCEC_RNA, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, HNSC_RNA, by.x="GeneSet",by.y="GeneSet",all=T)
GSEA_Protein <- merge(GSEA_Protein, LUAD_RNA, by.x="GeneSet",by.y="GeneSet",all=T)

write.table(GSEA_Protein, file="CPTAC_individual_Protein_RNA_commonGenes_C5BP_20211013.txt", quote=FALSE, sep="\t", row.names=F, col.names=T, na="")

### update module 1: compare the RNA and Protein enrichment

# RNA_up <- GSEA_Protein$GeneSet[!is.na(GSEA_Protein$FDR_CPTAC_RNA_AS_cancer) & GSEA_Protein$FDR_CPTAC_RNA_AS_cancer<0.1 & GSEA_Protein$NES_CPTAC_RNA_AS_cancer>0]
# RNA_down <- GSEA_Protein$GeneSet[!is.na(GSEA_Protein$FDR_CPTAC_RNA_AS_cancer) & GSEA_Protein$FDR_CPTAC_RNA_AS_cancer<0.01 & GSEA_Protein$NES_CPTAC_RNA_AS_cancer<0]
# protein_up <- GSEA_Protein$GeneSet[!is.na(GSEA_Protein$FDR_CPTAC_AS_cancer) & GSEA_Protein$FDR_CPTAC_AS_cancer<0.0001 & GSEA_Protein$NES_CPTAC_AS_cancer>0]
# protein_down <- GSEA_Protein$GeneSet[!is.na(GSEA_Protein$FDR_CPTAC_AS_cancer) & GSEA_Protein$FDR_CPTAC_AS_cancer<0.001 & GSEA_Protein$NES_CPTAC_AS_cancer<0]

rownames(GSEA_Protein) <- GSEA_Protein$GeneSet
GSEA_Protein <- GSEA_Protein[,-1]

# ## venn plot
# venn.diagram(list(RNA_up, protein_up),
#              category.names = c("Enriched in high aneuploidy (RNA)" , "Enriched in high aneuploidy (Protein)"),
#              filename = 'upreg_RNA_Protein_GSEA_CPTAC',
#              output=TRUE,
#              imagetype="png",
#              height = 3000,
#              width = 3000,
#              resolution = 300)
# venn.diagram(list(RNA_down, protein_down),
#              category.names = c("Enriched in low aneuploidy (RNA)" , "Enriched in low aneuploidy (Protein)"),
#              filename = 'downreg_RNA_Protein_GSEA_CPTAC',
#              output=TRUE,
#              imagetype="png",
#              height = 3000,
#              width = 3000,
#              resolution = 300)

up_geneset <- c("GO_DNA_REPLICATION",
                "GO_DNA_METABOLIC_PROCESS",
                "GO_CHROMATIN_ORGANIZATION",
                "GO_CHROMATIN_REMODELING",
                
                "GO_MRNA_METABOLIC_PROCESS",
                "GO_MRNA_3_END_PROCESSING",
                "GO_MRNA_TRANSPORT",
                "GO_RNA_SPLICING",
                "GO_SPLICEOSOMAL_SNRNP_ASSEMBLY",
                "GO_RNA_POLYADENYLATION",
                "GO_RNA_METHYLATION",
                "GO_DNA_TEMPLATED_TRANSCRIPTION_ELONGATION",
                "GO_DNA_TEMPLATED_TRANSCRIPTION_TERMINATION",
                
                "GO_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS",
                "GO_RIBOSOME_BIOGENESIS",
                "GO_NCRNA_METABOLIC_PROCESS",
                "GO_NCRNA_PROCESSING",
                "GO_NCRNA_TRANSCRIPTION",
                "GO_RRNA_METABOLIC_PROCESS",
                "GO_TRNA_METABOLIC_PROCESS",
                "GO_TRNA_PROCESSING",
                
                "GO_MITOCHONDRIAL_TRANSLATION",
                "GO_MITOCHONDRIAL_TRANSPORT",
                "GO_MITOCHONDRIAL_CALCIUM_ION_TRANSMEMBRANE_TRANSPORT",  
                
                "GO_ACTIVATION_OF_IMMUNE_RESPONSE",
                "GO_INNATE_IMMUNE_RESPONSE",
                "GO_ADAPTIVE_IMMUNE_RESPONSE",
                "GO_REGULATION_OF_CELL_KILLING",
                "GO_RESPONSE_TO_CYTOKINE",
                "GO_RESPONSE_TO_INTERFERON_GAMMA",
                
                "GO_ACTIN_FILAMENT_POLYMERIZATION",
                "GO_ACTIN_FILAMENT_DEPOLYMERIZATION",
                "GO_ACTIN_FILAMENT_ORGANIZATION",
                "GO_CELL_CELL_ADHESION",
                "GO_CELL_MATRIX_ADHESION",
                "GO_CELL_MOTILITY"
)

up_geneset <- up_geneset[up_geneset %in% rownames(GSEA_NES)]

process <- c(rep("DNA",4),
             rep("Transcription",9),
             rep("Translation",8),
             rep("Mitochondria",3),
             rep("Immune response",5),
             rep("Cytoskeleton",5))

process <- factor(process, levels=c("DNA","Transcription", "Translation", "Mitochondria", "Immune response", "Cytoskeleton"))

GSEA_FDR <- GSEA_Protein[seq(3, ncol(GSEA_Protein),3)]
GSEA_NES <- GSEA_Protein[seq(1, ncol(GSEA_Protein),3)]

GSEA_NES_pick <- GSEA_NES[up_geneset,]
GSEA_FDR_pick <- GSEA_FDR[up_geneset,]

up_geneset <- factor(up_geneset, levels=up_geneset)

GSEA_NES_pick[GSEA_FDR_pick>=0.1 & !is.na(GSEA_FDR_pick)] <- NA

row_an<-HeatmapAnnotation(pathway=process,
                          show_legend=T,
                          show_annotation_name=c(pathway = FALSE),
                          annotation_label=NULL,
                          which="row")

column_an<-HeatmapAnnotation(sample=c(rep("Protein",7),rep("RNA",7)),
                             which="column",
                             show_annotation_name=c(sample = FALSE),
                             col=list(sample=c("Protein"="#CC6677","RNA"="#88CCEE")))

mycol <- colorRamp2(c(-3.5,0,3.5), c("blue","white","red"))

pdf("/Users/pc2644/Desktop/selected_geneset_individual_cancer_CPTAC_Protein_RNA_Heatmap.pdf", width=8.5, height=7)
Heatmap(as.matrix(GSEA_NES_pick),col=mycol,name="GSEA NES",
        cluster_rows=FALSE,cluster_columns=FALSE,
        row_names_side="left",top_annotation=column_an,
        column_names_side="bottom",right_annotation=row_an,
        column_split=factor(c(rep("Protein",7),rep("RNA",7)),levels=c("Protein","RNA")),
        row_split=process,
        rect_gp=gpar(col = "gray40",lwd=0.3),
        row_names_rot=0,
        row_names_gp=gpar(fontsize = 8),
        column_names_gp=gpar(fontsize = 8),
        row_names_max_width=max_text_width(rownames(GSEA_NES_pick), gp = gpar(fontsize = 12))
)
dev.off()


# GSEA_filter1 <- na.omit(GSEA_Protein)
# rownames(GSEA_filter1) <- GSEA_filter1$GeneSet
# GSEA_filter1 <- GSEA_filter1[,-1]
# p_threshold <- 0.01
# GSEA_filter1_p <- GSEA_filter1[seq(2, ncol(GSEA_filter1),3)]
# GSEA_filter1_NES <- GSEA_filter1[seq(1, ncol(GSEA_filter1),3)]
# GSEA_filter1_NES_threshold <- GSEA_filter1_NES
# GSEA_filter1_NES_threshold[GSEA_filter1_p>p_threshold] <- "NS"
# 
# index <- grepl("polymerase",rownames(GSEA_filter1_NES_threshold),ignore.case=T) |
#          grepl("rna",rownames(GSEA_filter1_NES_threshold),ignore.case=T) |
#          grepl("transcript",rownames(GSEA_filter1_NES_threshold),ignore.case=T) |
#          grepl("ribosome",rownames(GSEA_filter1_NES_threshold),ignore.case=T) |
#          grepl("peptide",rownames(GSEA_filter1_NES_threshold),ignore.case=T) |
#          grepl("protein",rownames(GSEA_filter1_NES_threshold),ignore.case=T) |
#          grepl("translation",rownames(GSEA_filter1_NES_threshold),ignore.case=T) |
#          grepl("proteasome",rownames(GSEA_filter1_NES_threshold),ignore.case=T)
# 
# GSEA_filter1_NES_threshold_pick <- GSEA_filter1_NES_threshold[index,]
# index2 <- unlist(apply(GSEA_filter1_NES_threshold_pick, 1, function(x) sum(x=="NS")!=length(x)))
# GSEA_filter1_NES_threshold_pick <- GSEA_filter1_NES_threshold_pick[index2,]
# write.table(GSEA_filter1_NES_threshold_pick, file="CCLE_CPTAC_pan_commonGenes_C2_filter_p0.01_keywords.txt", quote=FALSE, sep="\t", row.names=T, col.names=T, na="")  
  
### visualize the result
# result <- read.table(common.geneset.summary,sep="\t",header=TRUE,stringsAsFactors=FALSE)
# c5bp_geneset_path <- "/Users/chengpan/Desktop/Genesets/c5.bp.v7.0.symbols.gmt"
# c5bp_geneset <- read.gmt(c5bp_geneset_path)
# common_geneset23 <- c5bp_geneset[rownames(common.geneset.summary23)]
# common_geneset123 <- c5bp_geneset[rownames(common.geneset.summary123)]
# distance.geneset23 <- matrix(data=NA, nrow=length(common_geneset23), ncol=length(common_geneset23))
# for (i in 1:length(common_geneset23)) {
#   for (j in 1:length(common_geneset23)) {
#     distance.geneset23[i,j] <- length(intersect(common_geneset23[[i]],common_geneset23[[j]]))/length(union(common_geneset23[[i]],common_geneset23[[j]]))
#   }
# }
# rownames(distance.geneset23) <- rownames(common.geneset.summary23)
# colnames(distance.geneset23) <- rownames(common.geneset.summary23)
# h <- Heatmap(distance.geneset23,cluster_rows=TRUE,cluster_columns=TRUE)
# cluster <- column_dend(h)
# cluster.cut <- cutree(as.hclust(cluster),k=11)
# cluster.cut <- data.frame(geneset=names(cluster.cut),groups=cluster.cut)
# cluster.cut <- cluster.cut[order(cluster.cut$groups),]
# common.geneset.summary23.rank <- common.geneset.summary23[cluster.cut$geneset,]

# df <- data.frame(geneset=rownames(common.geneset.summary23),
#                  group=ifelse(grepl("MITOC", rownames(common.geneset.summary23)),1,5))
# df$group[df$group==1 & grepl("TRANSPORT", rownames(common.geneset.summary23))] <- 2
# df$group[df$group==1 & grepl("ORGANIZATION", rownames(common.geneset.summary23))] <- 3
# df$group[grepl("AMINO_ACID", rownames(common.geneset.summary23))] <- 4
# df$group[grepl("IMMUNE", rownames(common.geneset.summary23)) | grepl("HUMORAL", rownames(common.geneset.summary23))] <- 6
# df$group[grepl("VESICLE", rownames(common.geneset.summary23))] <- 7
# df <- df[order(df$group,decreasing=FALSE),]
# df$group <- factor(df$group,levels=1:7, labels=c("Mito Translation","Mito Transport","Mito Organization","Amino Acid","Translation","Immune","Vesicle"))
# common.geneset.summary23.rank <- common.geneset.summary23[as.character(df$geneset),]
# an <- HeatmapAnnotation(genese=df$group,which="row")
# Heatmap(as.matrix(common.geneset.summary23.rank),name="NES",
#         cluster_rows=FALSE,cluster_columns=FALSE,
#         row_names_side="left",
#         column_title="GLM-GSEA-Proteome",
#         rect_gp=gpar(col = "gray40",lwd=1),
#         column_names_gp=gpar(fontsize = 9),
#         row_names_gp=gpar(fontsize = 9),
#         left_annotation=an,
#         row_names_max_width = unit(10, "cm")
# )

# df <- data.frame(geneset=rownames(common.geneset.summary23),
#                  group=ifelse(grepl("MITOC", rownames(common.geneset.summary23)),1,6))
# df$group[df$group==1 & grepl("TRANSPORT", rownames(common.geneset.summary23))] <- 2
# df$group[grepl("AMINO_ACID", rownames(common.geneset.summary23))] <- 3
# df$group[grepl("DNA", rownames(common.geneset.summary23)) & grepl("REPLICATION", rownames(common.geneset.summary23))] <- 4
# df$group[grepl("DNA_REPAIR", rownames(common.geneset.summary23))] <- 5
# df$group[grepl("ACTIN", rownames(common.geneset.summary23)) | grepl("MYOFIBRIL", rownames(common.geneset.summary23)) |
#            grepl("CYTOSKELETON", rownames(common.geneset.summary23)) | grepl("MUSCLE", rownames(common.geneset.summary23)) ] <- 7
# df <- df[order(df$group,decreasing=FALSE),]
# df$group <- factor(df$group,levels=1:7, labels=c("Mito Translation","Mito Transport","Amino Acid","DNA Replication","DNA Repair","Translation","CytoSkeleton"))
# common.geneset.summary23.rank <- common.geneset.summary23[as.character(df$geneset),]
# an <- HeatmapAnnotation(geneset=df$group,which="row")
# Heatmap(as.matrix(common.geneset.summary23.rank),name="NES",
#         cluster_rows=FALSE,cluster_columns=FALSE,
#         row_names_side="left",
#         column_title="GLM-GSEA-Proteome",
#         rect_gp=gpar(col = "gray40",lwd=1),
#         column_names_gp=gpar(fontsize = 9),
#         row_names_gp=gpar(fontsize = 9),
#         left_annotation=an,
#         row_names_max_width = unit(10, "cm")
# )
