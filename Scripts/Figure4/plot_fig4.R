### this script is use to plot figure 4

library(circlize)
library(ComplexHeatmap)

rm(list=ls())

# ### check the ribosome genes in CPTAC and CCLE
# ribosome_cptac_path <- "/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis/CPTAC_pan_commonGenes_lm_Protein_AS_statistic.txt"
# ribosome_ccle_path <- "/Users/pc2644/Documents/DM_Aneuploidy/Compensation/CCLE/CCLE_pan_commonGenes_lm_Protein_AS_statistic.txt"
# 
# ribosome_cptac <- read.table(ribosome_cptac_path, sep="\t", stringsAsFactors=F, header=T)
# colnames(ribosome_cptac) <- paste0(colnames(ribosome_cptac),"_CPTAC")
# ribosome_ccle <- read.table(ribosome_ccle_path, sep="\t", stringsAsFactors=F, header=T)
# colnames(ribosome_ccle) <- paste0(colnames(ribosome_ccle),"_CCLE")
# 
# table <- merge(ribosome_cptac,ribosome_ccle,by.x="gene_CPTAC",by.y="gene_CCLE",sort=F)
# table_filter <- table[,c(1,4,6,9,11)]
# colnames(table_filter)[1] <- "gene"
# rownames(table_filter) <- table_filter$gene
# table_filter <- table_filter[,-1]
# 
# ### read ribosome genes
# ribo_path <- "/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis/ribosome_genelist.txt"
# ribo <- read.table(ribo_path,header=T,sep="\t",stringsAsFactors=F)
# ribo_common <- intersect(ribo$gene,rownames(table_filter))
# table_ribo <- table_filter[ribo_common,]
# plot(table_ribo$AS_Tvalue_CPTAC,table_ribo$AS_Tvalue_CCLE)

### fig. 4b heatmap, genesets related to gene expression, different model
geneset <- read.table("/Users/pc2644/Dropbox (NYU Langone Health)/Davoli LAB/Manuscripts/Proteogenomic Analysis of Aneuploidy/Manuscript_2021/figure_2021/fig4/GSEA/selected_geneset_pan_cancer_CPTAC_CCLE.txt",
                      sep="\t", stringsAsFactors=F, header=T)

rownames(geneset) <- geneset$geneset
geneset$process <- factor(geneset$process, levels=unique(geneset$process))
geneset[geneset=="NS"] <- NA
geneset[,-1:-3] <- apply(geneset[,-1:-3],2,as.numeric)

row_an<-HeatmapAnnotation(pathway=geneset$process,
                          show_legend=T,
                          show_annotation_name=c(pathway = FALSE),
                          annotation_label=NULL,
                          which="row")

column_an<-HeatmapAnnotation(sample=c(rep("Tissue",7),rep("Cell",5)),
                             which="column",
                             show_annotation_name=c(sample = FALSE),
                             col=list(sample=c("Cell"="#CC6677","Tissue"="#88CCEE")))

mycol <- colorRamp2(c(-3.5,0,3.5), c("blue","white","red"))

# pdf("/Users/pc2644/Desktop/selected_geneset_pan_cancer_CPTAC_CCLE_Heatmap.pdf", width=8.5, height=5)
Heatmap(as.matrix(geneset[,-1:-3]),col=mycol,name="GSEA NES",
        cluster_rows=FALSE,cluster_columns=FALSE,
        row_names_side="left",top_annotation=column_an,
        column_names_side="bottom",right_annotation=row_an,
        column_split=factor(c(rep("Tissue",7),rep("Cell",5)),levels=c("Tissue","Cell")),
        row_split=geneset$process,
        rect_gp=gpar(col = "gray40",lwd=0.3),
        row_names_rot=0,
        row_names_gp=gpar(fontsize = 8),
        column_names_gp=gpar(fontsize = 8),
        row_names_max_width=max_text_width(rownames(geneset), gp = gpar(fontsize = 12))
)
# dev.off()

### fig. 4d heatmap, genesets related to gene expression, different cancers
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/GeneExpressionPathway/CPTAC_CCLE_GSEA_NES_threshold_LM_Protein_AS.RData")
if (sum(rownames(geneset)!=rownames(GSEA_filter_NES_threshold))!=0) {stop("check the ordering of genesets!")}
geneset2 <- cbind(geneset[,1:3],GSEA_filter_NES_threshold)
geneset2[geneset2=="NS"] <- NA
geneset2[,-1:-3] <- apply(geneset2[,-1:-3],2,as.numeric)

column_an2<-HeatmapAnnotation(sample=c(rep("Tissue",3),rep("Cell",3)),
                              which="column",
                              show_annotation_name=c(sample = FALSE),
                              col=list(sample=c("Cell"="#CC6677","Tissue"="#88CCEE")))
mycol2 <- colorRamp2(c(-3,0,3), c("blue","white","red"))

# pdf("/Users/pc2644/Desktop/selected_geneset_specific_cancer_CPTAC_CCLE_Heatmap.pdf", width=8, height=5)
Heatmap(as.matrix(geneset2[,-1:-3]),col=mycol2,name="GSEA NES",
        cluster_rows=FALSE,cluster_columns=FALSE,
        row_names_side="left",top_annotation=column_an2,
        column_names_side="bottom",right_annotation=row_an,
        column_split=factor(c(rep("Tissue",3),rep("Cell",3)),levels=c("Tissue","Cell")),
        row_split=geneset$process,
        rect_gp=gpar(col = "gray40",lwd=0.3),
        row_names_rot=0,
        row_names_gp=gpar(fontsize = 8),
        column_names_gp=gpar(fontsize = 8),
        row_names_max_width=max_text_width(rownames(geneset), gp = gpar(fontsize = 12))
)
# dev.off()

### fig. S4 heatmap, genesets related to gene expression, different model (RNA)
load("/Users/pc2644/Documents/DM_Aneuploidy/Compensation/PanAnalysis_CPTAC/CPTAC_CCLE_GSEA_NES_threshold_LM_RNA_AS.RData")
if (sum(rownames(geneset)!=rownames(GSEA_filter1_NES_threshold_pick))!=0) {stop("check the ordering of genesets!")}
geneset3 <- cbind(geneset[,1:3],GSEA_filter1_NES_threshold_pick)
geneset3[geneset3=="NS"] <- NA
geneset3[,-1:-3] <- apply(geneset3[,-1:-3],2,as.numeric)

mycol3 <- colorRamp2(c(-2.5,0,2.5), c("blue","white","red"))

pdf("/Users/pc2644/Desktop/selected_geneset_specific_cancer_CPTAC_CCLE_RNA_Heatmap.pdf", width=8, height=5)
Heatmap(as.matrix(geneset3[,-1:-3]),col=mycol2,name="GSEA NES",
        cluster_rows=FALSE,cluster_columns=FALSE,
        row_names_side="left",top_annotation=column_an,
        column_names_side="bottom",right_annotation=row_an,
        column_split=factor(c(rep("Tissue",7),rep("Cell",5)),levels=c("Tissue","Cell")),
        row_split=geneset$process,
        rect_gp=gpar(col = "gray40",lwd=0.3),
        row_names_rot=0,
        row_names_gp=gpar(fontsize = 8),
        column_names_gp=gpar(fontsize = 8),
        row_names_max_width=max_text_width(rownames(geneset), gp = gpar(fontsize = 12))
)
dev.off()
