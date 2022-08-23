library(preprocessCore)
setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/211011.heatmap.fig3/")
data1<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/211011.heatmap.fig3/comprehensive.table.211011.txt",
                  sep="\t",header = T)
rownames(data1)<-data1$item
data2<-data1[,c("Pan.DNA.RNA.Corr","CCLE.DNA.RNA.Corr","NCI60.DNA.RNA.Corr","hCEC.all.DNA.RNA.Corr",
                "Pan.RNA.Pro.Corr","Normal.RNA.Pro.Corr","GTEx.RNA.Pro.Corr",
                "CCLE.RNA.Pro.Corr" ,                 
               "NCI60.RNA.Pro.Corr",               
               "hCEC.all.RNA.Pro.Corr"  )]

colnames(data2)[1]<-"CPTAC.DNA.RNA.Corr"
colnames(data2)[5]<-"CPTAC.RNA.Pro.Corr"
#colnames(data2)[6]<-"Wang, et al.Normal.RNA.Pro.Corr"
#data3<-normalize.quantiles(as.matrix(data2),copy=TRUE)
data3<-as.data.frame(scale(data2))
data3$median.dna.rna<-apply(data3[,1:4],1,mean)
data3$median.rna.pro<-apply(data3[,5:10],1,mean)
data3$delta<-data3$median.dna.rna-data3$median.rna.pro
#rownames(data3)<-rownames(data2)
#colnames(data3)<-colnames(data2)
#data3<-scale(data2)
data3[data3>1]<-1
data3[data3<(-1)]<-(-1)
#data4<-data3[order( -data3[,5]),]
data3<-data3[order( data3[,13]),]
data4<-data3[,-c(11:13)]

library(circlize)
library(ComplexHeatmap)
library(ggcorrplot)
library(colorspace)
ht_opt(RESET = TRUE)
ht_opt(heatmap_column_names_gp = gpar(fontface = "italic",fontsize=8), 
       heatmap_row_names_gp= gpar(fontsize = 8),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE
)    

hcl_palettes(plot = T)
col_fun = diverging_hcl("Blue-Red2",n=5)

mat<-as.matrix(data4)
ha<-Heatmap(mat,name = "Z-score", col = col_fun,
            cluster_columns = FALSE, show_row_dend = F, rect_gp = gpar(col= "white"), 
            show_column_names = T,cluster_rows = F,)

mat_cor<-as.matrix(data2[,-11:-13])
mat1<-round(cor(mat_cor), 1)
col_fun1 = diverging_hcl("Green-Orange",n=5)
ha1<-Heatmap(mat1,name = "Corr", col = col_fun1,
            cluster_columns = FALSE, show_row_dend = FALSE, rect_gp = gpar(col= "white"), 
            show_column_names = F,cluster_rows = F,clustering_distance_rows = "pearson",
            clustering_method_rows = "complete")

ht_list = ha1 %v% ha
draw(ht_list)

pdf(paste0("combine.TD_dna_rna_211014.pdf"),width = 6,height=10)
draw(ht_list)
dev.off()











