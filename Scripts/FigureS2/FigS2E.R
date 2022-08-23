setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210917.degradation/")
library(tidyr)
library(ggridges)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(data.table)
corum <- as.data.frame(read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/database/coreComplexesv3.0.txt", sep = "\t"))
corum<-corum[corum$Organism=="Human",]
master_list <- strsplit(as.character(corum$subunits.Entrez.IDs.), split = ";")
master_list <- unique(as.numeric(as.character(unlist(master_list))))
# Create master list of CORUM UniProt IDs
master_list_uniprot <- strsplit(as.character(corum$subunits.UniProt.IDs.), split = ";")
master_list_uniprot <- unique(as.character(unlist(master_list_uniprot)))
# Create master list of CORUM gene names (use this)
master_list_names <- strsplit(as.character(corum$subunits.Gene.name.), split = ";")
master_list_names <- unique(as.character(unlist(master_list_names)))

degradation<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/dataset_outside/degradation/RPE-1.txt",sep="\t",header = T)
deg<-degradation[,c("Gene.names","Degradation.profile")]
deg1<-separate_rows(deg,Gene.names)
para<-c("ED","NED")
deg2<-deg1[deg1$Degradation.profile %in% para,]
score<-deg2
score$Corum<-ifelse(score$Gene.names%in%master_list_names,"Corum","NoCorum")
colnames(score)[1]<-"genes"

colon<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/Colon.corr.table.txt",sep="\t")
colnames(colon)[2]<-"colon.cor.rna.dna"
colnames(colon)[5]<-"colon.cor.rna.pro"
colnames(colon)[7]<-"colon.cor.dna.pro"
colon<-colon[colon$zero_percent.x<0.7,]
breast<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/BREAST.corr.table.txt",sep="\t")
colnames(breast)[2]<-"breast.cor.rna.dna"
colnames(breast)[5]<-"breast.cor.rna.pro"
colnames(breast)[7]<-"breast.cor.dna.pro"
breast<-breast[breast$zero_percent.x<0.7,]
ovarian<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/OV.corr.table.txt",sep="\t")
colnames(ovarian)[2]<-"ovarian.cor.rna.dna"
colnames(ovarian)[5]<-"ovarian.cor.rna.pro"
colnames(ovarian)[7]<-"ovarian.cor.dna.pro"
ovarian<-ovarian[ovarian$zero_percent.x<0.7,]
ccRCC<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/ccRCC.corr.table.txt",sep = "\t")
colnames(ccRCC)[2]<-"ccRCC.cor.rna.dna"
colnames(ccRCC)[5]<-"ccRCC.cor.rna.pro"
colnames(ccRCC)[7]<-"ccRCC.cor.dna.pro"
ccRCC<-ccRCC[ccRCC$zero_percent.x<0.7,]
endometrial<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/Endometrial.corr.table.txt",sep = "\t")
colnames(endometrial)[2]<-"endometrial.cor.rna.dna"
colnames(endometrial)[5]<-"endometrial.cor.rna.pro"
colnames(endometrial)[7]<-"endometrial.cor.dna.pro"
endometrial<-endometrial[endometrial$zero_percent.x<0.7,]
luad<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/LUAD.corr.table.txt",sep = "\t")
colnames(luad)[2]<-"luad.cor.rna.dna"
colnames(luad)[5]<-"luad.cor.rna.pro"
colnames(luad)[7]<-"luad.cor.dna.pro"
luad<-luad[luad$zero_percent.x<0.7,]
hnsc<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/HNSC.corr.table.txt",sep = "\t")
colnames(hnsc)[2]<-"hnsc.cor.rna.dna"
colnames(hnsc)[5]<-"hnsc.cor.rna.pro"
colnames(hnsc)[7]<-"hnsc.cor.dna.pro"
hnsc<-hnsc[hnsc$zero_percent.x<0.7,]

merge.data<-merge(merge(merge(merge(merge(merge(colon,breast,by="gene",all = T),ovarian,by="gene",all = T),ccRCC,by="gene",all = T),
                              endometrial,by="gene",all = T),luad,by="gene",all = T),hnsc,by="gene",all = T)
Pan.dna.rna<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% "cor.rna.dna"])
Pan.rna.pro<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% "cor.rna.pro"])

Pan.dna.rna$mean.rna.dna<-apply(Pan.dna.rna[,2:8],1,function(x)mean(x,na.rm=T))
Pan.rna.pro$mean.rna.pro<-apply(Pan.rna.pro[,2:8],1,function(x)mean(x,na.rm=T))

Pan.data<-merge(Pan.dna.rna,Pan.rna.pro,by="merge.data$gene",all = T)
colnames(Pan.data)[1]<-"gene"

names<-list(colon,breast,ovarian,ccRCC,endometrial,luad,hnsc,Pan.data)
names1<-c("Colon","Breast","OV","ccRCC","Endometrial","LUAD","HNSC","PAN-Cancer")

col1<-"#AAA900" ###budddha gold
col2<-"#6500aa" ## purple
col8.nc<-"#777600" ### er
col7.nc<-"#c4c34d" ### pm
col8<-"#3d0066" ## er
col7<-"#a366cc" ###pm

for (j in 1:8){
  #j<-3
  cancer<-"Degradation"
  CPTAC<-names1[j]
  merge.data<-names[[j]]
  if(names1[j]=="PAN-Cancer"){
    colnames(merge.data)[1]<-"genes" 
    score.data<-merge(score,merge.data[,c("genes","mean.rna.dna","mean.rna.pro")],by="genes")
    colnames(score.data)[4]<-"corr.rna.dna"
    colnames(score.data)[5]<-"corr.rna.pro"
  }else{
    colnames(merge.data)[1]<-"genes" 
    score.data<-merge(score,merge.data[,c(1,2,5)],by="genes")
    colnames(score.data)[4]<-"corr.rna.dna"
    colnames(score.data)[5]<-"corr.rna.pro"
  }

  corum<-score.data[score.data$Corum=="Corum",]
  corum$group<-"Corum:all"
  corum.ed<-score.data[score.data$Degradation.profile=="ED" & score.data$Corum=="Corum",]
  corum.ed$group<-"Corum:ED"
  corum.ned<-score.data[score.data$Degradation.profile=="NED" & score.data$Corum=="Corum",]
  corum.ned$group<-"Corum:NED"
  nocorum<-score.data[score.data$Corum=="NoCorum",]
  nocorum$group<-"NoCorum:all"
  nocorum.ed<-score.data[score.data$Degradation.profile=="ED" & score.data$Corum=="NoCorum",]
  nocorum.ed$group<-"NoCorum:ED"
  nocorum.ned<-score.data[score.data$Degradation.profile=="NED" & score.data$Corum=="NoCorum",]
  nocorum.ned$group<-"NoCorum:NED"
  
 # plot.data<-rbind(corum,corum.ed,corum.ned,nocorum,nocorum.ed,nocorum.ned)
#  plot.data$group<-factor(plot.data$group,levels = c("Corum:NED","NoCorum:NED","Corum:ED","NoCorum:ED",
#                                                     "Corum:all", "NoCorum:all"))  
  
  plot.data<-rbind(corum.ed,corum.ned,nocorum.ed,nocorum.ned)
  plot.data$group<-factor(plot.data$group,levels = c("Corum:NED","Corum:ED","NoCorum:NED","NoCorum:ED"))  
  
  g1<-ggplot(plot.data, aes(x = corr.rna.dna, y = group)) +
    geom_density_ridges(scale = 1.5,quantile_lines=TRUE,
                        quantile_fun=function(x,...)median(x),
                        aes(fill = group)) + 
    scale_fill_manual(values=c(alpha(col7,0.5),alpha(col7.nc,0.5),alpha(col8,0.5),alpha(col8.nc,0.5)))+
    scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
    scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
    coord_cartesian(clip = "off") +
    xlim(0, 1)+
    ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
    theme_classic()+
    geom_text(data=plot.data %>% group_by(group) %>% 
                summarise(corr.rna.dna=median(corr.rna.dna)),
              aes(label=round(corr.rna.dna,3)), 
              position=position_nudge(y=0.5,x=-0.1), colour="black", size=3.5)+# to avoid clipping of the very top of the top ridgeline
    theme_ridges()+
    theme(legend.position = "none")
  
  g2<-ggplot(plot.data, aes(x = corr.rna.pro, y = group)) +
    geom_density_ridges(scale = 1.5,quantile_lines=TRUE,
                        quantile_fun=function(x,...)median(x),
                        aes(fill = group)) + 
    scale_fill_manual(values=c(alpha(col7,0.5),alpha(col7.nc,0.5),alpha(col8,0.5),alpha(col8.nc,0.5)))+
    scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
    scale_x_continuous(expand = c(0, 0)) + # for both axes to remove unneeded padding
    coord_cartesian(clip = "off") +
    xlim(0, 1)+
    ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
    theme_classic()+
    geom_text(data=plot.data %>% group_by(group) %>% 
                summarise(corr.rna.pro=median(corr.rna.pro)),
              aes(label=round(corr.rna.pro,3)), 
              position=position_nudge(y=0.5,x=-0.1), colour="black", size=3.5)+# to avoid clipping of the very top of the top ridgeline
    theme_ridges()+
    theme(legend.position = "none")
  
  
  figure<-ggarrange(g1,g2,nrow=1,ncol=2)
  pdf(paste0(names1[j],".degradation.pdf"),width = 12)
  print(figure)
  dev.off()
}
  
  
  
  
  
  
  













