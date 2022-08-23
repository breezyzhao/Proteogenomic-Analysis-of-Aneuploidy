predictvals <- function(model, xvar, yvar, xrange = NULL, samples = 100, ...) {
  
  # If xrange isn't passed in, determine xrange from the models.
  # Different ways of extracting the x range, depending on model type
  if (is.null(xrange)) {
    if (any(class(model) %in% c("lm", "glm")))
      xrange <- range(model$model[[xvar]])
    else if (any(class(model) %in% "loess"))
      xrange <- range(model$x)
  }
  
  newdata <- data.frame(x = seq(xrange[1], xrange[2], length.out = samples))
  names(newdata) <- xvar
  newdata[[yvar]] <- predict(model, newdata = newdata, ...)
  newdata
}


setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/211026.control.test.fig3a.b/")
library(ggplot2)
library(biomaRt)
library(ggpointdensity)
library(viridis)
library(ggrepel)
library(ggpubr)
library(data.table)
library(plyr)
library(dplyr)
library(ggpmisc)
corum <- as.data.frame(read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/database/coreComplexesv3.0.txt", sep = "\t"))
corum<-corum[corum$Organism=="Human",]

master_list <- strsplit(as.character(corum$subunits.Entrez.IDs.), split = ";")
master_list <- unique(as.numeric(as.character(unlist(master_list))))
# Create master list of CORUM UniProt IDs
master_list_uniprot <- strsplit(as.character(corum$subunits.UniProt.IDs.), split = ";")
master_list_uniprot <- unique(as.character(unlist(master_list_uniprot)))
# Create master list of CORUM gene names (use this)
master_list_names <- strsplit(as.character(corum$subunits.Gene.name.), split = ";")
master_list_names <- unique(as.character(unlist(master_list_names)))


colon<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct/Colon.corr.table.txt",sep="\t")
colnames(colon)[2]<-"colon.cor.rna.dna"
colnames(colon)[5]<-"colon.cor.rna.pro"
colnames(colon)[7]<-"colon.cor.dna.pro"
breast<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct/BREAST.corr.table.txt",sep="\t")
colnames(breast)[2]<-"breast.cor.rna.dna"
colnames(breast)[5]<-"breast.cor.rna.pro"
colnames(breast)[7]<-"breast.cor.dna.pro"
ovarian<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct/OV.corr.table.txt",sep="\t")
colnames(ovarian)[2]<-"ovarian.cor.rna.dna"
colnames(ovarian)[5]<-"ovarian.cor.rna.pro"
colnames(ovarian)[7]<-"ovarian.cor.dna.pro"
ccRCC<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct/ccRCC.corr.table.txt",sep = "\t")
colnames(ccRCC)[2]<-"ccRCC.cor.rna.dna"
colnames(ccRCC)[5]<-"ccRCC.cor.rna.pro"
colnames(ccRCC)[7]<-"ccRCC.cor.dna.pro"
endometrial<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct/Endometrial.corr.table.txt",sep = "\t")
colnames(endometrial)[2]<-"endometrial.cor.rna.dna"
colnames(endometrial)[5]<-"endometrial.cor.rna.pro"
colnames(endometrial)[7]<-"endometrial.cor.dna.pro"
luad<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct/LUAD.corr.table.txt",sep = "\t")
colnames(luad)[2]<-"luad.cor.rna.dna"
colnames(luad)[5]<-"luad.cor.rna.pro"
colnames(luad)[7]<-"luad.cor.dna.pro"
hnsc<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct/HNSC.corr.table.txt",sep = "\t")
colnames(hnsc)[2]<-"hnsc.cor.rna.dna"
colnames(hnsc)[5]<-"hnsc.cor.rna.pro"
colnames(hnsc)[7]<-"hnsc.cor.dna.pro"
score<-na.omit(read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/36.evolution/hg19.gene.score.txt",sep="\t",header=T))
score$corum<-ifelse(score$genes%in%master_list_names,TRUE,FALSE)

Pan.data<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct/Pan-Cancer.corr.table.txt",
                     sep="\t",header = T)
colnames(Pan.data)[2]<-"pan.cor.rna.dna"
colnames(Pan.data)[4]<-"pan.cor.rna.pro"
colnames(Pan.data)[5]<-"pan.cor.dna.pro"

colnames(Pan.data)[1]<-"gene"
cof<-list()
names<-list(colon,breast,ovarian,ccRCC,endometrial,luad,hnsc,Pan.data)
names1<-c("Colon","Breast","OV","ccRCC","Endometrial","LUAD","HNSC","PAN-Cancer")
for (k in 1:8){
data1<-names[[k]]
if(names1[k]=="PAN-Cancer"){
  colnames(data1)[2]<-"corr.rna.dna"
  colnames(data1)[4]<-"corr.rna.pro"
}else{
  colnames(data1)[2]<-"corr.rna.dna"
  colnames(data1)[5]<-"corr.rna.pro"
}
rna.dna<-quantile(data1$corr.rna.dna, probs = c(0.35, 0.65))
rna.dna.low<-rna.dna[1]
rna.dna.high<-rna.dna[2]
rna.prot<-quantile(data1$corr.rna.pro, probs = c(0.35, 0.65))
rna.prot.low<-rna.prot[1]
rna.prot.high<-rna.prot[2]
cof1<-as.data.frame(cbind(rna.dna.high,rna.dna.low,rna.prot.high,rna.prot.low))
cof1$cancer<-names1[k]
cof[[k]]<-cof1
for (i in 1:nrow(data1)){
  #i<-56
if (data1$corr.rna.dna[i]>rna.dna.high & data1$corr.rna.pro[i]>rna.prot.high){
  data1$group[i]<-c("Corr_RNA_DNA_High & Corr_RNA_Pro_High")
}else if(data1$corr.rna.dna[i]>rna.dna.high & data1$corr.rna.pro[i]<rna.prot.low){
  data1$group[i]<-c("Corr_RNA_DNA_High & Corr_RNA_Pro_Low")
}else if(data1$corr.rna.dna[i]<rna.dna.low & data1$corr.rna.pro[i]<rna.prot.low){
  data1$group[i]<-c("Corr_RNA_DNA_Low & Corr_RNA_Pro_Low")
}else if (data1$corr.rna.dna[i]<rna.dna.low & data1$corr.rna.pro[i]>rna.prot.high){
  data1$group[i]<-c("Corr_RNA_DNA_Low & Corr_RNA_Pro_High")
}else{
  data1$group[i]<-c("Others")
}
}

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   version = 75)
genes <- data1$gene
G_list <- getBM(filters = "hgnc_symbol", 
                attributes = c("hgnc_symbol", "uniprot_swissprot","uniprot_swissprot_accession"
                ),
                values = genes, mart = mart)
gene.test1<-merge(G_list,data1,by.x="hgnc_symbol",by.y="gene")
gene.test2<- na.omit(gene.test1 %>% dplyr::mutate_all(na_if,""))

#write.table(gene.test2,paste0(names1[k],"-cat-0.35.txt"),sep="\t",row.names = F,quote = F)


data1$group<-factor(data1$group,level=c("Corr_RNA_DNA_High & Corr_RNA_Pro_High",
                                        "Corr_RNA_DNA_High & Corr_RNA_Pro_Low",
                                        "Corr_RNA_DNA_Low & Corr_RNA_Pro_Low",
                                        "Corr_RNA_DNA_Low & Corr_RNA_Pro_High"))


counts <- ddply(data1, .(round_any(data1$corr.rna.dna,0.05), round_any(data1$corr.rna.pro,0.05)), nrow)
names(counts) <- c("corr.rna.dna", "corr.rna.pro", "Freq")
testlist<-list()
#plot(counts$corr.rna.dna,counts$corr.rna.pro)
len.dna.rna<-unique(counts$corr.rna.dna)
for (m in 1:length(len.dna.rna)){
  corr.rna.dna<-len.dna.rna[m]
  count.temp<-counts[counts$corr.rna.dna%in%corr.rna.dna,]
  max.temp<-max(count.temp$Freq)
  rna.pro.temp1<-count.temp[count.temp$Freq %in%max.temp,]
  if(nrow(rna.pro.temp1)==1){
  corr.rna.pro<-rna.pro.temp1$corr.rna.pro[1]
  }else{
    corr.rna.pro<-mean(rna.pro.temp1$corr.rna.pro)
  }
  freq<-max.temp
  data.temp<-as.data.frame(cbind(corr.rna.dna,corr.rna.pro,freq))
  testlist[[m]]<-data.temp
}
test1<-do.call(rbind,testlist)
test<-test1[test1$corr.rna.dna>-0.5,]
#plot(test$dna.rna,test$rna.pro)
#cor.test(test$dna.rna,test$rna.pro)
model <- lm(corr.rna.pro ~ corr.rna.dna, test)
cor1<-cor.test(test$corr.rna.dna,test$corr.rna.pro,method = "pearson")
cor2<-cor.test(test$corr.rna.dna,test$corr.rna.pro,method = "spearman")
model.use <- predictvals(model, "corr.rna.dna","corr.rna.pro",xrange = c(min(model$model$corr.rna.dna),max(model$model$corr.rna.dna)))

g1<-ggplot(data1, aes(x = corr.rna.dna, y = corr.rna.pro)) +
  geom_pointdensity(adjust = 0.1,size=0.01) +
  geom_density_2d_filled(alpha = 0.8) +
  #geom_density_2d(size = 0.01, colour = "black")+
  scale_color_viridis()+
  xlab("DNA-RNA Corr")+
  ylab("RNA-Pro Corr")+
  xlim(-0.5,1)+ylim(-0.6,1)+
  theme_classic()+
  theme(axis.text=element_text(size=12,face = "bold"),
         axis.title=element_text(size=14,face="bold"),
         plot.title = element_text(size = 14, face = "bold"),
         legend.text=element_text(size=10,face = "bold"),
         legend.title=element_text(size=12,face = "bold"))+
  theme(legend.position="right")+
  geom_line(data = model.use)+
  annotate("text", x = -0.4, y = -0.4, label = paste0("r=",round(cor1$estimate,3),",p=",signif(cor1$p.value,3)),
           hjust=0,vjust=0.5)+
  annotate("text", x = -0.4, y = -0.3, label = paste0("rho=",round(cor2$estimate,3),",p=",signif(cor2$p.value,3)),
           hjust=0,vjust=0.5)
  
  #geom_smooth(method='loess', formula= y~x,colour="black",se = FALSE)+
  #stat_cor(method = "pearson", label.x.npc = "left", label.y.npc ="bottom" ,size = 5)

g2<-ggplot(test, aes(x = corr.rna.dna, y = corr.rna.pro)) +
  geom_pointdensity(adjust = 0.1,size=1) +
  #geom_density_2d(size = 0.01, colour = "black")+
  scale_color_viridis()+
  xlab("DNA-RNA Corr")+
  ylab("RNA-Pro Corr")+
  xlim(-0.5,1)+ylim(-0.6,1)+
  theme_classic()+
  theme(axis.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=10,face = "bold"),
        legend.title=element_text(size=12,face = "bold"))+
  theme(legend.position="right")+
  geom_line(data = model.use)

g3<-ggscatter(test, x = "corr.rna.dna", y="corr.rna.pro", add = "reg.line") +
  stat_cor(label.x = 0.2, label.y = 1) +
  stat_regline_equation(label.x = -0.5, label.y = 1)+
  xlim(-0.5,1)+ylim(-0.6,1)



  
pdf(paste0(names1[k],".test-0.35.pdf"),width = 10,height = 8)
print(g1)
dev.off()
pdf(paste0(names1[k],".test-control.pdf"),width = 10,height = 8)
print(g3)
dev.off()

}
 #cof.data<-do.call(rbind,cof)
 #write.table(cof.data,"cof.cutoff.txt",sep="\t",row.names = F,quote = F)

# library(devtools)
# install_github("ProcessMiner/nlcor")
# library(nlcor)
# #data1
# x<-data1$corr.rna.dna
# y<-data1$corr.rna.pro
# cor(x,y)
# # [1] 6.488616e-17
# # nonlinear correlation is more representative
# nlcor(x,y, plt = T)
