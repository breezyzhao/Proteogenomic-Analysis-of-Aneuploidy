cancer<-c("Colon")
setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210910.complex/loss.test/")

data_cnv<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/Colon-Prospective-Data/raw_from_python/colon_CNV.tsv",
                     sep="\t",header = T)
data_rna<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/Colon-Prospective-Data/raw_from_python/colon_transcriptomics.tsv",
                     sep="\t",header = T)
data_pro<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/Colon-Prospective-Data/raw_from_python/colon_proteomics.tsv",
                     sep="\t",header = T)
mutation<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/Colon/COADREAD_210528_mutation.txt",
                     sep="\t",header = T)


library(ggplot2)
library(boot)
if(cancer=="HNSC"| cancer=="OV"|cancer=="Colon"|cancer=="Endometrial"|cancer=="LUAD"){
  rownames(data_cnv)<-data_cnv$Patient_ID #BREAST,CCRCC:Name; HNSC,LUAD,Endometrial,OV,Colon:Patient_ID
}else if(cancer=="BREAST"| cancer=="ccRCC"){
  rownames(data_cnv)<-data_cnv$Name
}
data_cnv.use<-as.data.frame(t(data_cnv[,-1]))
cnv_by_gene1<-data_cnv.use[order(row.names(data_cnv.use)) , order(colnames(data_cnv.use))]
cnv_by_gene1<-cnv_by_gene1[,-3:-4]
cnv_by_gene<-as.data.frame(apply(cnv_by_gene1,2,as.numeric))
rownames(cnv_by_gene)<-rownames(cnv_by_gene1)


rownames(data_rna)<-data_rna$Patient_ID
data_rna.use<-as.data.frame(t(scale(data_rna[,-1]))) ##normalize based on sample
rna_by_gene1<-data_rna.use[order(row.names(data_rna.use)) , order(colnames(data_rna.use))]
median<-as.data.frame(apply(rna_by_gene1,1,function (x) median(x,na.rm=T)))
median$group<-ifelse(median$`apply(rna_by_gene1, 1, function(x) median(x, na.rm = T))`>=quantile(median$`apply(rna_by_gene1, 1, function(x) median(x, na.rm = T))`,0.0),1,0)
rna_by_gene<-rna_by_gene1[rownames(rna_by_gene1)%in% rownames(median[median$group==1,]),]

if(cancer=="HNSC"| cancer=="Endometrial"|cancer=="Colon"){
  rownames(data_pro)<-data_pro$Patient_ID#BREAST,OV,CCRCC,LUAD:Name hnsc,Endometrial,Colon:Patient_ID
}else if(cancer=="BREAST"| cancer=="OV"|cancer=="ccRCC"|cancer=="LUAD"){
  rownames(data_pro)<-data_pro$Name
}
data_pro.use<-as.data.frame(t(data_pro[,-1]))
pro_by_gene1<-data_pro.use[order(row.names(data_pro.use)) , order(colnames(data_pro.use))]
pro_by_gene1<-pro_by_gene1[,-3:-4]
pro_by_gene<-as.data.frame(apply(pro_by_gene1,2,as.numeric))
rownames(pro_by_gene)<-rownames(pro_by_gene1)

if(cancer=="Colon"){
  mutation.use<-mutation[mutation$Total_Num_Mutations<600,]
}else{
  mutation.use<-mutation[mutation$Total_Num_Mutations>=0,]
}

########### deal with data
datanames1<-intersect(rownames(cnv_by_gene),rownames(rna_by_gene))
datanames1.1<-intersect(row.names(pro_by_gene),datanames1)
datanames2<-intersect(colnames(cnv_by_gene),colnames(rna_by_gene))
datanames2.1<-intersect(colnames(pro_by_gene),datanames2)
datanames2.2<-intersect(mutation.use$Sample_Barcode,datanames2.1)
rna_by_gene_small<-rna_by_gene[rownames(rna_by_gene) %in% datanames1.1, colnames(rna_by_gene) %in% datanames2.2]
rna_by_gene_small1<-as.data.frame(scale(rna_by_gene_small))
cnv_by_gene_small<-cnv_by_gene[rownames(cnv_by_gene) %in% datanames1.1, colnames(cnv_by_gene) %in% datanames2.2]
pro_by_gene_small<-pro_by_gene[rownames(pro_by_gene) %in% datanames1.1,colnames(pro_by_gene) %in% datanames2.2]

gene.mean<-as.data.frame(cbind(rownames(cnv_by_gene_small),apply(cnv_by_gene_small,1,mean)))
gene.mean$V2<-as.numeric(gene.mean$V2)

corr<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/Colon.corr.table.txt",
                 sep="\t",header = T)

corum.gene<-corr[corr$corum.x==TRUE,]$gene

data.temp<-merge(gene.mean,corr,by.x="V1",by.y="gene") ### genes we need to find high cor.rna.dna low cor.rna.pro

corum.loss<-data.temp[data.temp$V2<(-0.2) & data.temp$V1%in% corum.gene,]$V1
corum.gain<-data.temp[data.temp$V2>(0.2) & data.temp$V1%in% corum.gene,]$V1

for (k in 1:length(corum.loss)){
#gene1<-c("FAM210B",corum.gain[k]) #FAM210B nocorum PRPF6 Courm
gene2<-c("CTDNEP1",corum.loss[k]) #INPPRK corum CTDNEP1 nocorum

#gene.list<-list(gene1)
gene.list<-list(gene2)
for (i in 1:length(gene.list)){

  gene<-gene.list[[i]]
cnv_by_gene.use<-as.data.frame(t(cnv_by_gene_small[rownames(cnv_by_gene_small) %in% gene,]))
rna_by_gene.use<-as.data.frame(t(rna_by_gene_small[rownames(rna_by_gene_small) %in% gene,]))
pro_by_gene.use<-as.data.frame(t(pro_by_gene_small[rownames(pro_by_gene_small) %in% gene,]))

cnv<-melt(cnv_by_gene.use)
rna<-melt(rna_by_gene.use)
pro<-melt(pro_by_gene.use)

merge.data<-cbind(cnv,rna,pro)
colnames(merge.data)[2]<-"cnv"
colnames(merge.data)[4]<-"rna"
colnames(merge.data)[6]<-"pro"
plot.data<-merge.data[,c(1,2,4,6)]

plot.data$variable<-factor(plot.data$variable,levels = c("CTDNEP1",corum.loss[k]))

col1<-"#aaa900" # buddha gold#"#FF0033" ###red nocorum
col2<-"#6500aa" # purple #"#3300FF" ## blue corum

g1<-ggplot(plot.data, aes(x=cnv, y=rna, color=variable)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  ggtitle("DNA-RNA Corr")+
  theme_classic()+
  annotate(x = min(plot.data[plot.data$variable%in%gene[1],]$cnv)+0.4, y = max(plot.data[plot.data$variable%in%gene[1],]$rna), 
           label=paste("r = ", round(cor(plot.data[plot.data$variable%in%gene[1],]$cnv, plot.data[plot.data$variable%in%gene[1],]$rna),3)), 
           geom="text", size=5,col=col1,hjust=0)+
  annotate(x=min(plot.data[plot.data$variable%in%gene[1],]$cnv)+0.4, y=max(plot.data[plot.data$variable%in%gene[1],]$rna)-0.6,
           label=paste("r = ", round(cor(plot.data[plot.data$variable%in%gene[2],]$cnv, plot.data[plot.data$variable%in%gene[2],]$rna),3)), 
           geom="text", size=5,col=col2,hjust=0)+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))

g2<-ggplot(plot.data, aes(x=pro, y=rna, color=variable)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  ggtitle("RNA-Pro Corr")+
  theme_classic()+
  annotate(x = min(plot.data[plot.data$variable%in%gene[1],]$pro,na.rm = T)+0.4, y = max(plot.data[plot.data$variable%in%gene[1],]$rna,na.rm = T), 
           label=paste("r = ", round(cor(plot.data[plot.data$variable%in%gene[1],]$pro, plot.data[plot.data$variable%in%gene[1],]$rna,use="na.or.complete"),3)), 
           geom="text", size=5,col=col1,hjust=0)+
  annotate(x=min(plot.data[plot.data$variable%in%gene[1],]$pro,na.rm = T)+0.4, y=max(plot.data[plot.data$variable%in%gene[1],]$rna,na.rm = T)-0.6,
           label=paste("r = ", round(cor(plot.data[plot.data$variable%in%gene[2],]$pro, plot.data[plot.data$variable%in%gene[2],]$rna,use="na.or.complete"),3)), 
           geom="text", size=5,col=col2,hjust=0)+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))

g3<-ggplot(plot.data, aes(x=cnv, y=pro, color=variable)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  ggtitle("DNA-Pro Corr")+
  theme_classic()+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))

figure<-ggarrange(g1,g2,g3, ncol=3,nrow=1,common.legend = TRUE,legend = "right")
library(gridExtra)
library(ggpubr)
pdf(paste0(corum.loss[k],"_DNA_RNA_PRO-211029.pdf"),onefile=FALSE,height=3,width=10)
print(figure)
dev.off() 
}
}

