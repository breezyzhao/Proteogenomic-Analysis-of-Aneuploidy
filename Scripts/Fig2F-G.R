bootstrap<-function(dataset,iterations,index,groups,para1,para2){
  s.size.g1 <- length(dataset[,index][dataset[,groups]%in%para1])
  s.size.g2 <- length(dataset[,index][dataset[,groups]%in%para2])
  pool <- dataset[,index]
  obs.diff.b1 <- median (dataset[,index][dataset[,groups]==para1]) - median (dataset[,index][dataset[,groups]==para2])
  iterations <- iterations
  sampl.dist.b1 <- NULL
  for (i in 1 : iterations) {
    resample <- sample (c(1:length (pool)), length(pool), replace = TRUE) 
    g1.perm = pool[resample][1 : s.size.g1]
    g2.perm = pool[resample][(s.size.g1+1) : length(pool)]
    sampl.dist.b1[i] = median (g1.perm) - median (g2.perm) 
  }
  p.boot1 <- (sum ( abs(sampl.dist.b1) >= abs(obs.diff.b1)) + 1)/ (iterations+1)
  return(p.boot1)
}
set.seed(1234)
setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210806.phylopscore/")
#################################################################################################
#
#
#         COLON,breast,ovarian,CCRCC,LUAD,ENDOMETRIAL
#         Lung squamous, embargo date is Dec 1 2021
#         Head and neck is Aug 3, 2021
#         GBM is March 1 2021
#
#
#################################################################################################

corum <- as.data.frame(read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/database/coreComplexesv3.0.txt", sep = "\t"))
corum<-corum[corum$Organism=="Human",]

score<-na.omit(read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/36.evolution/hg19.gene.score.txt",sep="\t",header=T))
library(ggplot2)
library(boot)
library(ggpubr)
library(rstatix)
master_list <- strsplit(as.character(corum$subunits.Entrez.IDs.), split = ";")
master_list <- unique(as.numeric(as.character(unlist(master_list))))
# Create master list of CORUM UniProt IDs
master_list_uniprot <- strsplit(as.character(corum$subunits.UniProt.IDs.), split = ";")
master_list_uniprot <- unique(as.character(unlist(master_list_uniprot)))
# Create master list of CORUM gene names (use this)
master_list_names <- strsplit(as.character(corum$subunits.Gene.name.), split = ";")
master_list_names <- unique(as.character(unlist(master_list_names)))

score$corum<-ifelse(score$genes%in%master_list_names,TRUE,FALSE)

corums_cors <- score[score$corum == TRUE, "median.PhyloP.score"]
noncorums_cors <- score[score$corum == FALSE, "median.PhyloP.score"]

########################################################################################################
## generate pan-cancer data
###############################3


colon<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig1/data/Colon.corr.table.txt",sep="\t")
colnames(colon)[2]<-"colon.cor.rna.dna"
colnames(colon)[5]<-"colon.cor.rna.pro"
breast<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig1/data/BREAST.corr.table.txt",sep="\t")
colnames(breast)[2]<-"breast.cor.rna.dna"
colnames(breast)[5]<-"breast.cor.rna.pro"
ovarian<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig1/data/OV.corr.table.txt",sep="\t")
colnames(ovarian)[2]<-"ovarian.cor.rna.dna"
colnames(ovarian)[5]<-"ovarian.cor.rna.pro"
ccRCC<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig1/data/ccRCC.corr.table.txt",sep = "\t")
colnames(ccRCC)[2]<-"ccRCC.cor.rna.dna"
colnames(ccRCC)[5]<-"ccRCC.cor.rna.pro"
endometrial<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig1/data/Endometrial.corr.table.txt",sep = "\t")
colnames(endometrial)[2]<-"endometrial.cor.rna.dna"
colnames(endometrial)[5]<-"endometrial.cor.rna.pro"
luad<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig1/data/LUAD.corr.table.txt",sep = "\t")
colnames(luad)[2]<-"luad.cor.rna.dna"
colnames(luad)[5]<-"luad.cor.rna.pro"
hnsc<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig1/data/HNSC.corr.table.txt",sep = "\t")
colnames(hnsc)[2]<-"hnsc.cor.rna.dna"
colnames(hnsc)[5]<-"hnsc.cor.rna.pro"
score<-na.omit(read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/36.evolution/hg19.gene.score.txt",sep="\t",header=T))
score$corum<-ifelse(score$genes%in%master_list_names,TRUE,FALSE)


merge.data<-merge(merge(merge(merge(merge(merge(colon,breast,by="gene",all = T),ovarian,by="gene",all = T),ccRCC,by="gene",all = T),
                 endometrial,by="gene",all = T),luad,by="gene",all = T),hnsc,by="gene",all = T)
Pan.dna.rna<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% "cor.rna.dna"])
Pan.rna.pro<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% "cor.rna.pro"])

Pan.dna.rna$max.rna.dna<-apply(Pan.dna.rna[,2:8],1,function(x)max(x,na.rm=T))
Pan.rna.pro$max.rna.pro<-apply(Pan.rna.pro[,2:8],1,function(x)max(x,na.rm=T))

Pan.data<-merge(Pan.dna.rna,Pan.rna.pro,by="merge.data$gene",all = T)
colnames(Pan.data)[1]<-"gene"

names<-list(colon,breast,ovarian,ccRCC,endometrial,luad,hnsc,Pan.data)
names1<-c("Colon","Breast","OV","ccRCC","Endometrial","LUAD","HNSC","PAN-Cancer")
col1<-"#13502f" ###Dark green 
col2<-"#79b695" ## light green  ##76b493
cancer<-"PhyloP.score"

for (i in 1:length(names)){
#i<-1
  
merge.data<-names[[i]]
CPTAC<-names1[i]
  if(names1[i]=="PAN-Cancer"){
data.dna.rna<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% ".rna.dna"]) 
 colnames(data.dna.rna)[1]<-"genes" 
score.dna.rna<-merge(score,data.dna.rna[,c("genes","max.rna.dna")],by="genes")
score.dna.rna1<-score.dna.rna
score.dna.rna$group<-ifelse(score.dna.rna$median.PhyloP.score>quantile(score.dna.rna$median.PhyloP.score,0.7),"TOP30%",
                            ifelse(score.dna.rna$median.PhyloP.score<quantile(score.dna.rna$median.PhyloP.score,0.3),"BOT30%",NA))

score.dna.rna$Corum<-ifelse(score.dna.rna$corum==T,"Corum","No_Corum")
corum.top30<-quantile(score.dna.rna[score.dna.rna$Corum=="Corum",]$median.PhyloP.score,0.7,na.rm=T)
corum.bot30<-quantile(score.dna.rna[score.dna.rna$Corum=="Corum",]$median.PhyloP.score,0.3,na.rm=T)
nocorum.top30<-quantile(score.dna.rna[score.dna.rna$Corum=="No_Corum",]$median.PhyloP.score,0.7,na.rm=T)
nocorum.bot30<-quantile(score.dna.rna[score.dna.rna$Corum=="No_Corum",]$median.PhyloP.score,0.3,na.rm=T)
score.dna.rna$group_1<-ifelse(score.dna.rna$Corum=="Corum" & score.dna.rna$median.PhyloP.score>corum.top30,"Corum.TOP30%",
                              ifelse(score.dna.rna$Corum=="Corum" & score.dna.rna$median.PhyloP.score<corum.bot30,"Corum.BOT30%",
                                    ifelse(score.dna.rna$Corum=="No_Corum" & score.dna.rna$median.PhyloP.score>nocorum.top30,"No_Corum.TOP30%",
                                           ifelse(score.dna.rna$Corum=="No_Corum" & score.dna.rna$median.PhyloP.score<corum.bot30,"No_Corum.BOT30%",NA))))

score.dna.rna<-na.omit(score.dna.rna)
score.dna.rna$group_2<-"Cor.DNA.RNA"
#### rna.pro
data.rna.pro<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% ".rna.pro"]) 
colnames(data.rna.pro)[1]<-"genes" 
score.rna.pro<-merge(score,data.rna.pro[,c("genes","max.rna.pro")],by="genes")
score.rna.pro1<-score.rna.pro
score.rna.pro$group<-ifelse(score.rna.pro$median.PhyloP.score>quantile(score.rna.pro$median.PhyloP.score,0.7),"TOP30%",
                            ifelse(score.rna.pro$median.PhyloP.score<quantile(score.rna.pro$median.PhyloP.score,0.3),"BOT30%",NA))

score.rna.pro$Corum<-ifelse(score.rna.pro$corum==T,"Corum","No_Corum")
corum.top30<-quantile(score.rna.pro[score.rna.pro$Corum=="Corum",]$median.PhyloP.score,0.7,na.rm=T)
corum.bot30<-quantile(score.rna.pro[score.rna.pro$Corum=="Corum",]$median.PhyloP.score,0.3,na.rm=T)
nocorum.top30<-quantile(score.rna.pro[score.rna.pro$Corum=="No_Corum",]$median.PhyloP.score,0.7,na.rm=T)
nocorum.bot30<-quantile(score.rna.pro[score.rna.pro$Corum=="No_Corum",]$median.PhyloP.score,0.3,na.rm=T)
score.rna.pro$group_1<-ifelse(score.rna.pro$Corum=="Corum" & score.rna.pro$median.PhyloP.score>corum.top30,"Corum.TOP30%",
                              ifelse(score.rna.pro$Corum=="Corum" & score.rna.pro$median.PhyloP.score<corum.bot30,"Corum.BOT30%",
                                     ifelse(score.rna.pro$Corum=="No_Corum" & score.rna.pro$median.PhyloP.score>nocorum.top30,"No_Corum.TOP30%",
                                            ifelse(score.rna.pro$Corum=="No_Corum" & score.rna.pro$median.PhyloP.score<corum.bot30,"No_Corum.BOT30%",NA))))

score.rna.pro<-na.omit(score.rna.pro)
score.rna.pro$group_2<-"Cor.rna.pro"
  }else{
    
    data.dna.rna<-merge.data[,1:3]
    score.dna.rna<-merge(score,data.dna.rna[,1:2],by.x="genes",by.y="gene")
    score.dna.rna1<-score.dna.rna
    score.dna.rna$group<-ifelse(score.dna.rna[,2]>quantile(score.dna.rna[,2],0.7),"TOP30%",
                                ifelse(score.dna.rna[,2]<quantile(score.dna.rna[,2],0.3),"BOT30%",NA))
    
    score.dna.rna$Corum<-ifelse(score.dna.rna$corum==T,"Corum","No_Corum")
    corum.top30<-quantile(score.dna.rna[score.dna.rna$Corum=="Corum",][,2],0.7,na.rm=T)
    corum.bot30<-quantile(score.dna.rna[score.dna.rna$Corum=="Corum",][,2],0.3,na.rm=T)
    nocorum.top30<-quantile(score.dna.rna[score.dna.rna$Corum=="No_Corum",][,2],0.7,na.rm=T)
    nocorum.bot30<-quantile(score.dna.rna[score.dna.rna$Corum=="No_Corum",][,2],0.3,na.rm=T)
    score.dna.rna$group_1<-ifelse(score.dna.rna$Corum=="Corum" & score.dna.rna[,2]>corum.top30,"Corum.TOP30%",
                                  ifelse(score.dna.rna$Corum=="Corum" & score.dna.rna[,2]<corum.bot30,"Corum.BOT30%",
                                         ifelse(score.dna.rna$Corum=="No_Corum" & score.dna.rna[,2]>nocorum.top30,"No_Corum.TOP30%",
                                                ifelse(score.dna.rna$Corum=="No_Corum" & score.dna.rna[,2]<corum.bot30,"No_Corum.BOT30%",NA))))
    
    score.dna.rna<-na.omit(score.dna.rna)
    score.dna.rna$group_2<-"Cor.DNA.RNA"
    
    ### rna pro
    data.rna.pro<-merge.data[,c(1,5)]
    score.rna.pro<-merge(score,data.rna.pro[,1:2],by.x="genes",by.y="gene")
    score.rna.pro1<-score.rna.pro
    score.rna.pro$group<-ifelse(score.rna.pro[,2]>quantile(score.rna.pro[,2],0.7),"TOP30%",
                                ifelse(score.rna.pro[,2]<quantile(score.rna.pro[,2],0.3),"BOT30%",NA))
    
    score.rna.pro$Corum<-ifelse(score.rna.pro$corum==T,"Corum","No_Corum")
    corum.top30<-quantile(score.rna.pro[score.rna.pro$Corum=="Corum",][,2],0.7,na.rm=T)
    corum.bot30<-quantile(score.rna.pro[score.rna.pro$Corum=="Corum",][,2],0.3,na.rm=T)
    nocorum.top30<-quantile(score.rna.pro[score.rna.pro$Corum=="No_Corum",][,2],0.7,na.rm=T)
    nocorum.bot30<-quantile(score.rna.pro[score.rna.pro$Corum=="No_Corum",][,2],0.3,na.rm=T)
    score.rna.pro$group_1<-ifelse(score.rna.pro$Corum=="Corum" & score.rna.pro[,2]>corum.top30,"Corum.TOP30%",
                                  ifelse(score.rna.pro$Corum=="Corum" & score.rna.pro[,2]<corum.bot30,"Corum.BOT30%",
                                         ifelse(score.rna.pro$Corum=="No_Corum" & score.rna.pro[,2]>nocorum.top30,"No_Corum.TOP30%",
                                                ifelse(score.rna.pro$Corum=="No_Corum" & score.rna.pro[,2]<corum.bot30,"No_Corum.BOT30%",NA))))
    
    score.rna.pro<-na.omit(score.rna.pro)
    score.rna.pro$group_2<-"Cor.rna.pro"
}

D1<-score.dna.rna
colnames(D1)[5]<-"corr"

C.top.dna.rna<-D1[D1$group_1=="Corum.TOP30%",][,5]
NOC.top.dna.rna<-D1[D1$group_1=="No_Corum.TOP30%",][,5]
all.top.dna.rna<-D1[D1$group=="TOP30%",][,5]
C.bot.dna.rna<-D1[D1$group_1=="Corum.BOT30%",][,5]
NOC.bot.dna.rna<-D1[D1$group_1=="No_Corum.BOT30%",][,5]
all.bot.dna.rna<-D1[D1$group=="BOT30%",][,5]

test1<-signif(bootstrap(D1[,c(1,5,6)],10000,"corr","group","TOP30%","BOT30%"),3)
test2<-signif(bootstrap(D1[,c(1,5,8)],10000,"corr","group_1","Corum.BOT30%","Corum.TOP30%"),3)
test3<-signif(bootstrap(D1[,c(1,5,8)],10000,"corr","group_1","No_Corum.BOT30%","No_Corum.TOP30%"),3)

D1$group<-factor(D1$group,levels = c("TOP30%","BOT30%"))
RNA_DNA1<-ggplot(D1, aes(corr, stat(density),color = group)) + geom_density(alpha=0.1,size = 2) + xlab("RNA correlation with DNA") + 
  theme(text = element_text(size = 20)) + geom_vline(data=D1[D1$group=="TOP30%",], aes(xintercept = median(corr)), color=col1, size = 1.5, linetype='dashed') + 
  geom_vline(data = D1[D1$group=="BOT30%",], aes(xintercept = median(corr)), color = col2, size = 1.5, linetype='dashed')+
  annotate("text", x = -0, y = 3, label =paste("FDR=",test1,""),size=5,hjust = 0)+
  annotate("text", x = -0, y = 2.6, label =paste(signif(median(all.top.dna.rna),3),""),size=5,color=col1,hjust = 0)+
  annotate("text", x = -0, y = 2.2, label =paste(signif(median(all.bot.dna.rna),3),""),size=5,color=col2,hjust = 0)+
  ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
  xlim(-0, 1)+ylim(0,3)+
  theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                        axis.title=element_text(size=14,face="bold"),
                        plot.title = element_text(size = 16, face = "bold"),
                        legend.text=element_text(size=10,face = "bold"),
                        legend.title=element_text(size=12,face = "bold"))+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))+
  geom_hline(yintercept=0, colour="white", size=2)

D2<-D1[D1$Corum=="Corum",]
D2$group_1<-ifelse(D2$group_1=="Corum.TOP30%","TOP30%","BOT30%")
D2$group_1<-factor(D2$group_1,levels = c("TOP30%","BOT30%"))
RNA_DNA2<-ggplot(D2, aes(corr, stat(density),color = group_1)) + geom_density(alpha=0.1,size = 2) + xlab("RNA correlation with DNA") + 
  theme(text = element_text(size = 20)) + geom_vline(data=D2[D2$group_1=="TOP30%",], aes(xintercept = median(corr)), color=col1, size = 1.5, linetype='dashed') + 
  geom_vline(data = D2[D2$group_1=="BOT30%",], aes(xintercept = median(corr)), color = col2, size = 1.5, linetype='dashed')+
  annotate("text", x = -0, y = 3, label =paste("FDR=",test2,""),size=5,hjust = 0)+
  annotate("text", x = -0, y = 2.6, label =paste(signif(median(C.top.dna.rna),3),""),size=5,color=col1,hjust = 0)+
  annotate("text", x = -0, y = 2.2, label =paste(signif(median(C.bot.dna.rna),3),""),size=5,color=col2,hjust = 0)+
  ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
  xlim(-0, 1)+ylim(0,3)+
  theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                        axis.title=element_text(size=14,face="bold"),
                        plot.title = element_text(size = 16, face = "bold"),
                        legend.text=element_text(size=10,face = "bold"),
                        legend.title=element_text(size=12,face = "bold"))+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))+
  geom_hline(yintercept=0, colour="white", size=2)

D3<-D1[D1$Corum=="No_Corum",]
D3$group_1<-ifelse(D3$group_1=="No_Corum.TOP30%","TOP30%","BOT30%")
D3$group_1<-factor(D3$group_1,levels = c("TOP30%","BOT30%"))
RNA_DNA3<-ggplot(D3, aes(corr, stat(density),color = group_1)) + geom_density(alpha=0.1,size = 2) + xlab("RNA correlation with DNA") + 
  theme(text = element_text(size = 20)) + geom_vline(data=D3[D3$group_1=="TOP30%",], aes(xintercept = median(corr)), color=col1, size = 1.5, linetype='dashed') + 
  geom_vline(data = D3[D3$group_1=="BOT30%",], aes(xintercept = median(corr)), color = col2, size = 1.5, linetype='dashed')+
  annotate("text", x = -0, y = 3, label =paste("FDR=",test3,""),size=5,hjust = 0)+
  annotate("text", x = -0, y = 2.6, label =paste(signif(median(NOC.top.dna.rna),3),""),size=5,color=col1,hjust = 0)+
  annotate("text", x = -0, y = 2.2, label =paste(signif(median(NOC.bot.dna.rna),3),""),size=5,color=col2,hjust = 0)+
  ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
  xlim(-0, 1)+ylim(0,3)+
  theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                        axis.title=element_text(size=14,face="bold"),
                        plot.title = element_text(size = 16, face = "bold"),
                        legend.text=element_text(size=10,face = "bold"),
                        legend.title=element_text(size=12,face = "bold"))+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))+
  geom_hline(yintercept=0, colour="white", size=2)



R1<-score.rna.pro
colnames(R1)[5]<-"corr"
C.top.rna.pro<-R1[R1$group_1=="Corum.TOP30%",][,5]
NOC.top.rna.pro<-R1[R1$group_1=="No_Corum.TOP30%",][,5]
all.top.rna.pro<-R1[R1$group=="TOP30%",][,5]
C.bot.rna.pro<-R1[R1$group_1=="Corum.BOT30%",][,5]
NOC.bot.rna.pro<-R1[R1$group_1=="No_Corum.BOT30%",][,5]
all.bot.rna.pro<-R1[R1$group=="BOT30%",][,5]

test1<-signif(bootstrap(R1[,c(1,5,6)],10000,"corr","group","TOP30%","BOT30%"),3)
test2<-signif(bootstrap(R1[,c(1,5,8)],10000,"corr","group_1","Corum.BOT30%","Corum.TOP30%"),3)
test3<-signif(bootstrap(R1[,c(1,5,8)],10000,"corr","group_1","No_Corum.BOT30%","No_Corum.TOP30%"),3)

R1$group<-factor(R1$group,levels = c("TOP30%","BOT30%"))
RNA_PRO1<-ggplot(R1, aes(corr, stat(density),color = group)) + geom_density(alpha=0.1,size = 2) + xlab("RNA correlation with Protein") + 
  theme(text = element_text(size = 20)) + geom_vline(data=R1[R1$group=="TOP30%",], aes(xintercept = median(corr)), color=col1, size = 1.5, linetype='dashed') + 
  geom_vline(data = R1[R1$group=="BOT30%",], aes(xintercept = median(corr)), color = col2, size = 1.5, linetype='dashed')+
  annotate("text", x = -0, y = 3, label =paste("FDR=",test1,""),size=5,hjust = 0)+
  annotate("text", x = -0, y = 2.6, label =paste(signif(median(all.top.rna.pro),3),""),size=5,color=col1,hjust = 0)+
  annotate("text", x = -0, y = 2.2, label =paste(signif(median(all.bot.rna.pro),3),""),size=5,color=col2,hjust = 0)+
  ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
  xlim(-0, 1)+ylim(0,3)+
  theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                        axis.title=element_text(size=14,face="bold"),
                        plot.title = element_text(size = 16, face = "bold"),
                        legend.text=element_text(size=10,face = "bold"),
                        legend.title=element_text(size=12,face = "bold"))+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))+
  geom_hline(yintercept=0, colour="white", size=2)

R2<-R1[R1$Corum=="Corum",]
R2$group_1<-ifelse(R2$group_1=="Corum.TOP30%","TOP30%","BOT30%")
R2$group_1<-factor(R2$group_1,levels = c("TOP30%","BOT30%"))
RNA_PRO2<-ggplot(R2, aes(corr, stat(density),color = group_1)) + geom_density(alpha=0.1,size = 2) + xlab("RNA correlation with Protein") + 
  theme(text = element_text(size = 20)) + geom_vline(data=R2[R2$group_1=="TOP30%",], aes(xintercept = median(corr)), color=col1, size = 1.5, linetype='dashed') + 
  geom_vline(data = R2[R2$group_1=="BOT30%",], aes(xintercept = median(corr)), color = col2, size = 1.5, linetype='dashed')+
  annotate("text", x = -0, y = 3, label =paste("FDR=",test2,""),size=5,hjust = 0)+
  annotate("text", x = -0, y = 2.6, label =paste(signif(median(C.top.rna.pro),3),""),size=5,color=col1,hjust = 0)+
  annotate("text", x = -0, y = 2.2, label =paste(signif(median(C.bot.rna.pro),3),""),size=5,color=col2,hjust = 0)+
  ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
  xlim(-0, 1)+ylim(0,3)+
  theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                        axis.title=element_text(size=14,face="bold"),
                        plot.title = element_text(size = 16, face = "bold"),
                        legend.text=element_text(size=10,face = "bold"),
                        legend.title=element_text(size=12,face = "bold"))+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))+
  geom_hline(yintercept=0, colour="white", size=2)

R3<-R1[R1$Corum=="No_Corum",]
R3$group_1<-ifelse(R3$group_1=="No_Corum.TOP30%","TOP30%","BOT30%")
R3$group_1<-factor(R3$group_1,levels = c("TOP30%","BOT30%"))
RNA_PRO3<-ggplot(R3, aes(corr, stat(density),color = group_1)) + geom_density(alpha=0.1,size = 2) + xlab("RNA correlation with Protein") + 
  theme(text = element_text(size = 20)) + geom_vline(data=R3[R3$group_1=="TOP30%",], aes(xintercept = median(corr)), color=col1, size = 1.5, linetype='dashed') + 
  geom_vline(data = R3[R3$group_1=="BOT30%",], aes(xintercept = median(corr)), color = col2, size = 1.5, linetype='dashed')+
  annotate("text", x = -0, y = 3, label =paste("FDR=",test3,""),size=5,hjust = 0)+
  annotate("text", x = -0, y = 2.6, label =paste(signif(median(NOC.top.rna.pro),3),""),size=5,color=col1,hjust = 0)+
  annotate("text", x = -0, y = 2.2, label =paste(signif(median(NOC.bot.rna.pro),3),""),size=5,color=col2,hjust = 0)+
  ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
  xlim(-0, 1)+ylim(0,3)+
  theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                        axis.title=element_text(size=14,face="bold"),
                        plot.title = element_text(size = 16, face = "bold"),
                        legend.text=element_text(size=10,face = "bold"),
                        legend.title=element_text(size=12,face = "bold"))+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))+
  geom_hline(yintercept=0, colour="white", size=2)

#box plot
score.com<-unique(rbind(score.dna.rna1[,c(1,2,4)],score.rna.pro1[,c(1,2,4)]))
score.com$Corum<-ifelse(score.com$corum==T,"Corum","No_Corum")
score.com$Corum<-factor(score.com$Corum,levels = c("No_Corum","Corum"))
stat.test1 <- score.com %>%
  wilcox_test(median.PhyloP.score ~ Corum) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test1 <- stat.test1 %>%
  add_xy_position(x = "Corum", dodge = 0.8)
stat.test1$y.position<-1.15

boot1<-signif(bootstrap(score.com,10000,"median.PhyloP.score","Corum","Corum","No_Corum"),3)
stat.test1$boot<-boot1


n_fun <- function(x){
  return(data.frame(y = 1*1.35,
                    label = paste0("N=",length(x))))
}

col3<-"#aaa900" #dark yellow
col4<-"#6500aa" #purple
g1<-ggplot(score.com, aes(x=Corum, y=median.PhyloP.score)) + 
  # geom_point(aes(fill = Corum, color=Corum), alpha=0.5,
  #            size = 2, shape = 21, position = position_jitterdodge()) +
  geom_violin(aes(fill = Corum))+
  geom_boxplot(aes(fill = Corum),alpha = 0.5,outlier.colour = NA,width=0.5)+
  stat_summary(fun.data = n_fun, geom = "text", 
               aes(group=Corum),
               hjust = 0.5, position = position_dodge(0.8),size=5) +
  scale_color_manual(values=c(alpha(col3,1),alpha(col4,1)))+
  scale_fill_manual(values=c(alpha(col3,1),alpha(col4,1))) +
  theme_classic2()+
  ggtitle(paste0(""))+
  ylab("median of PhyloP score")+xlab("")+scale_y_continuous(limits = c(-1, 1.4))+
  theme(axis.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text=element_text(size=10,face = "bold"),
        legend.title=element_text(size=12,face = "bold"))+
  theme(legend.position = "right")+stat_pvalue_manual(
    stat.test1,  label = "boot",tip.length = 0,size=5
  )
#####rects plot
rects <- data.frame(y = 1:3,
                    colors = c("grey","grey","grey"),
                    text = paste(c("NoCorum","Corum","ALL")))

p3 <- ggplot(rects[1,], aes(x=0, y-1 , fill = colors,label=text)) +
  geom_tile(width = .1, height = .6) + # make square tiles
  geom_text(color = "grey",size=0.1) + # add white text in the middle
  scale_fill_identity(guide = "none") +# color the tiles with the colors in the data frame
  coord_fixed() + # make sure tiles are square
  theme_void()+ # remove any axis markings
  annotate(geom = "text", x = 0, y = rects$y[1]-1, label = rects$text[1], color = "white",size=7,
           angle = 90)
  p2 <- ggplot(rects[2,], aes(x=0, y-1 , fill = colors,label=text)) +
  geom_tile(width = .1, height = .6) + # make square tiles
  geom_text(color = "grey",size=0.1) + # add white text in the middle
  scale_fill_identity(guide = "none") +# color the tiles with the colors in the data frame
  coord_fixed() + # make sure tiles are square
  theme_void()+ # remove any axis markings
  annotate(geom = "text", x = 0, y = rects$y[2]-1, label = rects$text[2], color = "white",size=7,
           angle = 90)
  p1 <- ggplot(rects[3,], aes(x=0, y-1 , fill = colors,label=text)) +
    geom_tile(width = .1, height = .6) + # make square tiles
    geom_text(color = "grey",size=0.1) + # add white text in the middle
    scale_fill_identity(guide = "none") +# color the tiles with the colors in the data frame
    coord_fixed() + # make sure tiles are square
    theme_void()+ # remove any axis markings
    annotate(geom = "text", x = 0, y = rects$y[3]-1, label = rects$text[3], color = "white",size=7,
             angle = 90)

figure<-ggarrange(ggarrange(RNA_DNA1,RNA_PRO1,p1,nrow=1,ncol=3,common.legend=T,legend = "right",widths=c(1,1,0.1),heights = c(1,1,0.75)),
                  ggarrange(RNA_DNA2,RNA_PRO2,p2,nrow=1,ncol=3,common.legend=T,legend = "right",widths=c(1,1,0.1),heights = c(1,1,0.75)),
                  ggarrange(RNA_DNA3,RNA_PRO3,p3,nrow=1,ncol=3,common.legend=T,legend = "right",widths=c(1,1,0.1),heights = c(1,1,0.75)),
                  ggarrange(g1,NULL,NULL,nrow=1,ncol=3,common.legend=F,widths=c(1,0.6,0.1),heights = c(1,1,0.75)),
                  nrow=4,ncol = 1,heights =c(1,1,1,1.2))

pdf(paste0(names1[i],".phylopscore.pdf"),height=12,width=9)
print(figure)
dev.off() 
}









