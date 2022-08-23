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
col1<-"#aaa900" # buddha gold#"#FF0033" ###red
col2<-"#6500aa" # purple #"#3300FF" ## blue
setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes")
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
colon<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/Colon.corr.table.txt",sep="\t")
colnames(colon)[2]<-"colon.cor.rna.dna"
colnames(colon)[5]<-"colon.cor.rna.pro"
colnames(colon)[7]<-"colon.cor.dna.pro"
breast<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/BREAST.corr.table.txt",sep="\t")
colnames(breast)[2]<-"breast.cor.rna.dna"
colnames(breast)[5]<-"breast.cor.rna.pro"
colnames(breast)[7]<-"breast.cor.dna.pro"
ovarian<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/OV.corr.table.txt",sep="\t")
colnames(ovarian)[2]<-"ovarian.cor.rna.dna"
colnames(ovarian)[5]<-"ovarian.cor.rna.pro"
colnames(ovarian)[7]<-"ovarian.cor.dna.pro"
ccRCC<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/ccRCC.corr.table.txt",sep = "\t")
colnames(ccRCC)[2]<-"ccRCC.cor.rna.dna"
colnames(ccRCC)[5]<-"ccRCC.cor.rna.pro"
colnames(ccRCC)[7]<-"ccRCC.cor.dna.pro"
endometrial<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/Endometrial.corr.table.txt",sep = "\t")
colnames(endometrial)[2]<-"endometrial.cor.rna.dna"
colnames(endometrial)[5]<-"endometrial.cor.rna.pro"
colnames(endometrial)[7]<-"endometrial.cor.dna.pro"
luad<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/LUAD.corr.table.txt",sep = "\t")
colnames(luad)[2]<-"luad.cor.rna.dna"
colnames(luad)[5]<-"luad.cor.rna.pro"
colnames(luad)[7]<-"luad.cor.dna.pro"
hnsc<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/HNSC.corr.table.txt",sep = "\t")
colnames(hnsc)[2]<-"hnsc.cor.rna.dna"
colnames(hnsc)[5]<-"hnsc.cor.rna.pro"
colnames(hnsc)[7]<-"hnsc.cor.dna.pro"
score<-na.omit(read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/36.evolution/hg19.gene.score.txt",sep="\t",header=T))
score$corum<-ifelse(score$genes%in%master_list_names,TRUE,FALSE)


merge.data<-merge(merge(merge(merge(merge(merge(colon,breast,by="gene",all = T),ovarian,by="gene",all = T),ccRCC,by="gene",all = T),
                              endometrial,by="gene",all = T),luad,by="gene",all = T),hnsc,by="gene",all = T)
Pan.dna.rna<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% "cor.rna.dna"])
Pan.rna.pro<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% "cor.rna.pro"])
Pan.dna.pro<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% "cor.dna.pro"])

Pan.dna.rna$corr<-apply(Pan.dna.rna[,2:8],1,function(x)mean(x,na.rm=T))
Pan.rna.pro$corr<-apply(Pan.rna.pro[,2:8],1,function(x)mean(x,na.rm=T))
Pan.dna.pro$corr<-apply(Pan.dna.pro[,2:8],1,function(x)mean(x,na.rm=T))

Pan.dna.rna$corum<-ifelse(Pan.dna.rna$`merge.data$gene`%in%master_list_names,TRUE,FALSE)
colnames(Pan.dna.rna)[1]<-"gene"
Pan.rna.pro$corum<-ifelse(Pan.rna.pro$`merge.data$gene`%in%master_list_names,TRUE,FALSE)
colnames(Pan.rna.pro)[1]<-"gene"
Pan.dna.pro$corum<-ifelse(Pan.dna.pro$`merge.data$gene`%in%master_list_names,TRUE,FALSE)
colnames(Pan.dna.pro)[1]<-"gene"



cancer<-"Pan-Cancer"
CPTAC<-"CPTAC"

#plot
corums_cors <- Pan.dna.rna[Pan.dna.rna$corum == TRUE, "corr"]
noncorums_cors <- Pan.dna.rna[Pan.dna.rna$corum == FALSE, "corr"]
p.boot1<-p.adjust(signif(bootstrap(Pan.dna.rna,10000,"corr","corum",TRUE,FALSE),3),n=3,method = "fdr")

RNA_DNA<-ggplot(Pan.dna.rna, aes(corr, stat(density),color = corum)) + geom_density(alpha=0.1,size = 2) + xlab("RNA correlation with DNA") + 
  theme(text = element_text(size = 20)) + geom_vline(data=Pan.dna.rna[Pan.dna.rna$corum==TRUE,], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
  geom_vline(data = Pan.dna.rna[Pan.dna.rna$corum==FALSE,], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
  annotate("text", x = -0.5, y = 3, label =paste("FDR=",p.boot1,""),size=5,hjust = 0)+
  annotate("text", x = -0.5, y = 2.6, label =paste(signif(median(noncorums_cors),3),""),size=5,color=col1,hjust = 0)+
  annotate("text", x = -0.5, y = 2.2, label =paste(signif(median(corums_cors),3),""),size=5,color=col2,hjust = 0)+
  ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
  xlim(-0.5, 1)+ylim(0,3)+
  theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                        axis.title=element_text(size=14,face="bold"),
                        plot.title = element_text(size = 16, face = "bold"),
                        legend.text=element_text(size=10,face = "bold"),
                        legend.title=element_text(size=12,face = "bold"))+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))+
  geom_hline(yintercept=0, colour="white", size=2)

#plot
corums_cors <- Pan.rna.pro[Pan.rna.pro$corum == TRUE, "corr"]
noncorums_cors <- Pan.rna.pro[Pan.rna.pro$corum == FALSE, "corr"]
p.boot1<-p.adjust(signif(bootstrap(Pan.rna.pro,10000,"corr","corum",TRUE,FALSE),3),n=3,method = "fdr")
RNA_Protein<-ggplot(Pan.rna.pro, aes(corr, stat(density),color = corum)) + geom_density(alpha=0.1,size = 2) + xlab("RNA correlation with Protein") + 
  theme(text = element_text(size = 20)) + geom_vline(data=Pan.rna.pro[Pan.rna.pro$corum==TRUE,], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
  geom_vline(data = Pan.rna.pro[Pan.rna.pro$corum==FALSE,], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
  annotate("text", x = -0.5, y = 3, label =paste("FDR=",p.boot1,""),size=5,hjust = 0)+
  annotate("text", x = -0.5, y = 2.6, label =paste(signif(median(noncorums_cors),3),""),size=5,color=col1,hjust = 0)+
  annotate("text", x = -0.5, y = 2.2, label =paste(signif(median(corums_cors),3),""),size=5,color=col2,hjust = 0)+
  ggtitle(paste(cancer,"_",CPTAC,sep=""))+
  xlim(-0.5, 1)+ylim(0,3)+
  theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                        axis.title=element_text(size=14,face="bold"),
                        plot.title = element_text(size = 16, face = "bold"),
                        legend.text=element_text(size=10,face = "bold"),
                        legend.title=element_text(size=12,face = "bold"))+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))+
  geom_hline(yintercept=0, colour="white", size=2)

###PLOT
corums_cors <- Pan.dna.pro[Pan.dna.pro$corum == TRUE, "corr"]
noncorums_cors <- Pan.dna.pro[Pan.dna.pro$corum == FALSE, "corr"]
p.boot1<-p.adjust(signif(bootstrap(Pan.dna.pro,10000,"corr","corum",TRUE,FALSE),3),n=3,method = "fdr")

DNA_Protein<-ggplot(Pan.dna.pro, aes(corr, stat(density),color = corum)) + geom_density(alpha=0.1,size = 2) + xlab("DNA correlation with Protein") + 
  theme(text = element_text(size = 20)) + geom_vline(data=Pan.dna.pro[Pan.dna.pro$corum==TRUE,], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
  geom_vline(data = Pan.dna.pro[Pan.dna.pro$corum==FALSE,], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
  annotate("text", x = -0.5, y = 3, label =paste("FDR=",p.boot1,""),size=5,hjust = 0)+
  annotate("text", x = -0.5, y = 2.6, label =paste(signif(median(noncorums_cors),3),""),size=5,color=col1,hjust = 0)+
  annotate("text", x = -0.5, y = 2.2, label =paste(signif(median(corums_cors),3),""),size=5,color=col2,hjust = 0)+
  ggtitle(paste(cancer,"_",CPTAC,sep=""))+
  xlim(-0.5, 1)+ylim(0,3)+
  theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                        axis.title=element_text(size=14,face="bold"),
                        plot.title = element_text(size = 16, face = "bold"),
                        legend.text=element_text(size=10,face = "bold"),
                        legend.title=element_text(size=12,face = "bold"))+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))+
  geom_hline(yintercept=0, colour="white", size=2)


merge.temp<-merge(Pan.dna.rna[,c(1,9,10)],Pan.rna.pro[,c(1,9)],by="gene") 
merge.data<-merge(merge.temp,Pan.dna.pro[,c(1,9)],by="gene")
colnames(merge.data)[2]<-"corr.rna.dna"
colnames(merge.data)[4]<-"corr.rna.pro"
colnames(merge.data)[5]<-"cor.dna.pro"
write.table(merge.data,paste0(cancer,".corr.table.txt"),sep = "\t",row.names = F,quote = F)

figure<-ggarrange(RNA_DNA,RNA_Protein,DNA_Protein, ncol=3,nrow=1,common.legend = TRUE,legend = "right")
library(gridExtra)
library(ggpubr)
pdf(paste0(cancer,"_",CPTAC,"_DNA_RNA_PRO.pdf"),onefile=FALSE,height=3,width=10)
print(figure)
dev.off() 





