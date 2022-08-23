`%notin%` <- Negate(`%in%`)
set.seed(1234)
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
library(tidyr)
library(data.table)
library(ggplot2)
library(boot)
library(ggpubr)
library(rstatix)
library(clusterProfiler)
library(msigdbr)
setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig3/210911.ribo.pr.mt.pm/")

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
#######################################################################################################
sublocation<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/global_database/location_cell/subcellular_location.tsv",sep="\t",header = T)
subloc<-separate_rows(sublocation, Main.location,sep=";")
mt.loc<-c("Mitochondria")
plasma_membrane.loc<-c("Plasma membrane")
mito.list<-as.data.frame(unique(subloc[subloc$Main.location %in% mt.loc,]))
plasma_membrane.list<-as.data.frame(unique(subloc[subloc$Main.location %in% plasma_membrane.loc,]))

ribo.list<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210824.ribo/ribosomal.txt",
                      sep="\t",header = T)
proteasome.list<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210824.ribo/proteasome.txt",
                            sep="\t",header = T)

#ribo.PR.list<-as.data.frame(rbind(proteasome.list,ribo.list))

########################################################################################################
## generate pan-cancer data
###############################3
path1<-c("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210812_all_exp_genes/")

colon<-read.delim(paste0(path1,"/Colon.corr.table.txt"),sep="\t")
colnames(colon)[2]<-"colon.cor.rna.dna"
colnames(colon)[5]<-"colon.cor.rna.pro"
breast<-read.delim(paste0(path1,"/BREAST.corr.table.txt"),sep="\t")
colnames(breast)[2]<-"breast.cor.rna.dna"
colnames(breast)[5]<-"breast.cor.rna.pro"
ovarian<-read.delim(paste0(path1,"/OV.corr.table.txt"),sep="\t")
colnames(ovarian)[2]<-"ovarian.cor.rna.dna"
colnames(ovarian)[5]<-"ovarian.cor.rna.pro"
ccRCC<-read.delim(paste0(path1,"/ccRCC.corr.table.txt"),sep = "\t")
colnames(ccRCC)[2]<-"ccRCC.cor.rna.dna"
colnames(ccRCC)[5]<-"ccRCC.cor.rna.pro"
endometrial<-read.delim(paste0(path1,"/Endometrial.corr.table.txt"),sep = "\t")
colnames(endometrial)[2]<-"endometrial.cor.rna.dna"
colnames(endometrial)[5]<-"endometrial.cor.rna.pro"
luad<-read.delim(paste0(path1,"/LUAD.corr.table.txt"),sep = "\t")
colnames(luad)[2]<-"luad.cor.rna.dna"
colnames(luad)[5]<-"luad.cor.rna.pro"
hnsc<-read.delim(paste0(path1,"/HNSC.corr.table.txt"),sep = "\t")
colnames(hnsc)[2]<-"hnsc.cor.rna.dna"
colnames(hnsc)[5]<-"hnsc.cor.rna.pro"
score<-na.omit(read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/36.evolution/hg19.gene.score.txt",sep="\t",header=T))
score$corum<-ifelse(score$genes%in%master_list_names,TRUE,FALSE)


merge.data<-merge(merge(merge(merge(merge(merge(colon,breast,by="gene",all = T),ovarian,by="gene",all = T),ccRCC,by="gene",all = T),
                              endometrial,by="gene",all = T),luad,by="gene",all = T),hnsc,by="gene",all = T)
Pan.dna.rna<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% "cor.rna.dna"])
Pan.rna.pro<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% "cor.rna.pro"])

Pan.dna.rna$max.rna.dna<-apply(Pan.dna.rna[,2:8],1,function(x)mean(x,na.rm=T))
Pan.rna.pro$max.rna.pro<-apply(Pan.rna.pro[,2:8],1,function(x)mean(x,na.rm=T))

Pan.data<-merge(Pan.dna.rna,Pan.rna.pro,by="merge.data$gene",all = T)
colnames(Pan.data)[1]<-"gene"

names<-list(colon,breast,ovarian,ccRCC,endometrial,luad,hnsc,Pan.data)
names1<-c("COAD","BRCA","OV","ccRCC","UCEC","LUAD","HNSC","PAN-Cancer")
col1<-"#AAA900" ###budddha gold
col2<-"#6500aa" ## purple
col3<-"#bfc2c3" 
col4<-"#acafb0" 
col5<-"#999b9c" 
col6<-"#868889" 
col8.nc<-"#777600" ### er
#col7.nc<-"#c4c34d" ### pm
#col8<-"#3d0066" ## er
col7<-"#a366cc" ###pm
col9<-"#d7dbdb"
cancer<-"CPTAC"


for (i in 1:length(names)){
  #i<-1
  
  merge.data<-names[[i]]
  CPTAC<-names1[i]
  if(names1[i]=="PAN-Cancer"){
    data.dna.rna<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% ".rna.dna"]) 
    colnames(data.dna.rna)[1]<-"genes" 
    data.dna.rna$corum<-ifelse(data.dna.rna$genes%in%master_list_names,TRUE,FALSE)
    data.dna.rna<-data.dna.rna[,c("genes","max.rna.dna","corum")]
    data.rna.pro<-cbind(merge.data$gene, merge.data[,colnames(merge.data) %like% ".rna.pro"]) 
    colnames(data.rna.pro)[1]<-"genes" 
    data.rna.pro$corum<-ifelse(data.rna.pro$genes%in%master_list_names,TRUE,FALSE)
    data.rna.pro<-data.rna.pro[,c("genes","max.rna.pro","corum")]
  }else{
    data.dna.rna<-merge.data[,1:3]
    colnames(data.dna.rna)[3]<-"corum"
    data.rna.pro<-merge.data[,c(1,5,3)]
    colnames(data.rna.pro)[3]<-"corum"
  }
    D<-data.dna.rna
    colnames(D)[2]<-"corr"
    
    D1<-D
    D1$location<-ifelse(D1$corum=="TRUE","Corum:all","NoCorum")
    ##for ER only
    DD.ribo<-D[D$gene %in% ribo.list$Gene.name & D$corum==TRUE,]
    DD.pr<-D[D$gene %in% proteasome.list$Gene.name & D$corum==TRUE,]
    DD.pm<-D[D$gene %in% plasma_membrane.list$Gene.name & D$corum==TRUE,]
    DD.mito<-D[D$gene %in% mito.list$Gene.name & D$corum==TRUE,]
    ##
    
    corum.ribo<-D1[D1$gene%in%ribo.list$Gene.name & D1$gene %in% master_list_names ,]
    corum.ribo$location<-"Corum:Ribosome"
    corum.pr<-D1[D1$gene%in%proteasome.list$Gene.name & D1$gene %in% master_list_names ,]
    corum.pr$location<-"Corum:Proteasome"
    corum.PM<-D1[D1$gene%in%plasma_membrane.list$Gene.name & D1$gene %in% master_list_names,]
    corum.PM$location<-"Corum:plasma_membrane"
    corum.mito<-D1[D1$gene%in%mito.list$Gene.name & D1$gene %in% master_list_names,]
    corum.mito$location<-"Corum:Mitochondria"
    
    cor_rna_dna_colon<-rbind(D1,corum.ribo,corum.mito,corum.PM,corum.pr)

    cor_rna_dna_colon$location<-factor(cor_rna_dna_colon$location, levels=c("NoCorum", "Corum:all",
                                                                            "Corum:Ribosome",
                                                                            "Corum:Proteasome",
                                                                            "Corum:Mitochondria",
                                                                            "Corum:plasma_membrane"))
    corums_cors <- cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:all",]
    noncorums_cors <- cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum",]
     ribo_cors<-cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:Ribosome",]
     pr_cors<-cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:Proteasome",]
     mt_cors<-cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:Mitochondria",]
     PM_cors<-cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:plasma_membrane",]
     
     dataset1<-as.data.frame(rbind(ribo_cors,corums_cors))
     dataset2<-as.data.frame(rbind(pr_cors,corums_cors))
     dataset3<-as.data.frame(rbind(mt_cors,corums_cors))
     dataset4<-as.data.frame(rbind(PM_cors,corums_cors))
     
     set.seed(1234)
     test1<-p.adjust(bootstrap(dataset1,10000,"corr","location","Corum:all","Corum:Ribosome"),method = "fdr",n = 2)
     set.seed(1234)
     test2<-p.adjust(bootstrap(dataset2,10000,"corr","location","Corum:all","Corum:Proteasome"),method = "fdr",n = 2)
     set.seed(1234)
     test3<-p.adjust(bootstrap(dataset3,10000,"corr","location","Corum:all","Corum:Mitochondria"),method = "fdr",n = 2)
     set.seed(1234)
     test4<-p.adjust(bootstrap(dataset4,10000,"corr","location","Corum:all","Corum:plasma_membrane"),method = "fdr",n = 2)
     
     cor_rna_dna_colon.1<-cor_rna_dna_colon[cor_rna_dna_colon$location %in% c("NoCorum","Corum:Ribosome","Corum:all"),]
     RNA_DNA_C.1<-ggplot(cor_rna_dna_colon.1, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with DNA") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:Ribosome",], aes(xintercept = median(corr)), color = col7, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       
    #   
     annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test1,3)),size=6,color=col7,hjust = 0)+
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col7))+
       scale_color_manual(values=c(col1,col2,col7))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     cor_rna_dna_colon.2<-cor_rna_dna_colon[cor_rna_dna_colon$location %in% c("NoCorum","Corum:Proteasome","Corum:all"),]
     RNA_DNA_C.2<-ggplot(cor_rna_dna_colon.2, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with DNA") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:Proteasome",], aes(xintercept = median(corr)), color = col7, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       
       #   
       annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test2,3)),size=6,color=col7,hjust = 0)+
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col7))+
       scale_color_manual(values=c(col1,col2,col7))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     
     
     
     cor_rna_dna_colon.3<-cor_rna_dna_colon[cor_rna_dna_colon$location %in%c("NoCorum","Corum:Mitochondria","Corum:all"),]
     RNA_DNA_C.3<-ggplot(cor_rna_dna_colon.3, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with DNA") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:Mitochondria",], aes(xintercept = median(corr)), color = col7, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       
       #   
       annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test3,3)),size=6,color=col7,hjust = 0)+
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col7))+
       scale_color_manual(values=c(col1,col2,col7))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     
     cor_rna_dna_colon.4<-cor_rna_dna_colon[cor_rna_dna_colon$location %in%c("NoCorum","Corum:plasma_membrane","Corum:all"),]
     RNA_DNA_C.4<-ggplot(cor_rna_dna_colon.4, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with DNA") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:plasma_membrane",], aes(xintercept = median(corr)), color = col7, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       
       #   
       annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test4,3)),size=6,color=col7,hjust = 0)+
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col7))+
       scale_color_manual(values=c(col1,col2,col7))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     

     D2<-D
     D2$location<-ifelse(D2$corum=="TRUE","Corum:all","NoCorum")
     n.corum.ribo<-D2[D2$gene%in%ribo.list$Gene.name & D2$gene %notin% master_list_names,]
     #n.corum.ribo$location<-"NoCorum:Ribosome"
     n.corum.pr<-D2[D2$gene%in%proteasome.list$Gene.name & D2$gene %notin% master_list_names,]
     n.corum.pr$location<-"NoCorum:Proteasome"
     n.corum.PM<-D2[D2$gene%in%plasma_membrane.list$Gene.name & D2$gene %notin% master_list_names,]
     n.corum.PM$location<-"NoCorum:plasma_membrane"
     n.corum.mito<-D2[D2$gene%in%mito.list$Gene.name & D2$gene %notin% master_list_names,]
     n.corum.mito$location<-"NoCorum:Mitochondria"
     
     
     
     cor_rna_dna_colon<-rbind(D2,n.corum.ribo,n.corum.PM,
                              n.corum.mito,n.corum.pr)
     cor_rna_dna_colon$location<-factor(cor_rna_dna_colon$location, levels=c("NoCorum", "Corum:all","NoCorum:nucleus",
                                                                             "NoCorum:Ribosome",
                                                                             "NoCorum:Proteasome",
                                                                             "NoCorum:Mitochondria", 
                                                                             "NoCorum:plasma_membrane"))
     
     n.corums_cors <- cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:all",]
     n.noncorums_cors <- cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum",]
     n.ribo_cors<-cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum:Ribosome",]
     n.pr_cors<-cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum:Proteasome",]
     n.mt_cors<-cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum:Mitochondria",]
     n.PM_cors<-cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum:plasma_membrane",]
     
     dataset5<-as.data.frame(rbind(n.ribo_cors,n.noncorums_cors))
     dataset6<-as.data.frame(rbind(n.pr_cors,n.noncorums_cors))
     dataset7<-as.data.frame(rbind(n.mt_cors,n.noncorums_cors))
     dataset8<-as.data.frame(rbind(n.PM_cors,n.noncorums_cors))
     
     set.seed(1234)
     test5<-p.adjust(bootstrap(dataset5,10000,"corr","location","NoCorum","NoCorum:Ribosome"),method = "fdr",n = 2)
     set.seed(1234)
     test6<-p.adjust(bootstrap(dataset6,10000,"corr","location","NoCorum","NoCorum:Proteasome"),method = "fdr",n = 2)
     set.seed(1234)
     test7<-p.adjust(bootstrap(dataset7,10000,"corr","location","NoCorum","NoCorum:Mitochondria"),method = "fdr",n = 2)
     set.seed(1234)
     test8<-p.adjust(bootstrap(dataset8,10000,"corr","location","NoCorum","NoCorum:plasma_membrane"),method = "fdr",n = 2)
     
     
     cor_rna_dna_colon.5<-cor_rna_dna_colon[cor_rna_dna_colon$location %in%c("NoCorum","NoCorum:Ribosome","Corum:all"),]
     RNA_DNA_NC.5<-ggplot(cor_rna_dna_colon.5, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with DNA") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum:Ribosome",], aes(xintercept = median(corr)), color = col8.nc, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       
       #   
     annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test5,3)),size=6,color=col8.nc,hjust = 0)+
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col8.nc))+
       scale_color_manual(values=c(col1,col2,col8.nc))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     cor_rna_dna_colon.6<-cor_rna_dna_colon[cor_rna_dna_colon$location %in%c("NoCorum","NoCorum:Proteasome","Corum:all"),]
     RNA_DNA_NC.6<-ggplot(cor_rna_dna_colon.6, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with DNA") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum:Proteasome",], aes(xintercept = median(corr)), color = col8.nc, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       
       #   
       annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test6,3)),size=6,color=col8.nc,hjust = 0)+
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col8.nc))+
       scale_color_manual(values=c(col1,col2,col8.nc))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     cor_rna_dna_colon.7<-cor_rna_dna_colon[cor_rna_dna_colon$location %in%c("NoCorum","NoCorum:Mitochondria","Corum:all"),]
     RNA_DNA_NC.7<-ggplot(cor_rna_dna_colon.7, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with DNA") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum:Mitochondria",], aes(xintercept = median(corr)), color = col8.nc, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       
       #   
       annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test7,3)),size=6,color=col8.nc,hjust = 0)+
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col8.nc))+
       scale_color_manual(values=c(col1,col2,col8.nc))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     cor_rna_dna_colon.8<-cor_rna_dna_colon[cor_rna_dna_colon$location %in%c("NoCorum","NoCorum:plasma_membrane","Corum:all"),]
     RNA_DNA_NC.8<-ggplot(cor_rna_dna_colon.8, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with DNA") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_dna_colon[cor_rna_dna_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_dna_colon[cor_rna_dna_colon$location=="NoCorum:plasma_membrane",], aes(xintercept = median(corr)), color = col8.nc, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       
       #   
       annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test8,3)),size=6,color=col8.nc,hjust = 0)+
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col8.nc))+
       scale_color_manual(values=c(col1,col2,col8.nc))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     
     #############    RNA PROTEIN
     R<-data.rna.pro
     colnames(R)[2]<-"corr"

     R1<-R
     R1$location<-ifelse(R1$corum=="TRUE","Corum:all","NoCorum")
     ##for ribo and proteasome only
     RR.ribo<-R[R$gene %in% ribo.list$Gene.name & R$corum==TRUE,]
     RR.ribo$location<-"Ribosome and proteasome"
     name<-colnames(DD.ribo)[1]
     DR.ribo<-merge(DD.ribo,RR.ribo,by=name)
     colnames(DR.ribo)[2]<-"cor.dna.rna"
     colnames(DR.ribo)[4]<-"cor.rna.pro"
     write.table(DR.ribo,paste0(names1[i],".ribo.corr.txt"),sep="\t",row.names = F,quote = F)
     genenames<-as.character(DR.ribo[,1])
     m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% #subcategory = "BP"
       dplyr::select(gs_name, gene_symbol)
     ribo.go<- enricher(genenames, TERM2GENE=m_t2g)
     results<-ribo.go@result
     write.table(results,paste0(names1[i],".ribo.enrichment.txt"),sep="\t",row.names = F,quote = F)
     #for MT pathway
     RR.mito<-R[R$gene %in% mito.list$Gene.name & R$corum==TRUE,]
     RR.mito$location<-"mito"
     name<-colnames(DD.mito)[1]
     DR.mito<-merge(DD.mito,RR.mito,by=name)
     colnames(DR.mito)[2]<-"cor.dna.rna"
     colnames(DR.mito)[4]<-"cor.rna.pro"
     write.table(DR.mito,paste0(names1[i],".mito.corr.txt"),sep="\t",row.names = F,quote = F)
     genenames<-as.character(DR.mito[,1])
     m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% #subcategory = "BP"
       dplyr::select(gs_name, gene_symbol)
     mito<- enricher(genenames, TERM2GENE=m_t2g)
     results<-mito@result
     write.table(results,paste0(names1[i],".mito.enrichment.txt"),sep="\t",row.names = F,quote = F)
     ## for proteasome only
     RR.pr<-R[R$gene %in% proteasome.list$Gene.name & R$corum==TRUE,]
     RR.pr$location<-"pr"
     name<-colnames(DD.pr)[1]
     DR.pr<-merge(DD.pr,RR.pr,by=name)
     colnames(DR.pr)[2]<-"cor.dna.rna"
     colnames(DR.pr)[4]<-"cor.rna.pro"
     write.table(DR.pr,paste0(names1[i],".pr.corr.txt"),sep="\t",row.names = F,quote = F)
     genenames<-as.character(DR.pr[,1])
     m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% #subcategory = "BP"
       dplyr::select(gs_name, gene_symbol)
     pr<- enricher(genenames, TERM2GENE=m_t2g)
     results<-pr@result
     write.table(results,paste0(names1[i],".pr.enrichment.txt"),sep="\t",row.names = F,quote = F)
     
     # for PM pathway
     RR.pm<-R[R$gene %in% plasma_membrane.list$Gene.name & R$corum==TRUE,]
     RR.pm$location<-"PM"
     name<-colnames(DD.pm)[1]
     DR.pm<-merge(DD.pm,RR.pm,by=name)
     colnames(DR.pm)[2]<-"cor.dna.rna"
     colnames(DR.pm)[4]<-"cor.rna.pro"
     write.table(DR.pm,paste0(names1[i],".pm.corr.txt"),sep="\t",row.names = F,quote = F)
     genenames<-as.character(DR.pm[,1])
     m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% #subcategory = "BP"
       dplyr::select(gs_name, gene_symbol)
     pm<- enricher(genenames, TERM2GENE=m_t2g)
     results<-pm@result
     write.table(results,paste0(names1[i],".PM.enrichment.txt"),sep="\t",row.names = F,quote = F)
     
     ##
     
     corum.ribo<-R1[R1$gene%in%ribo.list$Gene.name & R1$gene %in% master_list_names ,]
     corum.ribo$location<-"Corum:Ribosome"
     corum.pr<-R1[R1$gene%in%proteasome.list$Gene.name & R1$gene %in% master_list_names ,]
     corum.pr$location<-"Corum:Proteasome"
     corum.PM<-R1[R1$gene%in%plasma_membrane.list$Gene.name & R1$gene %in% master_list_names,]
     corum.PM$location<-"Corum:plasma_membrane"
     corum.mito<-R1[R1$gene%in%mito.list$Gene.name & R1$gene %in% master_list_names,]
     corum.mito$location<-"Corum:Mitochondria"
     
     cor_rna_pro_colon<-rbind(R1,corum.ribo,corum.mito,corum.PM,corum.pr)
     
     cor_rna_pro_colon$location<-factor(cor_rna_pro_colon$location, levels=c("NoCorum", "Corum:all",
                                                                             "Corum:Ribosome",
                                                                             "Corum:Proteasome",
                                                                             "Corum:Mitochondria",
                                                                             "Corum:plasma_membrane"))
     corums_cors <- cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:all",]
     noncorums_cors <- cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum",]
     ribo_cors<-cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:Ribosome",]
     pr_cors<-cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:Proteasome",]
     mt_cors<-cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:Mitochondria",]
     PM_cors<-cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:plasma_membrane",]
     
     dataset11<-as.data.frame(rbind(ribo_cors,corums_cors))
     dataset22<-as.data.frame(rbind(pr_cors,corums_cors))
     dataset33<-as.data.frame(rbind(mt_cors,corums_cors))
     dataset44<-as.data.frame(rbind(PM_cors,corums_cors))
     
     set.seed(1234)
     test11<-p.adjust(bootstrap(dataset11,10000,"corr","location","Corum:all","Corum:Ribosome"),method = "fdr",n = 2)
     set.seed(1234)
     test22<-p.adjust(bootstrap(dataset22,10000,"corr","location","Corum:all","Corum:Proteasome"),method = "fdr",n = 2)
     set.seed(1234)
     test33<-p.adjust(bootstrap(dataset33,10000,"corr","location","Corum:all","Corum:Mitochondria"),method = "fdr",n = 2)
     set.seed(1234)
     test44<-p.adjust(bootstrap(dataset44,10000,"corr","location","Corum:all","Corum:plasma_membrane"),method = "fdr",n = 2)
     
     cor_rna_pro_colon.1<-cor_rna_pro_colon[cor_rna_pro_colon$location %in%c("NoCorum","Corum:Ribosome","Corum:all"),]
     RNA_Protein_C.1<-ggplot(cor_rna_pro_colon.1, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with Protein") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:Ribosome",], aes(xintercept = median(corr)), color = col7, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       
     annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test11,3)),size=6,color=col7,hjust = 0)+
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col7))+
       scale_color_manual(values=c(col1,col2,col7))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     cor_rna_pro_colon.2<-cor_rna_pro_colon[cor_rna_pro_colon$location %in%c("NoCorum","Corum:Proteasome","Corum:all"),]
     RNA_Protein_C.2<-ggplot(cor_rna_pro_colon.2, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with Protein") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:Proteasome",], aes(xintercept = median(corr)), color = col7, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       
       annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test22,3)),size=6,color=col7,hjust = 0)+
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col7))+
       scale_color_manual(values=c(col1,col2,col7))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     cor_rna_pro_colon.3<-cor_rna_pro_colon[cor_rna_pro_colon$location %in%c("NoCorum","Corum:Mitochondria","Corum:all"),]
     RNA_Protein_C.3<-ggplot(cor_rna_pro_colon.3, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with Protein") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:Mitochondria",], aes(xintercept = median(corr)), color = col7, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       
       annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test33,3)),size=6,color=col7,hjust = 0)+
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col7))+
       scale_color_manual(values=c(col1,col2,col7))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     cor_rna_pro_colon.4<-cor_rna_pro_colon[cor_rna_pro_colon$location %in%c("NoCorum","Corum:plasma_membrane","Corum:all"),]
     RNA_Protein_C.4<-ggplot(cor_rna_pro_colon.4, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with Protein") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:plasma_membrane",], aes(xintercept = median(corr)), color = col7, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       
       annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test44,3)),size=6,color=col7,hjust = 0)+
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"))+
       scale_fill_manual(values=c(col1,col2,col7))+
       scale_color_manual(values=c(col1,col2,col7))+
       geom_hline(yintercept=0, colour="white", size=2)

     R2<-R
     R2$location<-ifelse(R2$corum=="TRUE","Corum:all","NoCorum")
     n.corum.ribo<-R2[R2$gene%in%ribo.list$Gene.name & R2$gene %notin% master_list_names,]
     #n.corum.ribo$location<-"NoCorum:Ribosome"
     n.corum.pr<-R2[R2$gene%in%proteasome.list$Gene.name & R2$gene %notin% master_list_names,]
     n.corum.pr$location<-"NoCorum:Proteasome"
     n.corum.PM<-R2[R2$gene%in%plasma_membrane.list$Gene.name & R2$gene %notin% master_list_names,]
     n.corum.PM$location<-"NoCorum:plasma_membrane"
     n.corum.mito<-R2[R2$gene%in%mito.list$Gene.name & R2$gene %notin% master_list_names,]
     n.corum.mito$location<-"NoCorum:Mitochondria"
     
     
     
     cor_rna_pro_colon<-rbind(R2,n.corum.ribo,n.corum.PM,
                              n.corum.mito,n.corum.pr)
     cor_rna_pro_colon$location<-factor(cor_rna_pro_colon$location, levels=c("NoCorum", "Corum:all","NoCorum:nucleus",
                                                                             "NoCorum:Ribosome",
                                                                             "NoCorum:Proteasome",
                                                                             "NoCorum:Mitochondria", 
                                                                             "NoCorum:plasma_membrane"))
     
     n.corums_cors <- cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:all",]
     n.noncorums_cors <- cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum",]
     n.ribo_cors<-cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum:Ribosome",]
     n.pr_cors<-cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum:Proteasome",]
     n.mt_cors<-cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum:Mitochondria",]
     n.PM_cors<-cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum:plasma_membrane",]
     
    
     dataset55<-as.data.frame(rbind(n.ribo_cors,n.noncorums_cors))
     dataset66<-as.data.frame(rbind(n.pr_cors,n.noncorums_cors))
     dataset77<-as.data.frame(rbind(n.mt_cors,n.noncorums_cors))
     dataset88<-as.data.frame(rbind(n.PM_cors,n.noncorums_cors))
     
     set.seed(1234)
     test55<-p.adjust(bootstrap(dataset55,10000,"corr","location","NoCorum","NoCorum:Ribosome"),method = "fdr",n = 2)
     set.seed(1234)
     test66<-p.adjust(bootstrap(dataset66,10000,"corr","location","NoCorum","NoCorum:Proteasome"),method = "fdr",n = 2)
     set.seed(1234)
     test77<-p.adjust(bootstrap(dataset77,10000,"corr","location","NoCorum","NoCorum:Mitochondria"),method = "fdr",n = 2)
     set.seed(1234)
     test88<-p.adjust(bootstrap(dataset88,10000,"corr","location","NoCorum","NoCorum:plasma_membrane"),method = "fdr",n = 2)
     
     cor_rna_pro_colon.5<-cor_rna_pro_colon[cor_rna_pro_colon$location %in%c("NoCorum","NoCorum:Ribosome","Corum:all"),]
     RNA_Protein_NC.5<-ggplot(cor_rna_pro_colon.5, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with Protein") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum:Ribosome",], aes(xintercept = median(corr)), color = col8.nc, size = 1, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       #
     annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test55,3)),size=6,color=col8.nc,hjust = 0)+
     
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col8.nc))+
       scale_color_manual(values=c(col1,col2,col8.nc))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     cor_rna_pro_colon.6<-cor_rna_pro_colon[cor_rna_pro_colon$location %in%c("NoCorum","NoCorum:Proteasome","Corum:all"),]
     RNA_Protein_NC.6<-ggplot(cor_rna_pro_colon.6, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with Protein") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum:Proteasome",], aes(xintercept = median(corr)), color = col8.nc, size = 1, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       #
       annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test66,3)),size=6,color=col8.nc,hjust = 0)+
       
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col8.nc))+
       scale_color_manual(values=c(col1,col2,col8.nc))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     cor_rna_pro_colon.7<-cor_rna_pro_colon[cor_rna_pro_colon$location %in%c("NoCorum","NoCorum:Mitochondria","Corum:all"),]
     RNA_Protein_NC.7<-ggplot(cor_rna_pro_colon.7, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with Protein") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum:Mitochondria",], aes(xintercept = median(corr)), color = col8.nc, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       #
       annotate("text", x = -0.2, y = 2.9, label =paste0("FDR=",signif(test77,3)),size=6,color=col8.nc,hjust = 0)+
       
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,3)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col8.nc))+
       scale_color_manual(values=c(col1,col2,col8.nc))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     cor_rna_pro_colon.8<-cor_rna_pro_colon[cor_rna_pro_colon$location %in%c("NoCorum","NoCorum:plasma_membrane","Corum:all"),]
     RNA_Protein_NC.8<-ggplot(cor_rna_pro_colon.8, aes(corr, stat(density),color = location,linetype=location)) + geom_density(alpha=0.1,size = 1.5) + xlab("RNA correlation with Protein") + 
       theme(text = element_text(size = 20)) + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum",], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
       geom_vline(data=cor_rna_pro_colon[cor_rna_pro_colon$location=="Corum:all",], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
       geom_vline(data = cor_rna_pro_colon[cor_rna_pro_colon$location=="NoCorum:plasma_membrane",], aes(xintercept = median(corr)), color = col8.nc, size = 1.5, linetype='dashed')+
       scale_linetype_manual(name = "location", values = c("solid", "solid","dashed"))+
       #
       annotate("text", x = -0.2, y = 4.9, label =paste0("FDR=",signif(test88,3)),size=6,color=col8.nc,hjust = 0)+
       
       ggtitle(paste0(cancer,"_",CPTAC,sep=""))+
       xlim(-0.2, 1)+ylim(0,8)+
       theme_classic()+theme(axis.text=element_text(size=12,face = "bold"),
                             axis.title=element_text(size=14,face="bold"),
                             plot.title = element_text(size = 14, face = "bold"),
                             legend.text=element_text(size=10,face = "bold"),
                             legend.title=element_text(size=12,face = "bold"),
                             legend.position = "none")+
       scale_fill_manual(values=c(col1,col2,col8.nc))+
       scale_color_manual(values=c(col1,col2,col8.nc))+
       geom_hline(yintercept=0, colour="white", size=2)
     
     
     rects <- data.frame(y = 1:2,
                         colors = c("grey","grey"),
                         text = paste(c("NoCorum","Corum")))
     
     p2 <- ggplot(rects[1,], aes(x=0, y-1 , fill = colors,label=text)) +
       geom_tile(width = .1, height = .6) + # make square tiles
       geom_text(color = "grey",size=0.1) + # add white text in the middle
       scale_fill_identity(guide = "none") +# color the tiles with the colors in the data frame
       coord_fixed() + # make sure tiles are square
       theme_void()+ # remove any axis markings
       annotate(geom = "text", x = 0, y = rects$y[1]-1, label = rects$text[1], color = "white",size=7,
                angle = 90)
     p1 <- ggplot(rects[2,], aes(x=0, y-1 , fill = colors,label=text)) +
       geom_tile(width = .1, height = .6) + # make square tiles
       geom_text(color = "grey",size=0.1) + # add white text in the middle
       scale_fill_identity(guide = "none") +# color the tiles with the colors in the data frame
       coord_fixed() + # make sure tiles are square
       theme_void()+ # remove any axis markings
       annotate(geom = "text", x = 0, y = rects$y[2]-1, label = rects$text[2], color = "white",size=7,
                angle = 90)
     
     
     
library(gridExtra)
library("ggpubr")
#figure<-ggarrange(ggarrange(RNA_DNA_NC.5,RNA_Protein_NC.5,p2,nrow = 1,ncol = 3,common.legend = T,legend = "top",widths=c(1,1,0.1)),
#                  ggarrange(RNA_DNA_C.1,RNA_Protein_C.1,p1,nrow = 1,ncol = 3,common.legend = F,widths=c(1,1,0.1)),
#                  ggarrange(RNA_DNA_NC.6,RNA_Protein_NC.6,p2,nrow = 1,ncol = 3,common.legend = F,widths=c(1,1,0.1)),
#                  ggarrange(RNA_DNA_C.2,RNA_Protein_C.2,p1,nrow = 1,ncol = 3,common.legend = F,widths=c(1,1,0.1)),
##                  ggarrange(RNA_DNA_NC.7,RNA_Protein_NC.7,p2,nrow = 1,ncol = 3,common.legend = F,widths=c(1,1,0.1)),
#                  ggarrange(RNA_DNA_C.3,RNA_Protein_C.3,p1,nrow = 1,ncol = 3,common.legend = F,widths=c(1,1,0.1)),
#                  ggarrange(RNA_DNA_NC.8,RNA_Protein_NC.8,p2,nrow = 1,ncol = 3,common.legend = F,widths=c(1,1,0.1)),
#                  ggarrange(RNA_DNA_C.4,RNA_Protein_C.4,p1,nrow = 1,ncol = 3,common.legend = T,legend = "bottom",widths=c(1,1,0.1)),
#                  nrow = 8,ncol = 1,heights = c(1.2,1,1,1,1,1,1,1.2))

figure<-ggarrange(ggarrange(RNA_DNA_NC.5,RNA_Protein_NC.5,p2,nrow = 1,ncol = 3,common.legend = T,legend = "top",widths=c(1,1,0.1)),
                  ggarrange(RNA_DNA_C.1,RNA_Protein_C.1,p1,nrow = 1,ncol = 3,common.legend = F,widths=c(1,1,0.1)),
                  ggarrange(RNA_DNA_NC.6,RNA_Protein_NC.6,p2,nrow = 1,ncol = 3,common.legend = F,widths=c(1,1,0.1)),
                  ggarrange(RNA_DNA_C.2,RNA_Protein_C.2,p1,nrow = 1,ncol = 3,common.legend = F,widths=c(1,1,0.1)),
                  ggarrange(RNA_DNA_NC.7,RNA_Protein_NC.7,p2,nrow = 1,ncol = 3,common.legend = F,widths=c(1,1,0.1)),
                  ggarrange(RNA_DNA_C.3,RNA_Protein_C.3,p1,nrow = 1,ncol = 3,common.legend = F,widths=c(1,1,0.1)),
                  ggarrange(RNA_DNA_NC.8,RNA_Protein_NC.8,p2,nrow = 1,ncol = 3,common.legend = F,widths=c(1,1,0.1)),
                  ggarrange(RNA_DNA_C.4,RNA_Protein_C.4,p1,nrow = 1,ncol = 3,common.legend = T,legend = "bottom",widths=c(1,1,0.1)),
                  nrow = 8,ncol = 1,heights = c(1.2,1,1,1,1,1,1,1.2))
     
     
pdf(paste0(names1[i],"_",CPTAC,".Corum-sublocation,ribotop.prnext.mitomiddle.PMbot.pdf"),height=21.5,width=12)
print(figure)
dev.off() 
}




































