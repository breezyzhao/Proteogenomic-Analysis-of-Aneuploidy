setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest")
`%notin%` <- Negate(`%in%`)
names1<-c("COAD","BRCA","OV","HNSC","LUAD","UCEC","ccRCC")
datalist<-list()
for (k in 1:length(names1)){
  cancer<-names1[k]
  CPTAC<-"CPTAC"
  corum <- as.data.frame(read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/database/coreComplexesv3.0.txt", sep = "\t"))
  corum<-corum[corum$Organism=="Human",]
  if(cancer=="OV"){
    data_cnv<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/Ovarian-Prospective-Data/raw_from_python/ovarian_CNV.tsv",
                         sep="\t",header = T)
    data_rna<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/Ovarian-Prospective-Data/raw_from_python/ovarian_transcriptomics.tsv",
                         sep="\t",header = T)
    data_pro<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/Ovarian-Prospective-Data/raw_from_python/ovarian_proteomics.tsv",
                         sep="\t",header = T)
    mutation<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/Ovarian/Ovarian_2021-05-27-mutation.txt",
                         sep="\t",header = T)
  }else if(cancer=="COAD"){
    data_cnv<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/Colon-Prospective-Data/raw_from_python/colon_CNV.tsv",
                         sep="\t",header = T)
    data_rna<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/Colon-Prospective-Data/raw_from_python/colon_transcriptomics.tsv",
                         sep="\t",header = T)
    data_pro<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/Colon-Prospective-Data/raw_from_python/colon_proteomics.tsv",
                         sep="\t",header = T)
    mutation<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/Colon/COADREAD_210528_mutation.txt",
                         sep="\t",header = T)
  }else if(cancer=="BRCA"){
    data_cnv<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/Breast-Prospective-Data/raw_from_python/breast_CNV.tsv",
                         sep="\t",header = T)
    data_rna<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/Breast-Prospective-Data/raw_from_python/breast_transcriptomics.tsv",
                         sep="\t",header = T)
    data_pro<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/Breast-Prospective-Data/raw_from_python/breast_proteomics.tsv",
                         sep="\t",header = T)
    mutation<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/Breast/Breast_2021-05-27-mutation.txt",
                         sep="\t",header = T)
  }else if(cancer=="HNSC"){
    data_cnv<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/data_hnscc/hnscc_CNV_use.txt",
                         sep="\t",header = T)
    data_rna<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/data_hnscc/hnscc_transcriptomics.tsv",
                         sep="\t",header = T)
    data_pro<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/data_hnscc/hnscc_proteomics.tsv",
                         sep="\t",header = T)
    mutation<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/HNSCC/HNSCC_2021-05-27-mutation.txt",
                         sep="\t",header = T)
  }else if(cancer=="UCEC"){
    data_cnv<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/data_endometrial/endometrial_CNV.tsv",
                         sep="\t",header = T)
    data_rna<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/data_endometrial/endometrial_transcriptomics.tsv",
                         sep="\t",header = T)
    data_pro<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/data_endometrial/endometrial_proteomics.tsv",
                         sep="\t",header = T)
    mutation<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/Endomertrial/Endomertrial_2021-05-27-mutation.txt",
                         sep="\t",header = T)
  }else if(cancer =="LUAD"){
    data_cnv<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/data_luad/luad_CNV.tsv",
                         sep="\t",header = T)
    data_rna<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/data_luad/luad_transcriptomics.tsv",
                         sep="\t",header = T)
    data_pro<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/data_luad/luad_proteomics.tsv",
                         sep="\t",header = T)
    mutation<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/LUAD/LUAD_2021-05-27-mutation.txt",
                         sep="\t",header = T)
  }else if(cancer=="ccRCC"){
    data_cnv<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/data_ccrcc/ccrcc_CNV.tsv",
                         sep="\t",header = T)
    data_rna<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/data_ccrcc/ccrcc_transcriptomics.tsv",
                         sep="\t",header = T)
    data_pro<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/data_ccrcc/ccrcc_proteomics.tsv",
                         sep="\t",header = T)
    mutation<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/CCRCC/CCRCC_2021-05-27-mutation.txt",
                         sep="\t",header = T)
  }
  library(ggplot2)
  library(boot)
  if(cancer=="HNSC"| cancer=="OV"|cancer=="COAD"|cancer=="UCEC"|cancer=="LUAD"){
    rownames(data_cnv)<-data_cnv$Patient_ID #BREAST,CCRCC:Name; HNSC,LUAD,Endometrial,OV,Colon:Patient_ID
  }else if(cancer=="BRCA"| cancer=="ccRCC"){
    rownames(data_cnv)<-data_cnv$Name
  }
  data_cnv.use<-as.data.frame(t(data_cnv[,-1]))
  cnv_by_gene2<-data_cnv.use[order(row.names(data_cnv.use)) , order(colnames(data_cnv.use))]
  if(cancer=="ccRCC"| cancer=="BRCA"){
    cnv_by_gene1<-dplyr::select(cnv_by_gene2,-c(Patient_ID,Database_ID))
    cnv_by_gene1<-as.data.frame(sapply(cnv_by_gene1, as.numeric))
    rownames(cnv_by_gene1)<-rownames(cnv_by_gene2)
  }else{
    cnv_by_gene1<-cnv_by_gene2
  }
  cnv_by_gene<-cnv_by_gene1
  rownames(cnv_by_gene)<-rownames(cnv_by_gene2)
  
  
  rownames(data_rna)<-data_rna$Patient_ID
  data_rna.use<-as.data.frame(t(data_rna[,-1]))
  rna_by_gene1<-data_rna.use[order(row.names(data_rna.use)) , order(colnames(data_rna.use))]
  median<-as.data.frame(apply(rna_by_gene1,1,function (x) median(x,na.rm=T)))
  median$group<-ifelse(median$`apply(rna_by_gene1, 1, function(x) median(x, na.rm = T))`>=quantile(median$`apply(rna_by_gene1, 1, function(x) median(x, na.rm = T))`,0.0),1,0)
  rna_by_gene<-rna_by_gene1[rownames(rna_by_gene1)%in% rownames(median[median$group==1,]),]
  
  if(cancer=="HNSC"| cancer=="UCEC"|cancer=="COAD"){
    rownames(data_pro)<-data_pro$Patient_ID#BREAST,OV,CCRCC,LUAD:Name hnsc,Endometrial,Colon:Patient_ID
  }else if(cancer=="BRCA"| cancer=="OV"|cancer=="ccRCC"|cancer=="LUAD"){
    rownames(data_pro)<-data_pro$Name
  }
  data_pro.use<-as.data.frame(t(data_pro[,-1]))
  pro_by_gene1<-data_pro.use[order(row.names(data_pro.use)) , order(colnames(data_pro.use))]
  pro_by_gene1<-pro_by_gene1[,-3:-4]
  pro_by_gene<-as.data.frame(apply(pro_by_gene1,2,as.numeric))
  rownames(pro_by_gene)<-rownames(pro_by_gene1)
  
  if(cancer=="COAD"){
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
  
  
  master_list <- strsplit(as.character(corum$subunits.Entrez.IDs.), split = ";")
  master_list <- unique(as.numeric(as.character(unlist(master_list))))
  # Create master list of CORUM UniProt IDs
  master_list_uniprot <- strsplit(as.character(corum$subunits.UniProt.IDs.), split = ";")
  master_list_uniprot <- unique(as.character(unlist(master_list_uniprot)))
  # Create master list of CORUM gene names (use this)
  master_list_names <- strsplit(as.character(corum$subunits.Gene.name.), split = ";")
  master_list_names <- unique(as.character(unlist(master_list_names)))
  
  cnv.gain<-cnv_by_gene_small
  cnv.gain1<-cnv.gain
  cnv.loss<-cnv_by_gene_small
  cnv.loss1<-cnv.loss
  cnv.variance<-cnv_by_gene_small
  
  cnv.gain[cnv.gain<0.3]<-0
  cnv.gain[cnv.gain>=0.3]<-1
  
  cnv.gain1[cnv.gain1<0.3]<-0
  
  cnv.loss[cnv.loss>(-0.3)]<-0
  cnv.loss[cnv.loss<=(-0.3)]<-1
  
  cnv.loss1[cnv.loss1>(-0.3)]<-0
  
  new.table<-data.frame(gene=rownames(cnv_by_gene_small),gain.freq=rep(NA,nrow(cnv_by_gene_small)),
                        gain.mean=rep(NA,nrow(cnv_by_gene_small)),
                        loss.freq=rep(NA,nrow(cnv_by_gene_small)),loss.mean=rep(NA,nrow(cnv_by_gene_small)),
                        var=rep(NA,nrow(cnv_by_gene_small)))
  
  new.table$gain.freq<-apply(cnv.gain,1,function(x)sum(x)/ncol(cnv.gain))
  new.table$gain.mean<-apply(cnv.gain1,1,function(x)mean(x[x!=0],na.rm=T))
  new.table$loss.freq<-apply(cnv.loss,1,function(x)sum(x)/ncol(cnv.loss))
  new.table$loss.mean<-apply(cnv.loss1,1,function(x)mean(x[x!=0],na.rm=T))
  new.table$var<-apply(cnv.variance,1,var)
  new.table$Cancer<-cancer
  new.table[is.na(new.table)] <- 0
  new.table$Corum<-ifelse(new.table$gene%in%master_list_names,"CORUM","NoCORUM")
  datalist[[k]]<-new.table
  write.table(new.table,paste0(cancer,".gain.loss.var.table.txt"),sep="\t",row.names = F,quote = F)
}
#data.use<-Reduce(function(x, y) merge(x, y, all=TRUE), datalist)
data.use<-do.call(rbind,datalist)
plot.data<-na.omit(data.use)

plot.data$Corum<-factor(plot.data$Corum,levels = c("NoCORUM","CORUM"))

col1<-"#aaa900" # buddha gold#"#FF0033" ###red
col2<-"#6500aa" # purple #"#3300FF" ## blue

library(ggpubr)
library(rstatix)
stat.test1 <- plot.data %>%
  group_by(Cancer) %>%
  wilcox_test(var ~ Corum) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test1 <- stat.test1 %>%
  add_xy_position(x = "Cancer", dodge = 0.8)
stat.test1$y.position<-0.3


n_fun <- function(x){
  return(data.frame(y = 0.25,
                    label = paste0("N=",length(x))))
}
g1<-ggplot(plot.data, aes(x=Cancer, y=var)) + 
  #geom_point(aes(fill = Corum, color=Corum), alpha=0.5,
  #           size = 2, shape = 21, position = position_jitterdodge()) +
  geom_boxplot(aes(fill = Corum),alpha = 0.5,outlier.colour = NA)+
  stat_summary(fun.data = n_fun, geom = "text", 
               aes(group=Corum),
               hjust = 0.5, position = position_dodge(0.8),size=5) +
  scale_color_manual(values=c(col1,col2))+
  scale_fill_manual(values=c(col1,col2)) +
  theme_classic2()+
  ggtitle(paste0("variance"))+
  ylab("variance")+xlab("")+scale_y_continuous(limits = c(0, 0.3))+
  theme(axis.title.y = element_text( size = 15, face = "plain"),
        axis.text.x =element_text( size = 12, face = "plain",angle = 0,hjust = 0.5),
        plot.title = element_text(size=15))+
  theme(legend.position = "right")+stat_pvalue_manual(
    stat.test1,  label = "p",tip.length = 0,size=5
  )


pdf("var.211013.pdf",onefile=FALSE,height=6,width=16)
g1
dev.off() 
  
#####pan-cancer
library(data.table)
coad<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/COAD.gain.loss.var.table.txt",
                 sep="\t",header = T)
brca<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/BRCA.gain.loss.var.table.txt",
                 sep="\t",header = T)
ov<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/OV.gain.loss.var.table.txt",
                 sep="\t",header = T)
hnsc<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/HNSC.gain.loss.var.table.txt",
                 sep="\t",header = T)
luad<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/LUAD.gain.loss.var.table.txt",
                 sep="\t",header = T)
ccrcc<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/ccRCC.gain.loss.var.table.txt",
                 sep="\t",header = T)
ucec<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/UCEC.gain.loss.var.table.txt",
                  sep="\t",header = T)

datalist<-list(coad,brca,ov,hnsc,luad,ccrcc,ucec)
data.use<-Reduce(function(x, y) merge(x, y,by="gene",all=T), datalist)
gene<-data.use$gene
gain.freq<-apply(data.use[,colnames(data.use)%like%"gain.freq"],1,function(x)mean(x,na.rm=T))
gain.mean<-apply(data.use[,colnames(data.use)%like%"gain.mean"],1,function(x)mean(x,na.rm=T))
loss.freq<-apply(data.use[,colnames(data.use)%like%"loss.freq"],1,function(x)mean(x,na.rm=T))
loss.mean<-apply(data.use[,colnames(data.use)%like%"loss.mean"],1,function(x)mean(x,na.rm=T))
var<-apply(data.use[,colnames(data.use)%like%"var"],1,function(x)mean(x,na.rm=T))

pan<-as.data.frame(cbind(gene,gain.freq,gain.mean,loss.freq,loss.mean,var))

pan$Cancer<-"Pan-cancer"
pan$Corum<-ifelse(pan$gene%in%master_list_names,"CORUM","NoCORUM")

write.table(pan,"Pan.gain.loss.var.table.txt",sep="\t",row.names = F,quote = F)
  
  
  
  
  