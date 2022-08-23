bootstrap<-function(dataset,iterations){
  s.size.g1 <- length(dataset$corr[dataset$corum=="TRUE"])
  s.size.g2 <- length(dataset$corr[dataset$corum=="FALSE"])
  pool <- dataset$corr
  obs.diff.b1 <- median (dataset$corr[dataset$corum=="TRUE"]) - median (dataset$corr[dataset$corum=="FALSE"])
  iterations <- iterations
  sampl.dist.b1 <- NULL
  for (i in 1 : iterations) {
    resample <- sample (c(1:length (pool)), length(pool), replace = TRUE) 
    # "replace = TRUE" is the only difference between bootstrap and permutations
    
    g1.perm = pool[resample][1 : s.size.g1]
    g2.perm = pool[resample][(s.size.g1+1) : length(pool)]
    sampl.dist.b1[i] = median (g1.perm) - median (g2.perm) 
  }
  p.boot1 <- (sum ( abs(sampl.dist.b1) >= abs(obs.diff.b1)) + 1)/ (iterations+1)
  return(p.boot1)
}

set.seed(1234)
setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/211029.rmEXP10pct")
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
col1<-"#aaa900" # buddha gold#"#FF0033" ###red
col2<-"#6500aa" # purple #"#3300FF" ## blue
names1<-c("Colon","BREAST","OV","HNSC","LUAD","Endometrial","ccRCC")
for (k in 1:length(names1)){
cancer<-names1[k]
CPTAC<-"CPTAC"
corum <- as.data.frame(read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/database/coreComplexesv3.0.txt", sep = "\t"))
corum<-corum[corum$Organism=="Human",]
if(cancer=="OV"){
data_cnv<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/Ovarian-Prospective-Data/raw_from_python/ovarian_CNV.tsv",
                     sep="\t",header = T)
data_rna<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/Ovarian-Prospective-Data/raw_from_python/ovarian_transcriptomics.tsv",
                     sep="\t",header = T)
data_pro<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/Ovarian-Prospective-Data/raw_from_python/ovarian_proteomics.tsv",
                      sep="\t",header = T)
mutation<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/Ovarian/Ovarian_2021-05-27-mutation.txt",
                     sep="\t",header = T)
}else if(cancer=="Colon"){
  data_cnv<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/Colon-Prospective-Data/raw_from_python/colon_CNV.tsv",
                       sep="\t",header = T)
  data_rna<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/Colon-Prospective-Data/raw_from_python/colon_transcriptomics.tsv",
                       sep="\t",header = T)
  data_pro<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/Colon-Prospective-Data/raw_from_python/colon_proteomics.tsv",
                       sep="\t",header = T)
  mutation<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/Colon/COADREAD_210528_mutation.txt",
                       sep="\t",header = T)
}else if(cancer=="BREAST"){
  data_cnv<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/Breast-Prospective-Data/raw_from_python/breast_CNV.tsv",
                       sep="\t",header = T)
  data_rna<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/Breast-Prospective-Data/raw_from_python/breast_transcriptomics.tsv",
                       sep="\t",header = T)
  data_pro<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/Breast-Prospective-Data/raw_from_python/breast_proteomics.tsv",
                       sep="\t",header = T)
  mutation<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/Breast/Breast_2021-05-27-mutation.txt",
                       sep="\t",header = T)
}else if(cancer=="HNSC"){
  data_cnv<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/data_hnscc/hnscc_CNV_use.txt",
                       sep="\t",header = T)
  data_rna<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/data_hnscc/hnscc_transcriptomics.tsv",
                       sep="\t",header = T)
  data_pro<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/data_hnscc/hnscc_proteomics.tsv",
                       sep="\t",header = T)
  mutation<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/HNSCC/HNSCC_2021-05-27-mutation.txt",
                       sep="\t",header = T)
}else if(cancer=="Endometrial"){
  data_cnv<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/data_endometrial/endometrial_CNV.tsv",
                       sep="\t",header = T)
  data_rna<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/data_endometrial/endometrial_transcriptomics.tsv",
                       sep="\t",header = T)
  data_pro<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/data_endometrial/endometrial_proteomics.tsv",
                       sep="\t",header = T)
  mutation<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/Endomertrial/Endomertrial_2021-05-27-mutation.txt",
                       sep="\t",header = T)
}else if(cancer =="LUAD"){
  data_cnv<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/data_luad/luad_CNV.tsv",
                       sep="\t",header = T)
  data_rna<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/data_luad/luad_transcriptomics.tsv",
                       sep="\t",header = T)
  data_pro<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/data_luad/luad_proteomics.tsv",
                       sep="\t",header = T)
  mutation<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/LUAD/LUAD_2021-05-27-mutation.txt",
                       sep="\t",header = T)
}else if(cancer=="ccRCC"){
  data_cnv<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/data_ccrcc/ccrcc_CNV.tsv",
                       sep="\t",header = T)
  data_rna<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/data_ccrcc/ccrcc_transcriptomics.tsv",
                       sep="\t",header = T)
  data_pro<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/data_ccrcc/ccrcc_proteomics.tsv",
                       sep="\t",header = T)
  mutation<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/10.protein/CPTAC/CPTAC2_3/2020June_mutation/output/CCRCC/CCRCC_2021-05-27-mutation.txt",
                       sep="\t",header = T)
}
library(ggplot2)
library(boot)
if(cancer=="HNSC"| cancer=="OV"|cancer=="Colon"|cancer=="Endometrial"|cancer=="LUAD"){
  rownames(data_cnv)<-data_cnv$Patient_ID #BREAST,CCRCC:Name; HNSC,LUAD,Endometrial,OV,Colon:Patient_ID
}else if(cancer=="BREAST"| cancer=="ccRCC"){
  rownames(data_cnv)<-data_cnv$Name
}
data_cnv.use<-as.data.frame(t(data_cnv[,-1]))
cnv_by_gene2<-data_cnv.use[order(row.names(data_cnv.use)) , order(colnames(data_cnv.use))]
if(cancer=="ccRCC"| cancer=="BREAST"){
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
median$group<-ifelse(median$`apply(rna_by_gene1, 1, function(x) median(x, na.rm = T))`>=quantile(median$`apply(rna_by_gene1, 1, function(x) median(x, na.rm = T))`,0.1),1,0)
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
mutation.use<-mutation[mutation$Total_Num_Mutations>=0,]
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

# colon_rna_dna -----------------------------------------------------------
#test1 <- apply(cnv_by_gene_colon[,2:106],1,function(x) subset(which(x>0.02)))

all(colnames(rna_by_gene_small) == colnames(cnv_by_gene_small)) # should be true, might be false
all(row.names(rna_by_gene_small) == row.names(cnv_by_gene_small)) # should be true, might be false

# Make correlation df
cor_rna_dna <- data.frame(gene=row.names(rna_by_gene_small), corr=rep(NA, dim(rna_by_gene_small)[1]), 
                          corum=rep(NA, dim(rna_by_gene_small)[1]),zero_percent=rep(NA, dim(rna_by_gene_small)[1]))
for(i in 1:dim(cor_rna_dna)[1]) {
  if(cor_rna_dna$gene[i] %in% master_list_names) {cor_rna_dna$corum[i]<-TRUE}
}
# If not in CORUM, then it's FALSE
cor_rna_dna$corum[is.na(cor_rna_dna$corum)] = FALSE

# Calculate correlation betwen RNA and DNA #########8.29 ADD CUTOFF OF 0.02 calculate the sample percentage of -0.02~0.02 for each gene
for(i in 1:dim(cor_rna_dna)[1]) {
  #i<-980
  rna<-as.data.frame(as.numeric(rna_by_gene_small[i,]))
  dna<-as.data.frame(as.numeric(cnv_by_gene_small[i,]))
  data_com<-cbind(rna,dna)
  colnames(data_com)<-c("rna","dna")
  data_com<-na.omit(data_com)
  data_com$dna[data_com$dna> -0.02 & 0.02 > data_com$dna]<-0
  if(nrow(data_com)>0){
  for (n in 1:nrow(data_com)){
    if (data_com$dna[n] ==0) {data_com$group[n]<-0}
    else{data_com$group[n]<-1}
  }
  data_used<-data_com[data_com$group==1,]
  }
  if(nrow(data_used)==0){
    cor_rna_dna$corr[i]<-0
    cor_rna_dna$zero_percent[i]<-100
  }else{
    cor_rna_dna$corr[i] <- cor(data_used$rna, data_used$dna, method = "spearman", use = "pairwise.complete.obs")
    cor_rna_dna$zero_percent[i]<-sum(data_com$dna>-0.02 & data_com$dna<0.02)/nrow(data_com)
  }
}

# Remove rows with NA correlation values, this would result from having no values for a particular gene in any sample
cor_rna_dna <- cor_rna_dna[complete.cases(cor_rna_dna),]
cor_rna_dna_new<-cor_rna_dna[cor_rna_dna$zero_percent<=0.7,]#####remove >50% patients has a cnv between -0.02~0.02
cor_rna_dna_colon<-cor_rna_dna_new
# Limit analysis to relevant genes, if you want
corums_cors <- cor_rna_dna_new[cor_rna_dna_new$corum == TRUE, "corr"]
noncorums_cors <- cor_rna_dna_new[cor_rna_dna_new$corum == FALSE, "corr"]
p.boot1<-p.adjust(signif(bootstrap(cor_rna_dna_new,10000),3),n=3,method = "fdr")
# g<-wilcox.test(corums_cors, noncorums_cors)
# pval_dna_rna<-signif(g$p.value,3)
# Calculate p-value with Mann-Whitney test
# Calculate effect size of the distribution shift
shift_dna_rna<-round(((median(corums_cors)-median(noncorums_cors))/median(noncorums_cors))*100,3)
#write.table(cor_rna_dna_new,"Prospective_rna_dna_cor.txt",sep='\t',row.names = F)
# Plot density
RNA_DNA<-ggplot(cor_rna_dna_new, aes(corr, stat(density),color = corum)) + geom_density(alpha=0.1,size = 2) + xlab("RNA correlation with DNA") + 
  theme(text = element_text(size = 20)) + geom_vline(data=cor_rna_dna_new[cor_rna_dna_new$corum==TRUE,], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
  geom_vline(data = cor_rna_dna_new[cor_rna_dna_new$corum==FALSE,], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
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
########################################################################################################
all(colnames(rna_by_gene_small) == colnames(pro_by_gene_small))
all(row.names(rna_by_gene_small) == row.names(pro_by_gene_small))

# Make correlation df
cor_rna_prot <- data.frame(gene=row.names(rna_by_gene_small), corr=rep(NA, dim(rna_by_gene_small)[1]), 
                           corum=rep(NA, dim(rna_by_gene_small)[1]))
for(i in 1:dim(cor_rna_prot)[1]){
  if(row.names(rna_by_gene_small)[i] %in% master_list_names) {cor_rna_prot$corum[i]<-TRUE}
  
}
# If not in CORUM, then it's FALSE
cor_rna_prot$corum[is.na(cor_rna_prot$corum)] = FALSE

# Calculate correlation between RNA and protein ####
for(i in 1:dim(cor_rna_prot)[1]) {
  cor_rna_prot$corr[i] <- cor(as.numeric(rna_by_gene_small[i,]), as.numeric(pro_by_gene_small[i,]), method = "spearman", use = "pairwise.complete.obs")
}

# Remove rows with NA correlation values, this would result from having no values for a particular gene in any sample
cor_rna_prot <- cor_rna_prot[complete.cases(cor_rna_prot),]
cor_rna_prot_colon<-cor_rna_prot
# Limit analysis to relevant genes, if you want
corums_cors <- cor_rna_prot[cor_rna_prot$corum == TRUE, "corr"]
noncorums_cors <- cor_rna_prot[cor_rna_prot$corum == FALSE, "corr"]
p.boot1<-p.adjust(signif(bootstrap(cor_rna_prot,10000),3),n=3,method = "fdr")
# g<-wilcox.test(corums_cors, noncorums_cors)
# pval_rna_prot<-signif(g$p.value,3)
# Calculate effect size of the distribution shift
#shift_rna_prot<-signif(((median(corums_cors)-median(corums_cors))/median(corums_cors))*100,3)
#write.table(cor_rna_dna,"Prospective_rna_prot_cor.txt",sep='\t',row.names = F)
# Plot density
RNA_Protein<-ggplot(cor_rna_prot, aes(corr, stat(density),color = corum)) + geom_density(alpha=0.1,size = 2) + xlab("RNA correlation with Protein") + 
  theme(text = element_text(size = 20)) + geom_vline(data=cor_rna_prot[cor_rna_prot$corum==TRUE,], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
  geom_vline(data = cor_rna_prot[cor_rna_prot$corum==FALSE,], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
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
#################################################################################################################3
all(colnames(pro_by_gene_small) == colnames(cnv_by_gene_small)) # should be true, might be false
all(row.names(pro_by_gene_small) == row.names(cnv_by_gene_small)) # should be true, might be false
# Order the dfs, then try the logical checks again
# Make correlation df
cor_dna_prot <- data.frame(gene=row.names(pro_by_gene_small), corr=rep(NA, dim(pro_by_gene_small)[1]), 
                           corum=rep(NA, dim(pro_by_gene_small)[1]),zero_percent=rep(NA, dim(pro_by_gene_small)[1]))
for(i in 1:dim(cor_dna_prot)[1]) {
  if(cor_dna_prot$gene[i] %in% master_list_names) {cor_dna_prot$corum[i]<-TRUE}
}
# If not in CORUM, then it's FALSE
cor_dna_prot$corum[is.na(cor_dna_prot$corum)] = FALSE
# Calculate correlation between DNA and protein ##### add cutoff 0.02
for(i in 1:dim(cor_dna_prot)[1]) {
  #i<-1
  prot<-as.data.frame(as.numeric(pro_by_gene_small[i,]))
  dna<-as.data.frame(as.numeric(cnv_by_gene_small[i,]))
  data_com<-cbind(prot,dna)
  colnames(data_com)<-c("prot","dna")
  data_com<-na.omit(data_com)
  data_com$dna[data_com$dna> -0.02 & 0.02 > data_com$dna]<-0
  if(nrow(data_com)>0){
    for (n in 1:nrow(data_com)){
      if (data_com$dna[n] ==0) {data_com$group[n]<-0}
      else{data_com$group[n]<-1}
    }
    data_used<-data_com[data_com$group==1,]
  }
  if(nrow(data_used)==0){
    cor_dna_prot$corr[i]<-0
    cor_dna_prot$zero_percent[i]<-100
  }else{
    cor_dna_prot$corr[i] <- cor(data_used$prot, data_used$dna, method = "spearman", use = "pairwise.complete.obs")
    cor_dna_prot$zero_percent[i]<-sum(data_com$dna>-0.02 & data_com$dna<0.02)/nrow(data_com)
  }
}
# Remove rows with NA correlation values, this would result from having no values for a particular gene in any sample
cor_dna_prot <- cor_dna_prot[complete.cases(cor_dna_prot),]
cor_dna_prot_new<-cor_dna_prot[cor_dna_prot$zero_percent<=0.7,] #####remove >50% patients has a cnv between -0.02~0.02
cor_dna_prot_colon<-cor_dna_prot_new
# Calculate p-value with Mann-Whitney test
corums_cors <- cor_dna_prot_new[cor_dna_prot_new$corum == TRUE, "corr"]
noncorums_cors <- cor_dna_prot_new[cor_dna_prot_new$corum == FALSE, "corr"]
p.boot1<-p.adjust(signif(bootstrap(cor_dna_prot_new,10000),3),n=3,method = "fdr")
# g<-wilcox.test(corums_cors, noncorums_cors)
# pval_dna_prot<-signif(g$p.value,3)
# Calculate effect size of the distribution shift
shift_dna_prot<-round(((median(noncorums_cors)-median(corums_cors))/median(corums_cors))*100,3)
# Plot density
DNA_Protein<-ggplot(cor_dna_prot_new, aes(corr, stat(density),color = corum)) + geom_density(alpha=0.1,size = 2) + xlab("DNA correlation with Protein") + 
  theme(text = element_text(size = 20)) + geom_vline(data=cor_dna_prot_new[cor_dna_prot_new$corum==TRUE,], aes(xintercept = median(corr)), color=col2, size = 1.5, linetype='solid') + 
  geom_vline(data = cor_dna_prot_new[cor_dna_prot_new$corum==FALSE,], aes(xintercept = median(corr)), color = col1, size = 1.5, linetype='solid')+
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
merge.temp<-merge(cor_rna_dna_colon,cor_rna_prot_colon,by="gene") 
merge.data<-merge(merge.temp,cor_dna_prot_colon,by="gene")
colnames(merge.data)[2]<-"corr.rna.dna"
colnames(merge.data)[5]<-"corr.rna.pro"
colnames(merge.data)[7]<-"cor.dna.pro"
write.table(merge.data,paste0(cancer,".corr.table.txt"),sep = "\t",row.names = F,quote = F)

figure<-ggarrange(RNA_DNA,RNA_Protein,DNA_Protein, ncol=3,nrow=1,common.legend = TRUE,legend = "right")
library(gridExtra)
library(ggpubr)
 pdf(paste0(cancer,"_",CPTAC,"_DNA_RNA_PRO.pdf"),onefile=FALSE,height=3,width=10)
 print(figure)
 dev.off() 
}

