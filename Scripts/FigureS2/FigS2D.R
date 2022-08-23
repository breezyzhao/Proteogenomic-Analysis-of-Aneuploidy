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
col1<-"#aaa900" # buddha gold#"#FF0033" ###red
col2<-"#6500aa" # purple #"#3300FF" ## blue
set.seed(1234)
cancer<-"Multiple_tissue"
CPTAC<-"Normal"


setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210830.normaldatabase")

data_protein<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/dataset_outside/normal_protein_rna/protein.txt",
                     sep="\t",header = T)

data_rna<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/dataset_outside/normal_protein_rna/rna.txt",
                         sep="\t",header = T)


#tissue<-c("Colon","Endometrium","Kidney","Liver","Ovary")
rna_by_gene<-log2(data_rna[,3:21] +1)
rownames(rna_by_gene)<-make.names(data_rna$Gene.name,unique=T)
#rna_by_gene<-rna_by_gene[,colnames(rna_by_gene)%in%tissue,]

data_protein[data_protein == -Inf] <- 0
data_protein<-na.omit(data_protein)
pro_by_gene<-na.omit(log2(data_protein[,3:21] +1))
rownames(pro_by_gene)<-make.names(data_protein$Gene.name,unique = T)

datanames1<-intersect(rownames(pro_by_gene),rownames(rna_by_gene))
datanames2<-intersect(colnames(pro_by_gene),colnames(rna_by_gene))

rna_by_gene_small<-rna_by_gene[rownames(rna_by_gene) %in% datanames1, colnames(rna_by_gene) %in% datanames2]
rna_by_gene_small<-rna_by_gene_small[order(rownames(rna_by_gene_small)),]
pro_by_gene_small<-pro_by_gene[rownames(pro_by_gene) %in% datanames1,colnames(pro_by_gene) %in% datanames2]
pro_by_gene_small<-pro_by_gene_small[order(rownames(pro_by_gene_small)),]

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
#p.boot1<-p.adjust(signif(bootstrap(cor_rna_prot,10000),3),n=3,method = "fdr")
p.boot1<-p.adjust(signif(bootstrap(cor_rna_prot,10000,"corr","corum",TRUE,FALSE),3),n=3,method = "fdr")
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
                        legend.position = "None")+
  scale_fill_manual(values=c(col1, col2))+
  scale_color_manual(values=c(col1, col2))+
  geom_hline(yintercept=0, colour="white", size=2)

library(gridExtra)
library(ggpubr)
pdf(paste0(cancer,"_",CPTAC,"_DNA_RNA_PRO.pdf"),onefile=FALSE,height=3,width=3)
print(RNA_Protein)
dev.off() 

  
  
  
  
  
  
  