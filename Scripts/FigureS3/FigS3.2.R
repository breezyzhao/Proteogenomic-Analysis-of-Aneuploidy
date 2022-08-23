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
library(ggridges)
setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210827.sublocation.all/")

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
nucleus.loc<-c("Nucleoplasm","Nuclear speckles","Nuclear bodies","Nuclear membrane")
cytoplasm.loc<-c("Microtubules","Cytosol","Actin filaments","Centrosome","Centriolar satellite","Cytoplasmic bodies","Intermediate filaments",
                 "Cytokinetic bridge","Mitotic spindle","Microtubule ends")
nucleoli.loc<-c("Nucleoli","Nucleoli fibrillar center","Nucleoli rim")
organelles.loc<-c("Mitochondria","Vesicles","Golgi apparatus","Peroxisomes",
                  "Endosomes","Lysosomes")
er.loc<-c("Endoplasmic reticulum")
plasma_membrane.loc<-c("Plasma membrane")
mito.loc<-c("Mitochondria")
ves.loc<-c("Vesicles")
golgi.loc<-c("Golgi apparatus")
pero.loc<-c("Peroxisomes")
endo.loc<-c("Endosomes")
lyso.loc<-c("Lysosomes")

ribo.list<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210824.ribo/ribosomal.txt",
                      sep="\t",header = T)
proteasome.list<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210824.ribo/proteasome.txt",
                            sep="\t",header = T)
nucleus.list<-as.data.frame(unique(subloc[subloc$Main.location %in% nucleus.loc,]))
cytoplasm.list<-as.data.frame(unique(subloc[subloc$Main.location %in% cytoplasm.loc,]))
nucleoli.list<-as.data.frame(unique(subloc[subloc$Main.location %in% nucleoli.loc,]))
organelles.list<-as.data.frame(unique(subloc[subloc$Main.location %in% organelles.loc,]))
plasma_membrane.list<-as.data.frame(unique(subloc[subloc$Main.location %in% plasma_membrane.loc,]))
er.list<-as.data.frame(unique(subloc[subloc$Main.location %in% er.loc,]))
mito.list<-as.data.frame(unique(subloc[subloc$Main.location %in% mito.loc,]))
ves.list<-as.data.frame(unique(subloc[subloc$Main.location %in% ves.loc,]))
golgi.list<-as.data.frame(unique(subloc[subloc$Main.location %in% golgi.loc,]))
pero.list<-as.data.frame(unique(subloc[subloc$Main.location %in% pero.loc,]))
endo.list<-as.data.frame(unique(subloc[subloc$Main.location %in% endo.loc,]))
lyso.list<-as.data.frame(unique(subloc[subloc$Main.location %in% lyso.loc,]))
er.ribo.list<-as.data.frame(unique(subloc[subloc$Main.location %in% er.loc & subloc $Gene.name %in%ribo.list$Gene.name,]))
er.other.list<-as.data.frame(unique(subloc[subloc$Main.location %in% er.loc & subloc $Gene.name %notin%er.ribo.list$Gene.name,]))



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
score$Corum<-ifelse(score$genes%in%master_list_names,"Corum","NoCorum")


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
# col1<-"#AAA900" ###budddha gold
# col2<-"#6500aa" ## purple
# col3<-"#bfc2c3" 
# col4<-"#acafb0" 
# col5<-"#999b9c" 
# col6<-"#868889" 
# col8.nc<-"#777600" ### er
# col7.nc<-"#c4c34d" ### pm
# col8<-"#3d0066" ## er
# col7<-"#a366cc" ###pm
# col9<-"#d7dbdb"
# cancer<-"CPTAC"

for (j in 1:8){
  #j<-8
  cancer<-"Subloc"
  CPTAC<-names1[j]
  merge.data<-names[[j]]
  if(names1[j]=="PAN-Cancer"){
    colnames(merge.data)[1]<-"genes" 
    score.data<-merge(score,merge.data[,c("genes","mean.rna.dna","mean.rna.pro")],by="genes")
    colnames(score.data)[5]<-"corr.rna.dna"
    colnames(score.data)[6]<-"corr.rna.pro"
  }else{
    colnames(merge.data)[1]<-"genes" 
    score.data<-merge(score,merge.data[,c(1,2,5)],by="genes")
    colnames(score.data)[5]<-"corr.rna.dna"
    colnames(score.data)[6]<-"corr.rna.pro"
  }
  corum<-score.data[score.data$Corum=="Corum",]
  corum$group<-"Corum:all"
  corum.ribo<-score.data[score.data$genes%in%ribo.list$Gene.name & score.data$Corum=="Corum",]
  corum.ribo$group<-"Corum:Ribosome"
  corum.PR<-score.data[score.data$genes%in%proteasome.list$Gene.name & score.data$Corum=="Corum",]
  corum.PR$group<-"Corum:Proteasome"
  corum.nuclues<-score.data[score.data$genes%in%nucleus.list$Gene.name & score.data$Corum=="Corum",]
  corum.nuclues$group<-"Corum:nucleus"
  corum.cytoplasm<-score.data[score.data$genes%in%cytoplasm.list$Gene.name & score.data$Corum=="Corum",]
  corum.cytoplasm$group<-"Corum:cytoplasm"
  corum.nucleoli<-score.data[score.data$genes%in%nucleoli.list$Gene.name & score.data$Corum=="Corum",]
  corum.nucleoli$group<-"Corum:nucleoli"
  corum.organelle<-score.data[score.data$genes%in%organelles.list$Gene.name & score.data$Corum=="Corum",]
  corum.organelle$group<-"Corum:organelles"
  corum.PM<-score.data[score.data$genes%in%plasma_membrane.list$Gene.name & score.data$Corum=="Corum",]
  corum.PM$group<-"Corum:plasma_membrane"
  corum.ER<-score.data[score.data$genes%in%er.list$Gene.name & score.data$Corum=="Corum",]
  corum.ER$group<-"Corum:ER"
  corum.ER.ribo<-score.data[score.data$genes%in%er.ribo.list$Gene.name & score.data$Corum=="Corum",]
  corum.ER.ribo$group<-"Corum:ER.Ribo"
  corum.ER.other<-score.data[score.data$genes%in%er.other.list$Gene.name & score.data$Corum=="Corum",]
  corum.ER.other$group<-"Corum:ER.other"
  corum.mito<-score.data[score.data$genes%in%mito.list$Gene.name & score.data$Corum=="Corum",]
  corum.mito$group<-"Corum:Mitochondria"
  corum.ves<-score.data[score.data$genes%in%ves.list$Gene.name & score.data$Corum=="Corum",]
  corum.ves$group<-"Corum:Vesicles"
  corum.golgi<-score.data[score.data$genes%in%golgi.list$Gene.name & score.data$Corum=="Corum",]
  corum.golgi$group<-"Corum:Golgi apparatus"
  corum.pero<-score.data[score.data$genes%in%pero.list$Gene.name & score.data$Corum=="Corum",]
  corum.pero$group<-"Corum:Peroxisomes"
  corum.endo<-score.data[score.data$genes%in%endo.list$Gene.name & score.data$Corum=="Corum",]
  corum.endo$group<-"Corum:Endosomes"
  corum.lyso<-score.data[score.data$genes%in%lyso.list$Gene.name & score.data$Corum=="Corum",]
  corum.lyso$group<-"Corum:Lysosomes"
  
  nocorum<-score.data[score.data$Corum=="NoCorum",]
  nocorum$group<-"NoCorum:all"
  #nocorum.ribo<-score.data[score.data$genes%in%ribo.list$Gene.name & score.data$Corum=="NoCorum",]
  #nocorum.ribo$group<-"NoCorum:Ribosome"
  nocorum.PR<-score.data[score.data$genes%in%proteasome.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.PR$group<-"NoCorum:Proteasome"
  nocorum.nuclues<-score.data[score.data$genes%in%nucleus.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.nuclues$group<-"NoCorum:nucleus"
  nocorum.cytoplasm<-score.data[score.data$genes%in%cytoplasm.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.cytoplasm$group<-"NoCorum:cytoplasm"
  nocorum.nucleoli<-score.data[score.data$genes%in%nucleoli.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.nucleoli$group<-"NoCorum:nucleoli"
  nocorum.organelle<-score.data[score.data$genes%in%organelles.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.organelle$group<-"NoCorum:organelles"
  nocorum.PM<-score.data[score.data$genes%in%plasma_membrane.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.PM$group<-"NoCorum:plasma_membrane"
  nocorum.ER<-score.data[score.data$genes%in%er.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.ER$group<-"NoCorum:ER"
  #nocorum.ER.ribo<-score.data[score.data$genes%in%er.ribo.list$Gene.name & score.data$Corum=="NoCorum",]
  #nocorum.ER.ribo$group<-"NoCorum:ER.Ribo"
  nocorum.ER.other<-score.data[score.data$genes%in%er.other.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.ER.other$group<-"NoCorum:ER.other"
  nocorum.mito<-score.data[score.data$genes%in%mito.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.mito$group<-"NoCorum:Mitochondria"
  nocorum.ves<-score.data[score.data$genes%in%ves.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.ves$group<-"NoCorum:Vesicles"
  nocorum.golgi<-score.data[score.data$genes%in%golgi.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.golgi$group<-"NoCorum:Golgi apparatus"
  nocorum.pero<-score.data[score.data$genes%in%pero.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.pero$group<-"NoCorum:Peroxisomes"
  nocorum.endo<-score.data[score.data$genes%in%endo.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.endo$group<-"NoCorum:Endosomes"
  nocorum.lyso<-score.data[score.data$genes%in%lyso.list$Gene.name & score.data$Corum=="NoCorum",]
  nocorum.lyso$group<-"NoCorum:Lysosomes"
  
  plot.data<-rbind(corum,corum.cytoplasm,corum.endo,corum.ER,corum.ER.other,corum.ER.ribo,
                   corum.golgi,corum.lyso,corum.mito,corum.nucleoli,corum.nuclues,corum.organelle,
                   corum.pero,corum.PM,corum.PR,corum.ribo,corum.ves,
                   nocorum,nocorum.cytoplasm,nocorum.endo,nocorum.ER,nocorum.ER.other,
                   nocorum.golgi,nocorum.lyso,nocorum.mito,nocorum.nucleoli,nocorum.nuclues,nocorum.organelle,
                   nocorum.pero,nocorum.PM,nocorum.PR,nocorum.ves)
  plot.data$group<-factor(plot.data$group,levels = c("Corum:Lysosomes","NoCorum:Lysosomes","Corum:Endosomes","NoCorum:Endosomes",
                                                     "Corum:Peroxisomes","NoCorum:Peroxisomes","Corum:Golgi apparatus","NoCorum:Golgi apparatus",
                                                     "Corum:Vesicles","NoCorum:Vesicles","Corum:Mitochondria","NoCorum:Mitochondria",
                                                     "Corum:ER","NoCorum:ER","Corum:ER.Ribo",
                                                     "Corum:ER.other","NoCorum:ER.other","Corum:Ribosome",
                                                     "Corum:plasma_membrane","NoCorum:plasma_membrane",
                                                     "Corum:organelles","NoCorum:organelles",
                                                     "Corum:nucleoli","NoCorum:nucleoli",
                                                     "Corum:cytoplasm","NoCorum:cytoplasm",
                                                     "Corum:nucleus","NoCorum:nucleus",
                                                     "Corum:Proteasome","NoCorum:Proteasome",
                                                     "Corum:all", "NoCorum:all"))  
  write.table(plot.data,paste0(names1[j],".subloc.corr.txt"),sep = "\t",row.names = F,quote = F)
  g1<-ggplot(plot.data, aes(x = corr.rna.dna, y = group)) +
    geom_density_ridges(scale = 2,quantile_lines=TRUE,
                        quantile_fun=function(x,...)median(x),
                        aes(fill = group)) + 
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
    geom_density_ridges(scale = 2,quantile_lines=TRUE,
                        quantile_fun=function(x,...)median(x),
                        aes(fill = group)) + 
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
  pdf(paste0(names1[j],".sublocation.pdf"),width = 12,height = 16)
  print(figure)
  dev.off()
}
  
  
  
  
  
  
  
  
  
   
     
     














