setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS3/figS3.1/211020/")
library(stringr)
colon<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig3/210813.fisher.pathway/Colon.fisher.txt",sep="\t")
colon$ID<-str_to_title(colon$ID) 
breast<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig3/210813.fisher.pathway/BREAST.fisher.txt",sep="\t")
breast$ID<-str_to_title(breast$ID) 
ovarian<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig3/210813.fisher.pathway/OV.fisher.txt",sep="\t")
ovarian$ID<-str_to_title(ovarian$ID) 
ccRCC<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig3/210813.fisher.pathway/ccRCC.fisher.txt",sep = "\t")
ccRCC$ID<-str_to_title(ccRCC$ID) 
endometrial<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig3/210813.fisher.pathway/Endometrial.fisher.txt",sep = "\t")
endometrial$ID<-str_to_title(endometrial$ID) 
luad<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig3/210813.fisher.pathway/LUAD.fisher.txt",sep = "\t")
luad$ID<-str_to_title(luad$ID) 
hnsc<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig3/210813.fisher.pathway/HNSC.fisher.txt",sep = "\t")
hnsc$ID<-str_to_title(hnsc$ID) 
pancancer<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig3/210813.fisher.pathway/Pan-cancer.fisher.txt",
                      sep="\t")
pancancer$ID<-str_to_title(pancancer$ID) 
names<-list(colon,breast,ovarian,ccRCC,endometrial,luad,hnsc,pancancer)
names1<-c("COAD","BRCA","OV","ccRCC","UCEC","LUAD","HNSC","Pan-cancer")
datalist<-list()

col1<-"#20854e"
col2<-"#104327"
col3<-"#4d9d71"
col4<-"#90c2a7"


for (j in 1:8){
  #j<-3
  data1<-names[[j]]
  rownames(data1)<-data1$ID
  data2<-as.data.frame(data1[,c(9)] )#fdr value
  data.temp<--1*log10(data2)
  #data.temp[data.temp>4]<-4
  data.temp$rank<-rank(data.temp[,1])
  data.temp$ID<-data1$ID
  colnames(data.temp)[1]<-names1[j]
  colnames(data.temp)[2]<-paste0(names1[j],".rank")
  datalist[[j]]<-data.temp
}
HL <- Reduce(
  function(x, y, ...) merge(x, y, all=T,...), 
  datalist
)
HL[is.na(HL)] <- 0

HL$rank.mean<-apply(HL[,c(3,5,7,9,11,13,15,17)],1,mean)
HL<-HL[order(HL$rank.mean,decreasing = T),]

HL1<-HL[1:10,]
HL1<-HL1[order(HL1$rank.mean,decreasing = F),]
datalist1<-list()
for (j in 1:8){
HL1.temp<-HL1[,c(1,j*2)]
HL1.temp$group<-colnames(HL1.temp)[2]
colnames(HL1.temp)[2]<-"logp"
datalist1[[j]]<-HL1.temp
}

data.use<-do.call(rbind,datalist1)
data.use$group<-factor(data.use$group,levels = c("Pan-cancer","COAD","BRCA","OV","ccRCC","UCEC","LUAD","HNSC"))
data.use$log10FDR<-ifelse(data.use$logp>=8,8,data.use$logp)
data.use$ID<-factor(data.use$ID,level=unique(data.use$ID))

g1.bubble<-ggplot(data.use, aes(x =group , y = ID)) + 
  geom_point(aes(size = log10FDR),color=col2,shape=16) +
  #scale_size(range = c(0.5, 8))+# Adjust the range of points size
  theme_classic()+
  scale_x_discrete(name ="")+ylab("")+scale_y_discrete(position = "right")+
  theme(text = element_text(size=12),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"))+
  ggtitle("Top10 enrichments in \n Group1 High DNA-RNA Corr & Low RNA-Pro Corr")
#### FOR LH
datalist<-list()
for (j in 1:8){
  #j<-3
  data1<-names[[j]]
  rownames(data1)<-data1$ID
  data2<-as.data.frame(data1[,c(13)] )#FDR value
  data.temp<--1*log10(data2)
  #data.temp[data.temp>4]<-4
  data.temp$rank<-rank(data.temp[,1])
  data.temp$ID<-data1$ID
  colnames(data.temp)[1]<-names1[j]
  colnames(data.temp)[2]<-paste0(names1[j],".rank")
  datalist[[j]]<-data.temp
}
LH <- Reduce(
  function(x, y, ...) merge(x, y, all=T,...), 
  datalist
)
LH[is.na(LH)] <- 0

LH$rank.mean<-apply(LH[,c(3,5,7,9,11,13,15,17)],1,mean)
LH<-LH[order(LH$rank.mean,decreasing = T),]

LH1<-LH[1:10,]
LH1<-LH1[order(LH1$rank.mean,decreasing = F),]
datalist1<-list()
for (j in 1:8){
  LH1.temp<-LH1[,c(1,j*2)]
  LH1.temp$group<-colnames(LH1.temp)[2]
  colnames(LH1.temp)[2]<-"logp"
  datalist1[[j]]<-LH1.temp
}

data.use2<-do.call(rbind,datalist1)
data.use2$group<-factor(data.use2$group,levels = c("Pan-cancer","COAD","BRCA","OV","ccRCC","UCEC","LUAD","HNSC"))
data.use2$log10FDR<-ifelse(data.use2$logp>=8,8,data.use2$logp)
data.use2$ID<-factor(data.use2$ID,level=unique(data.use2$ID))


g2.bubble<-ggplot(data.use2, aes(x =group , y = ID)) + 
  geom_point(aes(size = log10FDR),color=col4,shape=16) +
  #scale_size(range = c(0.5, 8))+# Adjust the range of points size
  theme_classic()+
  scale_x_discrete(name ="")+ylab("")+scale_y_discrete(position = "right")+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45,size = 12,hjust = 1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"))+
  ggtitle("Top10 enrichments in \n Group2:Low DNA-RNA Corr & High RNA-Pro Corr")

#### FOR HH
datalist<-list()
for (j in 1:8){
  #j<-3
  data1<-names[[j]]
  rownames(data1)<-data1$ID
  data2<-as.data.frame(data1[,c(5)] )#FDR value
  data.temp<--1*log10(data2)
  #data.temp[data.temp>4]<-4
  data.temp$rank<-rank(data.temp[,1])
  data.temp$ID<-data1$ID
  colnames(data.temp)[1]<-names1[j]
  colnames(data.temp)[2]<-paste0(names1[j],".rank")
  datalist[[j]]<-data.temp
}
HH <- Reduce(
  function(x, y, ...) merge(x, y, all=T,...), 
  datalist
)
HH[is.na(HH)] <- 0

HH$rank.mean<-apply(HH[,c(3,5,7,9,11,13,15,17)],1,mean)
HH<-HH[order(HH$rank.mean,decreasing = T),]

HH1<-HH[1:10,]
HH1<-HH1[order(HH1$rank.mean,decreasing = F),]
datalist1<-list()
for (j in 1:8){
  HH1.temp<-HH1[,c(1,j*2)]
  HH1.temp$group<-colnames(HH1.temp)[2]
  colnames(HH1.temp)[2]<-"logp"
  datalist1[[j]]<-HH1.temp
}

data.use3<-do.call(rbind,datalist1)
data.use3$group<-factor(data.use3$group,levels = c("Pan-cancer","COAD","BRCA","OV","ccRCC","UCEC","LUAD","HNSC"))
data.use3$log10FDR<-ifelse(data.use3$logp>=8,8,data.use3$logp)
data.use3$ID<-factor(data.use3$ID,level=unique(data.use3$ID))


g3.bubble<-ggplot(data.use3, aes(x =group , y = ID)) + 
  geom_point(aes(size = log10FDR),color=col1,shape=16) +
  scale_size(range = c(0.5, 8))+# Adjust the range of points size
  theme_classic()+
  scale_x_discrete(name ="")+ylab("")+scale_y_discrete(position = "right")+
  theme(text = element_text(size=12),
        axis.text.x=element_blank(),legend.position = "none",
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"))+
  ggtitle("Top10 enrichments in \n Group 3 High DNA-RNA Corr & Low RNA-Pro Corr")

#### FOR LL
datalist<-list()
for (j in 1:8){
  #j<-3
  data1<-names[[j]]
  rownames(data1)<-data1$ID
  data2<-as.data.frame(data1[,c(17)] )#FDR value
  data.temp<--1*log10(data2)
  #data.temp[data.temp>4]<-4
  data.temp$rank<-rank(data.temp[,1])
  data.temp$ID<-data1$ID
  colnames(data.temp)[1]<-names1[j]
  colnames(data.temp)[2]<-paste0(names1[j],".rank")
  datalist[[j]]<-data.temp
}
LL <- Reduce(
  function(x, y, ...) merge(x, y, all=T,...), 
  datalist
)
LL[is.na(LL)] <- 0

LL$rank.mean<-apply(LL[,c(3,5,7,9,11,13,15,17)],1,mean)
LL<-LL[order(LL$rank.mean,decreasing = T),]

LL1<-LL[1:10,]
LL1<-LL1[order(LL1$rank.mean,decreasing = F),]
datalist1<-list()
for (j in 1:8){
  LL1.temp<-LL1[,c(1,j*2)]
  LL1.temp$group<-colnames(LL1.temp)[2]
  colnames(LL1.temp)[2]<-"logp"
  datalist1[[j]]<-LL1.temp
}

data.use4<-do.call(rbind,datalist1)
data.use4$group<-factor(data.use4$group,levels = c("Pan-cancer","COAD","BRCA","OV","ccRCC","UCEC","LUAD","HNSC"))
data.use4$log10FDR<-ifelse(data.use4$logp>=8,8,data.use4$logp)
data.use4$ID<-factor(data.use4$ID,level=unique(data.use4$ID))


g4.bubble<-ggplot(data.use4, aes(x =group , y = ID)) + 
  geom_point(aes(size = log10FDR),color=col3,shape=16) +
  scale_size(range = c(0.5, 8))+# Adjust the range of points size
  theme_classic()+
  scale_x_discrete(name ="")+ylab("")+scale_y_discrete(position = "right")+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45,size = 12,hjust = 1),legend.position = "none",
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"))+
  ggtitle("Top10 enrichments in \n Group 4 Low DNA-RNA Corr & Low RNA-Pro Corr")


library(ggpubr)
g1<-ggarrange(g1.bubble,g2.bubble,ncol=1,nrow=2,common.legend = F,legend = "none",align = "v",heights = c(1,1.2))
g2<-ggarrange(g3.bubble,g4.bubble,ncol=1,nrow=2,common.legend = F,legend = "none",align = "v",heights = c(1,1.2))
pdf("s3-1.CPTAC.pdf",width = 8,height = 8)
print(g1)
dev.off()
pdf("s3-2.CPTAC.pdf",width = 8,height = 8)
print(g2)
dev.off()








  
