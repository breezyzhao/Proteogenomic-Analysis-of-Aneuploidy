 setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/211026.control.test.fig3a.b/")
 data<-read.delim("comprehensive.table.211026.txt",
                  sep="\t",header = T)
 rownames(data)<-data$item
# #data1<-data[,c(-1,-11:-16)]
#  data1<-data[,c(-1)]
#  rownames(data1)<-data$item
#  data2<-as.data.frame(combn(colnames(data1),2))
#  datalist<-list()
#  for (i in 1:ncol(data2)){
#    data.temp<-data1[,colnames(data1) %in% data2[,i]]
#    names<-paste(colnames(data.temp),sep="",collapse = "_")
#    corr1<-cor.test(data.temp[,1],data.temp[,2],method = "pearson")
#    corr<-signif(corr1$estimate,3)
#    p.value<-signif(corr1$p.value,3)
#    corr.2<-cor.test(data.temp[,1],data.temp[,2],method = "spearman")
#    corr2<-signif(corr.2$estimate,3)
#    p.value2<-signif(corr.2$p.value,3)
#    data.new<-cbind(names,corr,p.value,corr2,p.value2)
#    colnames(data.new)[2]<-"pearson"
#    colnames(data.new)[4]<-"spearman"
#    datalist[[i]]<-data.new
#  }
#  data.new1<-do.call(rbind,datalist)
#  write.table(data.new1,"211014.corr.pathway.table.txt",sep="\t",row.names = F,quote = F)


#=================== control test
 library(ggcorrplot)
 library(colorspace)
 hcl_palettes(plot = T)
 col_fun = diverging_hcl("Blue-Red2",n=5)
setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/211026.control.test.fig3a.b/fig3b/")
#cancer<-c("Pan","COAD","BRCA","OV","ccRCC","LUAD","UCEC","HNSC","hCEC","hCEC.all")
# cancer<-"Pan"
# for (i in 1:length(cancer)){
# plot.data<-data
# require("ggrepel")
# set.seed(1234)
# g1<-#ggplot(plot.data, aes(x=paste0(cancer[i],".DNA.RNA.Corr"), y=paste0(cancer[i],".RNA.Pro.corr"),colour=mRNA.halflife)) +
#   ggplot(plot.data, aes_string(x=paste0(cancer[i],".DNA.RNA.Corr"), y=paste0(cancer[i],".RNA.Pro.Corr"))) +
#   geom_point(aes(size=protein.halflife,colour=mRNA.halflife),shape = 16) +
#    geom_errorbar(aes(ymin =Pan.RNA.Pro.Corr.SE.min ,ymax = Pan.RNA.Pro.Corr.SE.max),size=0.5,alpha=0.7,colour="darkgrey") + 
#    geom_errorbarh(aes(xmin = Pan.DNA.RNA.Corr.SE.min,xmax = Pan.DNA.RNA.Corr.SE.max),size=0.5,alpha=0.7,colour="darkgrey")+
#   scale_colour_gradient2(col_fun) +
#   theme_classic()+
#   scale_size(range = c(1,15),name = "Pro.halflife")+
#   #geom_text_repel(min.segment.length = 0,aes(label = paste0(rownames(plot.data),"(",round(plot.data$CORUM.pct,3),"_N=",plot.data$N.of.genes.pathway,")")),
# #                  size = 5,color="black") +
#    geom_text_repel(min.segment.length = 0,aes(label = paste0(rownames(plot.data))),
#                    size = 5,color="black") +
#   geom_smooth(method='lm', formula= y~x,colour="black",se = FALSE)+
#   stat_cor(method = "spearman", label.x.npc = "left", label.y.npc ="bottom" ,size = 5)
#    
# pdf(paste0(cancer[i],".corr.spearman.pdf"),width = 14,height=12)
# print(g1)
# dev.off()
# }


######
plot.data<-data
require(ggrepel)
library(ggpubr)
set.seed(1234)
#hcl_palettes(plot = T)
#col_fun = diverging_hcl("Blue-Red2",n=5)
g1<-ggplot(plot.data, aes(x=Pan.DNA.RNA.Corr,y=Pan.RNA.Pro.Corr,colour=mRNA.halflife)) +
   geom_point(aes(size=protein.halflife),shape = 16) +
   theme_classic()+
   scale_colour_gradient2(low = "#4A6FE3", high = "#D33F6A",mid = "#E2E2E2",midpoint = 12) +
   scale_size(range = c(1,15),name = "Pro.halflife")+
   geom_text_repel(min.segment.length = 0,aes(label = rownames(plot.data)),
                   size = 6.5,color="black") +
   geom_smooth(method='lm', formula= y~x,colour="black",se = FALSE)+
   stat_cor(method = "pearson", label.x.npc = "left", label.y.npc ="bottom" ,size = 5)

g2<-ggplot(plot.data, aes(x=Pan.DNA.RNA.Corr,y=Pan.RNA.Pro.Corr)) +
   geom_point(aes(size=mRNA.halflife),shape = 16,colour="dodgerblue2") +
   theme_classic()+
   #scale_colour_gradient2(low = "#4A6FE3", high = "#D33F6A",mid = "#E2E2E2",midpoint = 12) +
   scale_size(range = c(1,15),name = "mRNA.halflife")+
   geom_text_repel(min.segment.length = 0,aes(label = rownames(plot.data)),
                   size = 6.5,color="black") +
   geom_smooth(method='lm', formula= y~x,colour="black",se = FALSE)+
   stat_cor(method = "pearson", label.x.npc = "left", label.y.npc ="bottom" ,size = 5)
g3<-ggplot(plot.data, aes(x=Pan.DNA.RNA.Corr,y=Normal.RNA.Pro.Corr)) +
   geom_point(shape = 16,colour="black",size=5) +
   theme_classic()+
   #scale_colour_gradient2(low = "#4A6FE3", high = "#D33F6A",mid = "#E2E2E2",midpoint = 12) +
   #scale_size(range = c(1,15),name = "mRNA.halflife")+
   geom_text_repel(min.segment.length = 0,aes(label = rownames(plot.data)),
                   size = 6.5,color="black") +
   geom_smooth(method='lm', formula= y~x,colour="black",se = FALSE)+
   stat_cor(method = "spearman", label.x.npc = "left", label.y.npc ="bottom" ,size = 5)
pdf("211026.pathway.corr.CPTAC_rna_protein_vs_CPTAC_rna_pro.g1.pdf",width = 14,height=12)
print(g2)
dev.off()


