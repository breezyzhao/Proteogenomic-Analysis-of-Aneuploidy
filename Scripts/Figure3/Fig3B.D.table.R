setwd("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/211026.control.test.fig3a.b/")
corum <- as.data.frame(read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/database/coreComplexesv3.0.txt", sep = "\t"))
corum<-corum[corum$Organism=="Human",]
master_list <- strsplit(as.character(corum$subunits.Entrez.IDs.), split = ";")
master_list <- unique(as.numeric(as.character(unlist(master_list))))
master_list_uniprot <- strsplit(as.character(corum$subunits.UniProt.IDs.), split = ";")
master_list_uniprot <- unique(as.character(unlist(master_list_uniprot)))
master_list_names <- strsplit(as.character(corum$subunits.Gene.name.), split = ";")
master_list_names <- unique(as.character(unlist(master_list_names)))

colon<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct/Colon.corr.table.txt",sep="\t")
colnames(colon)[2]<-"colon.cor.rna.dna"
colnames(colon)[5]<-"colon.cor.rna.pro"
colnames(colon)[7]<-"colon.cor.dna.pro"

breast<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct//BREAST.corr.table.txt",sep="\t")
colnames(breast)[2]<-"breast.cor.rna.dna"
colnames(breast)[5]<-"breast.cor.rna.pro"
colnames(breast)[7]<-"breast.cor.dna.pro"

ovarian<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct//OV.corr.table.txt",sep="\t")
colnames(ovarian)[2]<-"ovarian.cor.rna.dna"
colnames(ovarian)[5]<-"ovarian.cor.rna.pro"
colnames(ovarian)[7]<-"ovarian.cor.dna.pro"

ccRCC<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct//ccRCC.corr.table.txt",sep = "\t")
colnames(ccRCC)[2]<-"ccRCC.cor.rna.dna"
colnames(ccRCC)[5]<-"ccRCC.cor.rna.pro"
colnames(ccRCC)[7]<-"ccRCC.cor.dna.pro"

endometrial<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct//Endometrial.corr.table.txt",sep = "\t")
colnames(endometrial)[2]<-"endometrial.cor.rna.dna"
colnames(endometrial)[5]<-"endometrial.cor.rna.pro"
colnames(endometrial)[7]<-"endometrial.cor.dna.pro"

luad<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct//LUAD.corr.table.txt",sep = "\t")
colnames(luad)[2]<-"luad.cor.rna.dna"
colnames(luad)[5]<-"luad.cor.rna.pro"
colnames(luad)[7]<-"luad.cor.dna.pro"

hnsc<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct//HNSC.corr.table.txt",sep = "\t")
colnames(hnsc)[2]<-"hnsc.cor.rna.dna"
colnames(hnsc)[5]<-"hnsc.cor.rna.pro"
colnames(hnsc)[7]<-"hnsc.cor.dna.pro"

coad.var<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/COAD.gain.loss.var.table.txt",
                     sep="\t",header = T)
brca.var<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/BRCA.gain.loss.var.table.txt",
                     sep="\t",header = T)
ov.var<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/OV.gain.loss.var.table.txt",
                     sep="\t",header = T)
hnsc.var<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/HNSC.gain.loss.var.table.txt",
                     sep="\t",header = T)
luad.var<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/LUAD.gain.loss.var.table.txt",
                     sep="\t",header = T)
ccrcc.var<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/ccRCC.gain.loss.var.table.txt",
                     sep="\t",header = T)
ucec.var<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/UCEC.gain.loss.var.table.txt",
                     sep="\t",header = T)
pan.var<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210903.controltest/210924.gain.loss.variance.score/Pan.gain.loss.var.table.txt",
                     sep="\t",header = T)

pathway<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210901.complex.table>>move_to_fig3//Picked.pathway.txt",
                    sep="\t",header = T)
path.names<-unique(pathway$group)
pan.corr.info<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.2/210819.rmEXP10pct/Pan-Cancer.corr.table.txt",
                      sep="\t",header = T)
halflife<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/dataset_outside/half-life/halflife-protein.trans.txt",
                     sep="\t",header = T)
degradation<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/dataset_outside/degradation/RPE-1.txt",
                        sep="\t",header = T)
phylop<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/36.evolution/hg19.gene.score.txt",
                   sep="\t",header = T)
tf.binding<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210908.TF/210915.additional/msigdb.TF.v7.4.txt",
                        sep="\t",header=T)
tf.GTRD<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210908.TF/210915.additional/msigdb.TF.GTRD.v7.4.txt",
                    sep="\t",header = T)
tf.legacy<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210908.TF/210915.additional/msigdb.TF.legacy.v7.4.txt",
                      sep="\t",header = T)

biophysical.data<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210907.Jason_biophysical.test/biophysical.txt",
                             sep="\t",header = T)

ttrust.data<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210908.TF/210915.additional/trrust.v2.txt",
                        sep="\t",header = T)

homology.perc<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210918.homology.percentage/homology.percentage.txt",
                                sep="\t",header = T)

normal.data<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210830.normaldatabase/Multiple_tissue.corr.table.txt",
                        sep="\t",header = T)

hcec.data<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210827.hCEC/6.corr.table.txt",
                      sep="\t",header = T)
hcec.all<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/fig2/210827.hCEC/15.corr.table.txt",
                     sep="\t",header = T)

yeast<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210930.yeast/yeast.corr.table.txt",sep="\t",header = T)

yeast.rnascore<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/210930.yeast/yeast.rnascore.table.txt",sep="\t",header = T)

gtex<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/potential_analysis/211010.GTEx/gtex.corr.txt",
                 sep="\t",header = T)

ccle<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.1/210813/CCLE.Pan-Cancer.corr.table.txt",
                 sep="\t",header = T)

nci60<-read.delim("/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/10.protein/11.12-For paper/Manuscript_2021/figure_2021/figS2/figS2.1/210813/NCI-60.Cell Lines.corr.table.txt",sep="\t",header = T)


new.table<-data.frame(item=path.names,mRNA.halflife=rep(NA,length(path.names)),protein.halflife=rep(NA,length(path.names)),
                      Pan.DNA.RNA.Corr=rep(NA,length(path.names)),
                      Pan.DNA.RNA.Corr.SE.min=rep(NA,length(path.names)),Pan.DNA.RNA.Corr.SE.max=rep(NA,length(path.names)),
                      Pan.RNA.Pro.Corr=rep(NA,length(path.names)),
                      Pan.RNA.Pro.Corr.SE.min=rep(NA,length(path.names)), Pan.RNA.Pro.Corr.SE.max=rep(NA,length(path.names)),
                      COAD.DNA.RNA.Corr=rep(NA,length(path.names)),COAD.RNA.Pro.Corr=rep(NA,length(path.names)),
                      BRCA.DNA.RNA.Corr=rep(NA,length(path.names)),BRCA.RNA.Pro.Corr=rep(NA,length(path.names)),
                      OV.DNA.RNA.Corr=rep(NA,length(path.names)),OV.RNA.Pro.Corr=rep(NA,length(path.names)),
                      ccRCC.DNA.RNA.Corr=rep(NA,length(path.names)),ccRCC.RNA.Pro.Corr=rep(NA,length(path.names)),
                      HNSC.DNA.RNA.Corr=rep(NA,length(path.names)),HNSC.RNA.Pro.Corr=rep(NA,length(path.names)),
                      LUAD.DNA.RNA.Corr=rep(NA,length(path.names)),LUAD.RNA.Pro.Corr=rep(NA,length(path.names)),
                      UCEC.DNA.RNA.Corr=rep(NA,length(path.names)),UCEC.RNA.Pro.Corr=rep(NA,length(path.names)),
                      ED.pct=rep(NA,length(path.names)),NED.pct=rep(NA,length(path.names)),
                      ED.score=rep(NA,length(path.names)),NED.score=rep(NA,length(path.names)),
                      phyloP=rep(NA,length(path.names)),
                      N.binding.sites.of.TF.median=rep(NA,length(path.names)),
                      N.binding.sites.of.TF.GTRD.median=rep(NA,length(path.names)),
                      N.binding.sites.of.TF.legacy.median=rep(NA,length(path.names)),
                      N.binding.sites.TTRUST.TF.median=rep(NA,length(path.names)),
                      homology.percentage=rep(NA,length(path.names)),
                      Normal.RNA.Pro.Corr=rep(NA,length(path.names)),
                      Normal.RNA.Pro.Corr.SE.min=rep(NA,length(path.names)), Normal.RNA.Pro.Corr.SE.max=rep(NA,length(path.names)),
                      GTEx.RNA.Pro.Corr=rep(NA,length(path.names)),
                      CCLE.DNA.RNA.Corr=rep(NA,length(path.names)),
                      CCLE.RNA.Pro.Corr=rep(NA,length(path.names)),
                      NCI60.DNA.RNA.Corr=rep(NA,length(path.names)),
                      NCI60.RNA.Pro.Corr=rep(NA,length(path.names)),
                      hCEC.DNA.RNA.Corr=rep(NA,length(path.names)),
                        hCEC.RNA.Pro.Corr=rep(NA,length(path.names)),
                      hCEC.all.DNA.RNA.Corr=rep(NA,length(path.names)),
                      hCEC.all.RNA.Pro.Corr=rep(NA,length(path.names)),
                      CORUM.pct=rep(NA,length(path.names)),
                      N.of.genes.Pan.pathway=rep(NA,length(path.names)),
                      N.of.genes.total.pathway=rep(NA,length(path.names)),
                      yeast.DNA.RNA.Corr=rep(NA,length(path.names)),yeast.RNA.Pro.Corr=rep(NA,length(path.names)),
                      Pan.yeast.DNA.RNA.Corr=rep(NA,length(path.names)),Pan.yeast.RNA.Pro.Corr=rep(NA,length(path.names)),
                      Yeast.RNA.score=rep(NA,length(path.names)),
                      yeast.pathway.pct=rep(NA,length(path.names)))

datalist<-list()
for (i in 1:length(path.names)){
  genes<-pathway[pathway$group==path.names[i],]
  if(nrow(genes)==0){
    new.table[i]<-"NA"
  }else{
  half.use<-halflife[halflife$Symbol.x%in%genes$genes,]
  new.table$mRNA.halflife[i]<-median(half.use$mRNA_half.life_average_.h.,na.rm = T)
  new.table$protein.halflife[i]<-median(half.use$Protein_half.life_average_.h.,na.rm = T)
  
  pan.corr.use<-pan.corr.info[pan.corr.info$gene %in%genes$genes,]
  new.table$Pan.DNA.RNA.Corr[i]<-median(pan.corr.use$corr.rna.dna,na.rm=T)
  new.table$Pan.RNA.Pro.Corr[i]<-median(pan.corr.use$corr.rna.pro,na.rm=T)
  
  se.rna.dna<-sqrt(var(pan.corr.use$corr.rna.dna)/length(pan.corr.use$corr.rna.dna))
  new.table$Pan.DNA.RNA.Corr.SE.min[i]<-median(pan.corr.use$corr.rna.dna,na.rm=T)-round(se.rna.dna,3)
  new.table$Pan.DNA.RNA.Corr.SE.max[i]<-median(pan.corr.use$corr.rna.dna,na.rm=T)+round(se.rna.dna,3)
  
  se.rna.pro<-sqrt(var(pan.corr.use$corr.rna.pro)/length(pan.corr.use$corr.rna.pro))
  new.table$Pan.RNA.Pro.Corr.SE.min[i]<-median(pan.corr.use$corr.rna.pro,na.rm=T)-round(se.rna.pro,3)
  new.table$Pan.RNA.Pro.Corr.SE.max[i]<-median(pan.corr.use$corr.rna.pro,na.rm=T)+round(se.rna.pro,3)
  
  
  colon.corr.use<-colon[colon$gene %in%genes$genes,]
  new.table$COAD.DNA.RNA.Corr[i]<-median(colon.corr.use$colon.cor.rna.dna,na.rm=T)
  new.table$COAD.RNA.Pro.Corr[i]<-median(colon.corr.use$colon.cor.rna.pro,na.rm=T)
  
  breast.corr.use<-breast[breast$gene %in%genes$genes,]
  new.table$BRCA.DNA.RNA.Corr[i]<-median(breast.corr.use$breast.cor.rna.dna,na.rm=T)
  new.table$BRCA.RNA.Pro.Corr[i]<-median(breast.corr.use$breast.cor.rna.pro,na.rm=T)
  
  ov.corr.use<-ovarian[ovarian$gene %in%genes$genes,]
  new.table$OV.DNA.RNA.Corr[i]<-median(ov.corr.use$ovarian.cor.rna.dna,na.rm=T)
  new.table$OV.RNA.Pro.Corr[i]<-median(ov.corr.use$ovarian.cor.rna.pro,na.rm=T)
  
  ccrcc.corr.use<-ccRCC[ccRCC$gene %in%genes$genes,]
  new.table$ccRCC.DNA.RNA.Corr[i]<-median(ccrcc.corr.use$ccRCC.cor.rna.dna,na.rm=T)
  new.table$ccRCC.RNA.Pro.Corr[i]<-median(ccrcc.corr.use$ccRCC.cor.rna.pro,na.rm=T)
  
  endometrial.corr.use<-endometrial[endometrial$gene %in%genes$genes,]
  new.table$UCEC.DNA.RNA.Corr[i]<-median(endometrial.corr.use$endometrial.cor.rna.dna,na.rm=T)
  new.table$UCEC.RNA.Pro.Corr[i]<-median(endometrial.corr.use$endometrial.cor.rna.pro,na.rm=T)
  
  luad.corr.use<-luad[luad$gene %in%genes$genes,]
  new.table$LUAD.DNA.RNA.Corr[i]<-median(luad.corr.use$luad.cor.rna.dna,na.rm=T)
  new.table$LUAD.RNA.Pro.Corr[i]<-median(luad.corr.use$luad.cor.rna.pro,na.rm=T)
  
  hnsc.corr.use<-hnsc[hnsc$gene %in%genes$genes,]
  new.table$HNSC.DNA.RNA.Corr[i]<-median(hnsc.corr.use$hnsc.cor.rna.dna,na.rm=T)
  new.table$HNSC.RNA.Pro.Corr[i]<-median(hnsc.corr.use$hnsc.cor.rna.pro,na.rm=T)
  
  degra.use<-degradation[degradation$Gene.names %in% genes$genes,]
  if(nrow(degra.use[degra.use$Degradation.profile=="ED",])==0){
    new.table$ED.pct[i]<-0
    new.table$ED.score[i]<-0
  }else{
    new.table$ED.pct[i]<-nrow(degra.use[degra.use$Degradation.profile=="ED",])/nrow(degra.use)
    new.table$ED.score[i]<-median(degra.use[degra.use$Degradation.profile=="ED",]$X_.score,na.rm=T)
  }
  if(nrow(degra.use[degra.use$Degradation.profile=="NED",])==0){
    new.table$NED.pct[i]<-0
    new.table$NED.score[i]<-0
  }else{
  new.table$NED.pct[i]<-nrow(degra.use[degra.use$Degradation.profile=="NED",])/nrow(degra.use)
  new.table$NED.score[i]<-median(degra.use[degra.use$Degradation.profile=="NED",]$X_.score,na.rm=T)
  }
  
  
  phylop.use<-phylop[phylop$genes%in%genes$genes,]
  new.table$phyloP[i]<-median(phylop.use$median.PhyloP.score,na.rm=T)
  
  ##TF BINGING SITES
  tf.use<-tf.binding[tf.binding$genes%in%genes$genes,]
  new.table$N.binding.sites.of.TF.median[i]<-median(tf.use$Freq,na.rm=T)
  #new.table[i,12]<-sum(tf.use$Freq,na.rm=T)
  
  tf.gtrd.use<-tf.GTRD[tf.GTRD$genes%in%genes$genes,]
  new.table$N.binding.sites.of.TF.GTRD.median[i]<-median(tf.gtrd.use$Freq,na.rm=T)
  
  tf.legacy.use<-tf.legacy[tf.legacy$genes%in%genes$genes,]
  new.table$N.binding.sites.of.TF.legacy.median[i]<-median(tf.legacy.use$Freq,na.rm=T)
  
  trust<-ttrust.data[ttrust.data$Var1%in%genes$genes,]
  new.table$N.binding.sites.TTRUST.TF.median[i]<-median(trust$Freq,na.rm=T)
  
  homology<-homology.perc[homology.perc$Gene_Symbol%in%genes$genes,]
  new.table$homology.percentage[i]<-median(homology$content_fraction,na.rm = T)
  
  normal<-normal.data[normal.data$gene%in%genes$genes,]
  new.table$Normal.RNA.Pro.Corr[i]<-median(normal$RNA.PRO.Corr,na.rm = T)
  
  normal.se.rna.pro<-sqrt(var(normal$RNA.PRO.Corr)/length(normal$RNA.PRO.Corr))
  new.table$Normal.RNA.Pro.Corr.SE.min[i]<-median(normal$RNA.PRO.Corr,na.rm=T)-round(normal.se.rna.pro,3)
  new.table$Normal.RNA.Pro.Corr.SE.max[i]<-median(normal$RNA.PRO.Corr,na.rm=T)+round(normal.se.rna.pro,3)
  
  gtex.use<-gtex[gtex$gene%in%genes$genes,]
  new.table$GTEx.RNA.Pro.Corr[i]<-median(gtex.use$corr,na.rm = T)
  
  ccle.use<-ccle[ccle$gene%in%genes$genes,]
  new.table$CCLE.DNA.RNA.Corr[i]<-median(ccle.use$corr.rna.dna,na.rm = T)
  new.table$CCLE.RNA.Pro.Corr[i]<-median(ccle.use$corr.rna.pro,na.rm = T)
  
  nci60.use<-nci60[nci60$gene%in%genes$genes,]
  new.table$NCI60.DNA.RNA.Corr[i]<-median(nci60.use$corr.rna.dna,na.rm = T)
  new.table$NCI60.RNA.Pro.Corr[i]<-median(nci60.use$corr.rna.pro,na.rm = T)
  
  hcec.use<-hcec.data[hcec.data$gene %in%genes$genes,]
  new.table$hCEC.DNA.RNA.Corr[i]<-median(hcec.use$corr.rna.dna,na.rm=T)
  new.table$hCEC.RNA.Pro.Corr[i]<-median(hcec.use$corr.rna.pro,na.rm=T)
  
  hcec.all.use<-hcec.all[hcec.all$gene %in%genes$genes,]
  new.table$hCEC.all.DNA.RNA.Corr[i]<-median(hcec.all.use$corr.rna.dna,na.rm=T)
  new.table$hCEC.all.RNA.Pro.Corr[i]<-median(hcec.all.use$corr.rna.pro,na.rm=T)
  
  corum.length<-intersect(genes$genes,master_list_names)
  new.table$CORUM.pct[i]<-length(corum.length)/nrow(genes)
  
  new.table$N.of.genes.total.pathway[i]<-nrow(genes)
  new.table$N.of.genes.Pan.pathway[i]<-nrow(pan.corr.use)
  
  
  yeast.use<-yeast[yeast$Human.gene.name %in%genes$genes,]
  new.table$yeast.DNA.RNA.Corr[i]<-median(yeast.use$DNA.RNA_Corr,na.rm=T)
  new.table$yeast.RNA.Pro.Corr[i]<-median(yeast.use$RNA.Pro_Corr,na.rm=T)
  
  new.table$yeast.pathway.pct[i]<-nrow(yeast.use)
  
  Pan.yeast<-pan.corr.info[pan.corr.info$gene %in%genes$genes,]
  Pan.yeast.use<-Pan.yeast[Pan.yeast$gene %in% yeast$Human.gene.name,]
  new.table$Pan.yeast.DNA.RNA.Corr[i]<-median(Pan.yeast.use$corr.rna.dna,na.rm=T)
  new.table$Pan.yeast.RNA.Pro.Corr[i]<-median(Pan.yeast.use$corr.rna.pro,na.rm=T)
  
  yeast.rnascore.use<-yeast.rnascore[yeast.rnascore$Human.gene.name %in%genes$genes,]
  new.table$Yeast.RNA.score[i]<-median(yeast.rnascore.use$score ,na.rm=T)
  
  
  biopysical.use<-biophysical.data[biophysical.data$gene%in%genes$genes,]
  
  data1<-as.data.frame(t(apply(biopysical.use[,5:34],2,function(x)median(x,na.rm=T))))
  data1$path<-path.names[i]
  
  coad.use<-coad.var[coad.var$gene%in%genes$genes,]
  coad.temp<-as.data.frame(t(apply(coad.use[,2:6],2,function(x)median(x,na.rm=T))))
  colnames(coad.temp)<-paste0("COAD.",colnames(coad.temp))
  
  brca.use<-brca.var[brca.var$gene%in%genes$genes,]
  brca.temp<-as.data.frame(t(apply(brca.use[,2:6],2,function(x)median(x,na.rm=T))))
  colnames(brca.temp)<-paste0("BRCA.",colnames(brca.temp))
  
  ov.use<-ov.var[ov.var$gene%in%genes$genes,]
  ov.temp<-as.data.frame(t(apply(ov.use[,2:6],2,function(x)median(x,na.rm=T))))
  colnames(ov.temp)<-paste0("OV.",colnames(ov.temp))
  
  hnsc.use<-hnsc.var[hnsc.var$gene%in%genes$genes,]
  hnsc.temp<-as.data.frame(t(apply(hnsc.use[,2:6],2,function(x)median(x,na.rm=T))))
  colnames(hnsc.temp)<-paste0("HNSC.",colnames(hnsc.temp))
  
  luad.use<-luad.var[luad.var$gene%in%genes$genes,]
  luad.temp<-as.data.frame(t(apply(luad.use[,2:6],2,function(x)median(x,na.rm=T))))
  colnames(luad.temp)<-paste0("LUAD.",colnames(luad.temp))
  
  ccrcc.use<-ccrcc.var[ccrcc.var$gene%in%genes$genes,]
  ccrcc.temp<-as.data.frame(t(apply(ccrcc.use[,2:6],2,function(x)median(x,na.rm=T))))
  colnames(ccrcc.temp)<-paste0("ccRCC.",colnames(ccrcc.temp))
  
  ucec.use<-ucec.var[ucec.var$gene%in%genes$genes,]
  ucec.temp<-as.data.frame(t(apply(ucec.use[,2:6],2,function(x)median(x,na.rm=T))))
  colnames(ucec.temp)<-paste0("UCEC.",colnames(ucec.temp))
  
  pan.use<-pan.var[pan.var$gene%in%genes$genes,]
  pan.temp<-as.data.frame(t(apply(pan.use[,2:6],2,function(x)median(x,na.rm=T))))
  colnames(pan.temp)<-paste0("Pan.",colnames(pan.temp))
  
  temp<-cbind(data1,coad.temp,brca.temp,ov.temp,hnsc.temp,luad.temp,ccrcc.temp,ucec.temp,pan.temp)
  datalist[[i]]<-temp
  
  }
} 

new.biopysical<-do.call(rbind,datalist)
new.use<-cbind(new.table,new.biopysical)

all(new.use$item==new.use$path)
new.use2<-dplyr::select(new.use,-c(path))

write.table(new.use2,"comprehensive.table.211026.txt",sep="\t",row.names = F,quote = F)






