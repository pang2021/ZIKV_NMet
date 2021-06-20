options(stringsAsFactors=F)
setwd("E:/ZIKA/results/GSEA")
annote_info<-read.table(file="E:/ZIKA/results/annotate_GO_select.txt",sep = "\t",header = T)
attach(annote_info)
annote_info$type_2d[mRNA_1>0&protein>0] <- 1
annote_info$type_2d[mRNA_1>0&protein<0] <- 2
annote_info$type_2d[mRNA_1<0&protein>0] <- 3
annote_info$type_2d[mRNA_1<0&protein<0] <- 4
detach(annote_info)

detach(annote_info)
mRNA_1_for_1<-read.table(file="GO_report_for_1_mRNA_1.txt",sep = "\t",header = T,quote="")
mRNA_1_for_0<-read.table(file="GO_report_for_0_mRNA_1.txt",sep = "\t",header = T,quote="")
protein_for_1<-read.table(file="GO_report_for_1_protein.txt",sep = "\t",header = T,quote="")
protein_for_0<-read.table(file="GO_report_for_0_protein.txt",sep = "\t",header = T,quote="")
##
mRNA_1_gsea<-rbind(mRNA_1_for_1,mRNA_1_for_0)
mRNA_1_gsea<-mRNA_1_gsea[mRNA_1_gsea$NOM.p.val<=0.01,]
protein_gsea<-rbind(protein_for_1,protein_for_0)

NES_file_2d<-data.frame(mRNA_1_gsea$NES,
protein_gsea[match(mRNA_1_gsea$NAME,protein_gsea$NAME),]$NES)
rownames(NES_file_2d)<-mRNA_1_gsea$NAME
colnames(NES_file_2d)<-c("mRNA_1","protein")


##2D_type
attach(NES_file_2d)
NES_file_2d$type_2d[mRNA_1>0&protein>0] <- 1
NES_file_2d$type_2d[mRNA_1>0&protein<0] <- 2
NES_file_2d$type_2d[mRNA_1<0&protein>0] <- 3
NES_file_2d$type_2d[mRNA_1<0&protein<0] <- 4
detach(NES_file_2d)

##画图
library(ggplot2)
gg<-ggplot(NES_file_2d,aes(x=mRNA_1,y=protein))+geom_point(size=3,color=NES_file_2d$type_2d+1)+annotate("text",
x=annote_info$mRNA_1,y=annote_info$protein,label=annote_info$pathway,size=4,color=annote_info$type_2d+1)
gg<-gg+ geom_hline(aes(yintercept=0), colour="black")+geom_vline(aes(xintercept=0), colour="black")
gg+theme_bw()+theme(axis.line = element_line(size=1, colour = "black"),panel.border = element_blank()
,panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none")
