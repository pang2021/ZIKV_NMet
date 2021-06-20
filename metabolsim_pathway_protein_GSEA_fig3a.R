options(stringsAsFactors=F)
setwd("E:/ZIKA/results/GSEA")
protein_for_1<-read.table(file="gsea_report_for_1_protein.txt",sep = "\t",header = T)
protein_for_0<-read.table(file="gsea_report_for_0_protein.txt",sep = "\t",header = T)
protein_gsea<-rbind(protein_for_1,protein_for_0)
metabolic_path<-read.table(file="metabolic_path.txt",sep = "\t",header = F)
metabolic_select<-read.table(file="E:/ZIKA/results/GSEA/protein/metabolism_select.txt",sep = "\t",header = T)
gsea_name<-c()
for(j in 1:dim(protein_gsea)[1]){
gsea_name[j]<-strsplit(protein_gsea$NAME," - M")[[j]][1]
}
index<-c()
for(i in 1:dim(metabolic_path)[1]){
index[i]<-match(toupper(metabolic_path[i,2]),gsea_name)
}
metabolic_gsea<-protein_gsea[index,]
metabolic_gsea<-na.omit(metabolic_gsea)
write.table(metabolic_gsea,file="metabolic_path_map.txt",sep = "\t",col.names = T,row.names = F,quote = F)
metabolic_gsea$NOM.p.val[metabolic_gsea$NOM.p.val<=0.0001]<-0.0001
metabolic_annotate<-metabolic_gsea[match(metabolic_select$path,metabolic_gsea$NAME),]
metabolic_gsea$NOM.p.val<- -log10(metabolic_gsea$NOM.p.val)
metabolic_annotate$NOM.p.val<- -log10(metabolic_annotate$NOM.p.val)




library(ggplot2)
library(ggrepel)
gg<-ggplot(data = metabolic_gsea, aes(x = NOM.p.val, y = NES,color=NES,label = NAME)) +
  geom_point(size=15*abs(metabolic_gsea$NES)) +
  scale_color_gradient2(low = "blue",high = "red",mid='white') +
  scale_size_area(max_size = 20)+
  labs(x="-log10 (p-value)",y="NES",title="metabolic_gsea") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    geom_text_repel(
    data =metabolic_annotate,
    aes(label = NAME),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = T )
gg+theme_bw()+theme(axis.line = element_line(size=1, colour = "black"),panel.border = element_blank()
,panel.grid.major=element_blank(),panel.grid.minor=element_blank())
