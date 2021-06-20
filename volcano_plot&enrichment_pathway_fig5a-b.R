options(stringsAsFactors=F)
setwd("E:/ZIKA/results/mRNA_1_deg")
#load('summary_data_pho_adjust.Rdata')
load("summary_data.Rdata")
pho_protein_file[pho_protein_file==0]<-0.01
temp_frame<-data.frame()
for(j in 1:dim(pho_protein_file)[1]){
if(mean(pho_protein_file[j,])!=0.01){
temp_value<-t.test(pho_protein_file[j,1:8],pho_protein_file[j,9:16],var.equal = T)
}else{
temp_value<-data.frame(1,c(1,1))
colnames(temp_value)<-c('p.value','estimate')
}
temp_frame<-rbind(temp_frame,c(temp_value$p.value[1],log2(temp_value$estimate[2]/temp_value$estimate[1])))
}
colnames(temp_frame)<-c("p_value","log2_FC")
temp_frame$p_value_fdr<-p.adjust(temp_frame$p_value,method='fdr')
temp_frame$p_value_fdr<- -log10(temp_frame$p_value_fdr)
temp_frame$p_value<- -log10(temp_frame$p_value)

fc=1
temp_frame$threshold<-0
temp_frame$threshold[temp_frame$log2_FC > fc & temp_frame$p_value_fdr> -log10(0.05)] = "up"
temp_frame$threshold[temp_frame$log2_FC < -fc & temp_frame$p_value_fdr> -log10(0.05)] = "down"
temp_frame$id<-rownames(protein_file)

new_frame<-cbind(protein_file,temp_frame)
#######volcano
library(ggplot2)
library(ggrepel)
gg<-ggplot(data = temp_frame, aes(x = log2_FC, y = p_value_fdr,color=threshold,label = id)) +
  geom_point(size=3) +
  geom_vline(xintercept=c(-fc,fc),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2_FC",y="log10_p",title="protein_deg") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    geom_text_repel(
    data = temp_frame[temp_frame$p_value_fdr>=2.193226,],
    aes(label = id),
    size = 5,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
gg+theme_bw()+theme(axis.line = element_line(size=1, colour = "black"),panel.border = element_blank()
,panel.grid.major=element_blank(),panel.grid.minor=element_blank())
###pathway enrichment
up_genes<-unique(temp_frame[with(temp_frame, (p_value > -log10(0.05) & log2_FC> 1)), ]$id)
down_genes<-unique(temp_frame[with(temp_frame, (p_value > -log10(0.05) & log2_FC< -1)), ]$id)

library(org.Mm.eg.db)
k=keys(org.Mm.eg.db,keytype = "ENSEMBL")
list=select(org.Mm.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
up_genes_entriz=na.omit(list[match(up_genes,list[,"SYMBOL"]),][,2])
down_genes_entriz=na.omit(list[match(down_genes,list[,"SYMBOL"]),][,2])

##up
library(clusterProfiler)
KEGG_gene<-enrichKEGG(up_genes_entriz,organism = "mouse", pvalueCutoff = 0.05)
GO_gene<-enrichGO(up_genes_entriz,'org.Mm.eg.db',ont = "BP", pvalueCutoff = 0.05)
View(KEGG_gene@result)
View(GO_gene@result)
write.table(KEGG_gene@result,file='up_pho_protein_kegg_p0.05&fc2.txt',sep = '\t',row.names = T,col.names = T,quote = F)
write.table(GO_gene@result,file='up_pho_protein_go_p0.05&fc2.txt',sep = '\t',row.names = T,col.names = T,quote = F)

##down
KEGG_gene<-enrichKEGG(down_genes_entriz,organism = "mouse", pvalueCutoff = 0.05)
GO_gene<-enrichGO(down_genes_entriz,'org.Mm.eg.db',ont = "BP", pvalueCutoff = 0.05)
View(KEGG_gene@result)
View(GO_gene@result)
write.table(KEGG_gene@result,file='down_pho_protein_kegg_p0.05&fc2.txt',sep = '\t',row.names = T,col.names = T,quote = F)
write.table(GO_gene@result,file='down_pho_protein_go_p0.05&fc2.txt',sep = '\t',row.names = T,col.names = T,quote = F)
