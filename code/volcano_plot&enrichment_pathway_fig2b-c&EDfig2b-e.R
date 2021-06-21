options(stringsAsFactors=F)
setwd("E:/ZIKA/results/mRNA_1_deg")
#load('summary_data_pho_adjust.Rdata')
load("summary_data.Rdata")
protein_min<-min(protein_file[protein_file!=0])
protein_file[protein_file==0]<-protein_min

kegg_select<-read.table(file='kegg_select.txt',sep='\t',header=T)
##protein
temp_frame<-data.frame()
for(j in 1:dim(protein_file)[1]){
if(mean(protein_file[j,])!=protein_min){
temp_value<-t.test(protein_file[j,1:8],protein_file[j,9:16],var.equal = T)
}else{
temp_value<-data.frame(1,c(1,1))
colnames(temp_value)<-c('p.value','estimate')
}
temp_frame<-rbind(temp_frame,c(temp_value$p.value[1],log2(temp_value$estimate[2]/temp_value$estimate[1])))
}
rownames(temp_frame)<-rownames(protein_file)
colnames(temp_frame)<-c("p_value","log2_FC")
temp_frame$p_value_fdr<-p.adjust(temp_frame$p_value,method='fdr')
temp_frame$p_value_fdr<- -log10(temp_frame$p_value_fdr)


fc=1
temp_frame$threshold<-0
temp_frame$threshold[temp_frame$log2_FC > fc & temp_frame$p_value_fdr> -log10(0.05)] = "up"
temp_frame$threshold[temp_frame$log2_FC < -fc & temp_frame$p_value_fdr> -log10(0.05)] = "down"
temp_frame$id<-rownames(protein_file)

new_frame<-cbind(protein_file,temp_frame)
write.table(new_frame,file='protein_deg.txt',sep = '\t',row.names = T,col.names = T,quote = F)
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
###pathway 
up_genes<-rownames(temp_frame[with(temp_frame, (p_value_fdr > -log10(0.05) & log2_FC> 1)), ])
down_genes<-rownames(temp_frame[with(temp_frame, (p_value_fdr > -log10(0.05) & log2_FC< -1)), ])

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
###kegg_select_heatmap
all_genes<-unlist(strsplit(GO_gene@result[match(kegg_select[4,1],GO_gene@result$Description),'geneID'],'/'))
all_genes_symbol=na.omit(list[match(all_genes,list[,"ENTREZID"]),][,3])
new_frame<-cbind(protein_file[all_genes_symbol,],temp_frame[all_genes_symbol,c(1:3)])
write.table(new_frame,file='axon development.txt',sep='\t',row.names = T,col.names = T,quote = F)
library(pheatmap)
a<-t(scale(t(new_frame[,1:16])))
pheatmap(a,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#113F8C","white","#AF1E23"))(50),
show_rownames = T,show_colnames = T,main="axon development")