options(stringsAsFactors=F)
setwd("E:/ZIKA/results/sample_interaction")
proteinGroups<-read.table(file="protein_file.txt",sep = "\t",header = T)
protein_names<-proteinGroups[,1]
proteinGroups<-proteinGroups[,-1]
proteinGroups<-as.matrix(proteinGroups)
rownames(proteinGroups)<-protein_names

phosphoGroups<-read.table(file="pho_protein_file.txt",sep = "\t",header = T)
phospho_names<-phosphoGroups[,1]
phosphoGroups<-phosphoGroups[,-1]
phosphoGroups<-as.matrix(phosphoGroups)
rownames(phosphoGroups)<-phospho_names

mRNA_1_file<-read.table(file="E:/ZIKA/data/datasets/mRNA_1.txt",sep = "\t",header = T)
temp_name<-mRNA_1_file[,1]
mRNA_1_file<-mRNA_1_file[,-1]
mRNA_1_file<-as.matrix(mRNA_1_file)
rownames(mRNA_1_file)<-temp_name

metabolism_file<-read.table(file="metabolism_1.txt",sep = "\t",header = T,quote="")
temp_name<-metabolism_file[,1]
metabolism_file<-metabolism_file[,-1]
metabolism_file<-as.matrix(metabolism_file)
rownames(metabolism_file)<-temp_name

df<-data.frame(c(metabolism_file[,1],metabolism_file[,2],metabolism_file[,3],metabolism_file[,4],metabolism_file[,5],
metabolism_file[,6],metabolism_file[,7],metabolism_file[,8],metabolism_file[,9],metabolism_file[,10],
metabolism_file[,11],metabolism_file[,12],metabolism_file[,13],metabolism_file[,14],metabolism_file[,15],
metabolism_file[,16]),
c(rep("Ctr_1",length(metabolism_file[,1])),rep("Ctr_2",length(metabolism_file[,2])),
rep("Ctr_3",length(metabolism_file[,3])),rep("Ctr_4",length(metabolism_file[,4])),
rep("Ctr_5",length(metabolism_file[,4])),rep("Ctr_6",length(metabolism_file[,4])),
rep("Ctr_7",length(metabolism_file[,4])),rep("Ctr_8",length(metabolism_file[,4])),
rep("ZIKA_1",length(metabolism_file[,4])),rep("ZIKA_2",length(metabolism_file[,4])),
rep("ZIKA_3",length(metabolism_file[,5])),rep("ZIKA_4",length(metabolism_file[,6])),
rep("ZIKA_5",length(metabolism_file[,6])),rep("ZIKA_6",length(metabolism_file[,6])),
rep("ZIKA_7",length(metabolism_file[,6])),rep("ZIKA_8",length(metabolism_file[,6]))))
colnames(df)<-c("num","group")
df$num<-log10(df$num+1)


boxplot( num ~ group, df, col = c("#00AFBB", "#E7B800", "#FC4E07","#0000CD","#FFFF00","#8A2BE2"))

library(pheatmap)
final_frame<-pho_protein_file
R_frame<-data.frame()
for(i in 1:16){
row_R<-c()
for(j in 1:16){
r_value<-cor.test(final_frame[,i],final_frame[,j])$estimate
row_R<-c(row_R,r_value)
}
R_frame<-rbind(R_frame,row_R)
}
rownames(R_frame)<-colnames(final_frame)
colnames(R_frame)<-colnames(final_frame)

get_upper_tri <- function(cormat){
  cormat[upper.tri(cormat,F)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(R_frame)


pheatmap(R_frame,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#113F8C","white","#AF1E23"))(50),display_numbers=T,
show_rownames = T,show_colnames = T,main="phoprotein sample_correlate")
