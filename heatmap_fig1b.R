options(stringsAsFactors=F)
load("GO_file.Rdata")
GO_mRNA_1_frame<-data.frame()
gaps_mRNA1<-c()
all_gene_mRNA1<-data.frame()
for(i in 1:length(GO_mRNA_1)){
GO_mRNA_1_frame<-rbind(GO_mRNA_1_frame,GO_mRNA_1[[i]])
gaps_mRNA1[i]<-dim(GO_mRNA_1_frame)[1]
gene_mRNA1<-data.frame(rep(names(GO_mRNA_1[i]),dim(GO_mRNA_1[[i]])[1]))
all_gene_mRNA1<-rbind(all_gene_mRNA1,gene_mRNA1)
}
rownames(all_gene_mRNA1)<-rownames(GO_mRNA_1_frame)

library(pheatmap)
pheatmap(scale(t(GO_mRNA_1_frame)),cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#113F8C","white","#AF1E23"))(50),
annotation_col = all_gene_mRNA1,gaps_row = 3,show_colnames = F,main="mRNA_1_pathway")


GO_protein_frame<-data.frame()
gaps_protein<-c()
all_gene_protein<-data.frame()
for(i in 1:length(GO_protein)){
GO_protein_frame<-rbind(GO_protein_frame,GO_protein[[i]])
gaps_protein[i]<-dim(GO_protein_frame)[1]
gene_protein<-data.frame(rep(names(GO_protein[i]),dim(GO_protein[[i]])[1]))
all_gene_protein<-rbind(all_gene_protein,gene_protein)
}
rownames(all_gene_protein)<-rownames(GO_protein_frame)
final_frame<-scale(t(GO_protein_frame))


pheatmap(final_frame,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#113F8C","white","#AF1E23"))(50),
gaps_row = 3,annotation_col = all_gene_protein,show_colnames = F,main="protein_pathway")


GO_pho_protein_frame<-data.frame()
gaps_pho_protein<-c()
all_gene_pho_protein<-data.frame()
for(i in 1:length(GO_pho_protein)){
GO_pho_protein_frame<-rbind(GO_pho_protein_frame,GO_pho_protein[[i]])
gaps_pho_protein[i]<-dim(GO_pho_protein_frame)[1]
gene_pho_protein<-data.frame(rep(names(GO_pho_protein[i]),dim(GO_pho_protein[[i]])[1]))
all_gene_pho_protein<-rbind(all_gene_pho_protein,gene_pho_protein)
}
rownames(all_gene_pho_protein)<-rownames(GO_pho_protein_frame)

pheatmap(scale(t(GO_pho_protein_frame)),cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#113F8C","white","#AF1E23"))(50),
gaps_row = 3,annotation_col = all_gene_pho_protein,show_colnames = F,main="pho_protein_pathway")

GO_metablism_frame<-data.frame()
gaps_metablism<-c()
all_gene_metablism<-data.frame()
for(i in 1:length(GO_metablism)){
temp_metabolism_name<-rownames(GO_metablism[[i]])
GO_metablism[[i]]<-apply(GO_metablism[[i]],2,as.numeric)
rownames(GO_metablism[[i]])<-temp_metabolism_name
GO_metablism_frame<-rbind(GO_metablism_frame,GO_metablism[[i]])
gaps_metablism[i]<-dim(GO_metablism_frame)[1]
gene_metablism<-data.frame(rep(names(GO_metablism[i]),dim(GO_metablism[[i]])[1]))
all_gene_metablism<-rbind(all_gene_metablism,gene_metablism)
}
rownames(all_gene_metablism)<-rownames(GO_metablism_frame)

final_frame<-scale(t(GO_metablism_frame))

pheatmap(final_frame,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#113F8C","white","#AF1E23"))(50),
gaps_row = 3,annotation_col = all_gene_metablism,show_colnames = F,main="metablism_pathway")
###integrate all heatmap
all_frame<-rbind(rbind(rbind(as.matrix(GO_mRNA_1_frame),as.matrix(GO_protein_frame)),
as.matrix(GO_pho_protein_frame)),as.matrix(GO_metablism_frame))
all_frame<-as.data.frame(all_frame)

all_gene_annotate<-rbind(rbind(rbind(as.matrix(all_gene_mRNA1),as.matrix(all_gene_protein)),
as.matrix(all_gene_pho_protein)),as.matrix(all_gene_metablism))
all_gene_annotate<-as.data.frame(all_gene_annotate)

all_gaps<-c(dim(all_gene_mRNA1)[1],
dim(all_gene_mRNA1)[1]+dim(all_gene_protein)[1],
dim(all_gene_mRNA1)[1]+dim(all_gene_protein)[1]+dim(all_gene_pho_protein)[1])

final_frame<-scale(t(all_frame))
final_frame[final_frame< -1]<- -1
final_frame[final_frame> 1]<- 1

pheatmap(final_frame,cluster_cols = F,cluster_rows = F,gaps_row = 3,
color = colorRampPalette(c("#113F8C","white","#AF1E23"))(50),
gaps_col = all_gaps,annotation_col = all_gene_annotate,show_colnames = F,main="multi-omics pathway")
