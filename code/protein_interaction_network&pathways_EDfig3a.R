options(stringsAsFactors=F)
options(stringsAsFactors=F)
setwd("E:/ZIKA/results/protein_cell_type")
load("summary_data.Rdata")
rownames(protein_file)<-toupper(rownames(protein_file))
load("cell_es.Rdata")
brain_cell<-c("Microglial cell","Neural stem cell","Astrocyte","Neuroblast","Neuron")
cell_file<-read.table(file="ZIKA_brain_cell.txt",sep = "\t",header = T,quote="")
cell_file$marker<-toupper(cell_file$marker)
cell_type<-names(table(cell_file$cell_type))

intersect_frame<-data.frame()
for(i in 1:length(cell_type)){
temp_frame<-data.frame(rep(cell_type[i],length(intersect(cell_list[[i]],rownames(protein_file)))),intersect(cell_list[[i]],rownames(protein_file)))
intersect_frame<-rbind(intersect_frame,temp_frame)
}
colnames(intersect_frame)<-c("cell_type","marker")

integrate_frame<-cbind(intersect_frame,protein_file[match(intersect_frame$marker,rownames(protein_file)),])
# integrate_frame<-rbind(integrate_frame[integrate_frame$cell_type%in%c("Microglial cell","Neural stem cell"),],
# integrate_frame[integrate_frame$cell_type%in%c("Astrocyte","Neuroblast","Neuron"),])
# annotate_rna<-integrate_frame[,1:2]

new_final_heat<-t(scale(t(integrate_frame[,3:8])))
new_final_heat<-cbind(new_final_heat[,4:6],new_final_heat[,1:3])

library(pheatmap)

pheatmap(t(new_final_heat),cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#113F8C","white","#AF1E23"))(50),
annotation_col =annotate_rna,show_colnames = F)

###
load("E:/ZIKA/data/KEGG_patheay_gmt/KEGG_list.Rdata")
metabolic_path<-read.table(file="E:/ZIKA/results/GSEA/metabolic_path.txt",sep = "\t",header = F)
observed_KEGG<-metabolic_path[,2]
all_kegg_name<-c()
for(i in 1:length(path_list)){
all_kegg_name[i]<-strsplit(names(path_list),' - M')[i][[1]]
}
observed_path<-path_list[match(observed_KEGG,all_kegg_name)]
observed_matrix<-list()
observed_profile<-data.frame()
for(i in 1:length(observed_path)){
observed_matrix[[i]]<-protein_file[intersect(rownames(protein_file),toupper(observed_path[[i]])),]
if(!is.null(dim(observed_matrix[[i]]))){
observed_profile<-rbind(observed_profile,data.frame(rep(observed_KEGG[i],dim(observed_matrix[[i]])[1]),
observed_matrix[[i]]))
}
}
colnames(observed_profile)[1]<-'path'
####
candidate_frame<-integrate_frame[integrate_frame$cell_type%in%brain_cell[5],]
R_value<-data.frame()
for(i in 1:dim(candidate_frame)[1]){
temp_row<-c()
for(j in 1:dim(observed_profile)[1]){
temp_row[j]<-cor.test(as.numeric(candidate_frame[i,3:18]),as.numeric(observed_profile[j,2:17]),
method = "spearman")$estimate
}
R_value<-rbind(R_value,temp_row)
}
R_value<-as.matrix(R_value)
rownames(R_value)<-candidate_frame$marker
colnames(R_value)<-rownames(observed_profile)
proper_value<-apply(R_value,1,function(x){return(x[x>=0.8])})
proper_dataframe<-data.frame()
for(i in 1:length(proper_value)){
temp_dataframe<-data.frame(rep(names(proper_value[i]),length(proper_value[[i]])),names(proper_value[[i]]),
as.numeric(proper_value[[i]]))
proper_dataframe<-rbind(proper_dataframe,temp_dataframe)
}
colnames(proper_dataframe)<-c("protein_1","protein_2","R_value")
proper_dataframe<-proper_dataframe[proper_dataframe$R_value!=1,]


library(org.Mm.eg.db)
k=keys(org.Mm.eg.db,keytype = "ENSEMBL")
list=select(org.Mm.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
list$SYMBOL<-toupper(list$SYMBOL)
final_entriz=na.omit(list[match(proper_dataframe$protein_2,list[,"SYMBOL"]),][,2])

library(clusterProfiler)
KEGG_gene<-enrichKEGG(final_entriz,organism = "mouse", pvalueCutoff = 0.05)
View(KEGG_gene@result)
write.table(KEGG_gene@result,file='E:/ZIKA/results/protein_cell_type/network/function/Neuron_KEGG.txt',
sep='\t',row.names=T,col.names=T,quote=F)

proper_KEGG<-KEGG_gene@result
KEGG_list<-list()
for(i in 1:dim(proper_KEGG)[1]){
KEGG_list[[i]]<-list[match(strsplit(proper_KEGG[i,]$geneID,"/")[[1]],list$ENTREZID),]$SYMBOL
}
names(KEGG_list)<-proper_KEGG$Description
Neuron_frame<-proper_dataframe[proper_dataframe$protein_2%in%unlist(KEGG_list),]
save(Neuron_frame,KEGG_list,KEGG_gene,file="Neuron_frame.Rdata")
write.table(Neuron_frame,file="network/Neuron_frame.txt",sep="\t",col.names=T,row.names = F,quote = F)

Neuron_annotate<-data.frame(c(names(table(Neuron_frame$protein_1)),
setdiff(names(table(Neuron_frame$protein_2)),names(table(Neuron_frame$protein_1)))),
c(rep(1,length(names(table(Neuron_frame$protein_1)))),
rep(0,length(setdiff(names(table(Neuron_frame$protein_2)),
names(table(Neuron_frame$protein_1)))))))
write.table(Neuron_annotate,file="network/Neuron_annotate.txt",sep="\t",col.names=T,row.names = F,quote = F)

