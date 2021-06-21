options(stringsAsFactors=F)
load("summary_data.Rdata")
rownames(protein_file)<-toupper(rownames(protein_file))
rownames(mRNA_1_file)<-toupper(rownames(mRNA_1_file))

###make cell type gene lists
cell_file<-read.table(file="ZIKA_brain_cell.txt",sep = "\t",header = T,quote="")
cell_file$marker<-toupper(cell_file$marker)
cell_type<-names(table(cell_file$cell_type))
cell_list<-list()
for(i in 1:length(cell_type)){
cell_list[[i]]<-unique(cell_file[cell_file$cell_type%in%cell_type[i],]$marker)
}
names(cell_list)<-cell_type


library(GSVA)
cell_es_1 <- gsva(mRNA_1_file, cell_list,mx.diff=F, verbose=FALSE, parallel.sz=0)
cell_es_2 <- gsva(protein_file, cell_list,mx.diff=F, verbose=FALSE, parallel.sz=0)
save(cell_es_1,cell_es_2,cell_list,file="cell_es.Rdata")

###buble plot show active cell type
options(stringsAsFactors=F)
setwd("E:/ZIKA/results/protein_cell_type")
load("cell_es.Rdata")
brain_cell<-c("Microglial cell","Neural stem cell","Astrocyte","Neuroblast","Neuron")
##RNA
cell_es_1<-cell_es_1[brain_cell,]
buble_frame<-data.frame()
for(i in 1:dim(cell_es_1)[2]){
temp_frame<-data.frame(rownames(cell_es_1),cell_es_1[,i],rep(colnames(cell_es_1)[i],dim(cell_es_1)[1]),
c(rep("brain_cell",5)))
buble_frame<-rbind(buble_frame,temp_frame)
}
colnames(buble_frame)<-c("cell_type","es","sample","class")
buble_frame$cell_type<-factor(buble_frame$cell_type,levels=c("Microglial cell","Neural stem cell",
"Astrocyte","Neuroblast","Neuron"))
index<-c()
for(i in 1:dim(buble_frame)[1]){
if(buble_frame[i,2]>0){
temp_index<-"es > 0"
}else{
temp_index<-"es < 0"
}
index<-c(index,temp_index)
}
buble_frame$index<-index
buble_frame$index<-factor(buble_frame$index,levels=c("es > 0","es < 0"))


##protein
cell_es_2<-cell_es_2[intersect(brain_cell,rownames(cell_es_2)),]
buble_frame<-data.frame()
for(i in 1:dim(cell_es_2)[2]){
temp_frame<-data.frame(rownames(cell_es_2),cell_es_2[,i],rep(colnames(cell_es_2)[i],dim(cell_es_2)[1]),
c(rep("brain_cell",5)))
buble_frame<-rbind(buble_frame,temp_frame)
}
colnames(buble_frame)<-c("cell_type","es","sample","class")
buble_frame$cell_type<-factor(buble_frame$cell_type,levels=c("Microglial cell","Neural stem cell",
"Endothelial cell","Astrocyte","Neuroblast","Neuron"))
index<-c()
for(i in 1:dim(buble_frame)[1]){
if(buble_frame[i,2]>0){
temp_index<-"es > 0"
}else{
temp_index<-"es < 0"
}
index<-c(index,temp_index)
}
buble_frame$index<-index
buble_frame$index<-factor(buble_frame$index,levels=c("es > 0","es < 0"))



library(ggplot2)
p = ggplot(buble_frame,aes(cell_type,sample))
pr = p+ geom_point(aes(size=es,color=index))

pr + theme_bw()

