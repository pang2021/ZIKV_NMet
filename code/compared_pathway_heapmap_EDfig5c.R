options(stringsAsFactors=F)
setwd('E:/ZIKA/results/datasets_compare')
load('E:/ZIKA/data/datasets/summary_data.Rdata')
load("E:/ZIKA/data/datasets/pre_summary_data.Rdata")

######pathway
load("E:/ZIKA/data/KEGG_patheay_gmt/KEGG_list.Rdata")
annotate_metablism<-read.table(file="E:/ZIKA/results/heatmap/annotate_metablism.txt",sep = "\t",header = T,quote="")
path_required<-c('Citrate cycle (TCA cycle)','Nicotinate and nicotinamide metabolism',
'Oxidative phosphorylation',
'Tryptophan metabolism')
path_list<-list()
for(i in 1:length(unique(annotate_metablism$pathway))){
path_list[[i]]<-annotate_metablism[grep(unique(annotate_metablism$pathway)[i],
annotate_metablism$pathway,fixed = T),]$KEGG
}
names(path_list)<-unique(annotate_metablism$pathway)


######calculate p and fc
calculate_result<-data.frame()
for(j in 1:length(path_required)){
temp_genes<-path_list[[grep(path_required[j],names(path_list),fixed = T)]]

temp_frame<-metabolism_file
all_frame_new<-data.frame()
for(i in 1:dim(temp_frame)[1]){
if(mean(temp_frame[i,])!=0.01){
temp_frame_new<-t.test(temp_frame[i,1:8],temp_frame[i,9:16],var.equal = T)
}else{
temp_frame_new<-data.frame(1,c(1,1))
colnames(temp_frame_new)<-c('p.value','estimate')
}
all_frame_new<-rbind(all_frame_new,c(temp_frame_new$p.value,temp_frame_new$estimate[2]/temp_frame_new$estimate[1],
log2(temp_frame_new$estimate[2]/temp_frame_new$estimate[1])))
}
colnames(all_frame_new)<-c('p_value','FC','log2_FC')
all_frame_new$p_fdr<-p.adjust(all_frame_new$p_value,method='fdr')
all_frame_new<-as.matrix(all_frame_new)
rownames(all_frame_new)<-rownames(temp_frame)

temp_frame<-pre_metabolism_file
all_frame_old<-data.frame()
for(i in 1:dim(temp_frame)[1]){
if(mean(temp_frame[i,])!=0.01){
temp_frame_new<-t.test(temp_frame[i,1:8],temp_frame[i,9:15],var.equal = T)
}else{
temp_frame_new<-data.frame(1,c(1,1))
colnames(temp_frame_new)<-c('p.value','estimate')
}
all_frame_old<-rbind(all_frame_old,c(temp_frame_new$p.value,temp_frame_new$estimate[2]/temp_frame_new$estimate[1],
log2(temp_frame_new$estimate[2]/temp_frame_new$estimate[1])))
}
colnames(all_frame_old)<-c('p_value','FC','log2_FC')
all_frame_old$p_fdr<-p.adjust(all_frame_old$p_value,method='fdr')
all_frame_old<-as.matrix(all_frame_old)
rownames(all_frame_old)<-rownames(temp_frame)

two_frame_intersect<-union(rownames(all_frame_new),rownames(all_frame_old))
##
# temp_two_frame_intersect<-unlist(lapply(strsplit(two_frame_intersect,'_'),function(x){x[[1]]}))
# temp_a<-data.frame()
# for(i in 1:length(intersect(temp_genes,temp_two_frame_intersect))){
# a<-all_frame_new[grep(intersect(temp_genes,temp_two_frame_intersect)[i],rownames(all_frame_new)),]
# temp_a<-rbind(temp_a,a)
# }
if(length(intersect(temp_genes,two_frame_intersect))>1){
all_frame<-cbind(intersect(temp_genes,two_frame_intersect),as.data.frame(all_frame_new)[intersect(temp_genes,two_frame_intersect),],as.data.frame(all_frame_old)[intersect(temp_genes,two_frame_intersect),])
}else if(length(intersect(temp_genes,
two_frame_intersect))==1){
all_frame<-cbind(intersect(temp_genes,two_frame_intersect),t(all_frame_new[rownames(all_frame_new)%in%intersect(temp_genes,two_frame_intersect),]),t(as.data.frame(all_frame_new)[intersect(temp_genes,two_frame_intersect),],as.data.frame(all_frame_old)[intersect(temp_genes,two_frame_intersect),]))
}else if(length(intersect(temp_genes,two_frame_intersect))==0){
all_frame<-data.frame(t(rep(NA,9)))
}
colnames(all_frame)<-c('gene','new_p_value','new_FC','new_log2_FC','new_p_fdr',
'old_p_value','old_FC','old_log2_FC','old_p_fdr')
rownames(all_frame)<-intersect(temp_genes,two_frame_intersect)

calculate_result<-rbind(calculate_result,data.frame(rep(path_required[j],dim(all_frame)[1]),all_frame))
}

colnames(calculate_result)<-c('path','gene','new_p_value','new_FC','new_log2_FC','new_p_fdr',
'old_p_value','old_FC','old_log2_FC','old_p_fdr')

#calculate_result$metabolite<-metabolism_id[match(calculate_result$gene,metabolism_id$KEGG),1]
#calculate_result<-as.matrix(calculate_result)
#rownames(calculate_result)<-calculate_result[,11]
save(calculate_result,file='metabolism_calculate_result.Rdata')

##########
options(stringsAsFactors=F)
setwd('E:/ZIKA/results/datasets_compare')
load('protein_calculate_result.Rdata')
write.table(calculate_result,file='protein_calculate_result.txt',sep='\t',row.names = T,col.names = T,quote = F)
##heatmap
calculate_result<-read.table(file='metabolism_select.txt',sep='\t',header=T,na.strings = 'NA')
heatmatrix<-calculate_result[,c(5,9)]
heatmatrix[is.na(heatmatrix)]<-0
row_name<-matrix(calculate_result[,1])
rownames(row_name)<-calculate_result[,11]
heatmatrix<-as.matrix(heatmatrix)
heatmatrix<-apply(heatmatrix,2,as.numeric)
rownames(heatmatrix)<-calculate_result[,11]
library(pheatmap)
heatmatrix[heatmatrix>3.41]<-3.41
heatmatrix[heatmatrix< -3.41]<- -3.41
pheatmap(heatmatrix,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#113F8C","white","#AF1E23"))(50),
annotation_row = as.data.frame(row_name),show_colnames = T,main="metabolism_heatmap")
