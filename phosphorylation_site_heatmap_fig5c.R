options(stringsAsFactors=F)
setwd('E:/ZIKA/results/datasets_compare')
pho_protein_file<-read.table(file="E:/ZIKA/data/datasets/pho_protein_post.txt",sep = "\t",header = T,quote="")
pho_protein_file[pho_protein_file==0]<-0.01
#
for(i in 1:dim(pho_protein_file)[1]){
index<-match(pho_protein_file$Gene.Symbol[i],rownames(protein_file))

if(!is.na(index)){
pho_protein_file[i,c(4:19)]<-pho_protein_file[i,c(4:19)]/protein_file[index,]
}else{
pho_protein_file[i,c(4:19)]<-pho_protein_file[i,c(4:19)]
}
}


load("E:/ZIKA/data/KEGG_patheay_gmt/KEGG_list.Rdata")
path_required<-c('MAPK signaling pathway','cGMP-PKG signaling pathway',
'cAMP signaling pathway','Rap1 signaling pathway','Ras signaling pathway',
'Calcium signaling pathway')
calculate_result<-data.frame()
for(j in 1:length(path_required)){
temp_genes<-path_list[[grep(path_required[j],names(path_list),fixed = T)]]

temp_frame<-pho_protein_file[pho_protein_file$Gene.Symbol%in%temp_genes,]
temp_name<-paste(temp_frame$Protein,temp_frame$Gene.Symbol,temp_frame$Phospho..STY..Site,sep='_')
temp_frame<-temp_frame[,c(4:19)]
temp_frame<-apply(temp_frame,2,as.numeric)
rownames(temp_frame)<-temp_name
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

calculate_result<-rbind(calculate_result,data.frame(rep(path_required[j],dim(all_frame_new)[1]),all_frame_new))
}

colnames(calculate_result)<-c('path','p_value','FC','log2_FC','p_fdr')

#calculate_result$metabolite<-metabolism_id[match(calculate_result$gene,metabolism_id$KEGG),1]
#calculate_result<-as.matrix(calculate_result)
#rownames(calculate_result)<-calculate_result[,11]
save(calculate_result,file='six_pathway_pho_normalized.Rdata')
write.table(calculate_result,file='six_pathway_pho_normalized.txt',sep='\t',row.names = T,col.names = T,quote = F)

##heatmap
options(stringsAsFactors=F)
setwd('E:/ZIKA/results/datasets_compare')
calculate_result<-read.table(file='heatmap_adjust.txt',sep='\t',header = T)
#calculate_result<-read.table(file='heatmap_unadjust.txt',sep='\t',header = T)
heatmatrix<-data.frame(calculate_result$log2_FC)
rownames(heatmatrix)<-calculate_result$pho
library(pheatmap)
pheatmap(heatmatrix,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white","#AF1E23"))(50),
show_colnames = T,main="pho_unadjust_heatmap")
