options(stringsAsFactors=F)
setwd("E:/ZIKA/results/mRNA_1_deg")
#load('summary_data_pho_adjust.Rdata')
load("summary_data.Rdata")
protein_min<-min(protein_file[protein_file!=0])
protein_file[protein_file==0]<-protein_min

#######pho_protein vs protein
unadjust_pho_frame<-read.table(file='mapk_unadjust.txt',sep='\t',header=T)

for(i in 1:16){
temp_protein_file<-protein_file[unadjust_pho_frame$Gene.Symbol[i],]
new_protein_file<-data.frame(temp_protein_file,c(rep(0,8),rep(1,8)))
colnames(new_protein_file)<-c('num','group')
temp_pho_file<-unadjust_pho_frame[i,]
temp_pho_file<-as.numeric(temp_pho_file[,4:19])
new_pho_file<-data.frame(temp_pho_file,c(rep(0,8),rep(1,8)))
colnames(new_pho_file)<-c('num','group')
all_frame<-rbind(new_protein_file,new_pho_file)
all_frame$color_group<-c(rep('protein',16),rep('pho_protein',16))

setwd('E:/ZIKA/results/mRNA_1_deg/mapk_protein_pho')
library(ggplot2)
ggplot(data=all_frame, aes(x=group, y=num,col=as.factor(all_frame$color_group)))+geom_point(size=3)+stat_smooth(method="lm")+
labs(title=paste('pho_protein ',unadjust_pho_frame$Gene.Symbol[i],unadjust_pho_frame$Phospho..STY..Site[i]))+
theme_bw()+theme(axis.line = element_line(size=1, colour = "black"),panel.border = element_blank()
,panel.grid.major=element_blank(),panel.grid.minor=element_blank())
paste(unadjust_pho_frame$Gene.Symbol[i],unadjust_pho_frame$Phospho..STY..Site[i])
}