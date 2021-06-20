library(xlsx)
library(xlsxjars)
library(rJava)
library(pheatmap)
options(stringsAsFactors = F)

## 
df <-read.xlsx("E:/R/PHH/heatmap_metabolomics_EDfig5a.xlsx", 
              sheetIndex = 1, startRow = 1, as.data.frame = T, header = T)
rownames(df) <- df$Metabolite
df <- df[,-1]


## 
#rg = colorRampPalette(c("lightcoral", "white", "lightblue4"))(1000) #navy firebrick3

rg = c(rep('firebrick', times = 50),  #"#4575B4" brown lightcoral orangered4 red4 firebrick3 firebrick
       colorRampPalette(c("firebrick", "white", "dodgerblue4"))(45),
       rep("dodgerblue4", times = 50)) #"#D73027" dodgerblue4 lightblue4 mediumblue steelblue4

## 
p <- pheatmap(df, legend = T, scale = "row", 
              
              #
              color = rev(rg), border_color = 'grey90', #'grey80'
              # 
              cellwidth = 2.6, cellheight = 2.25, #cex = 2.5,
              
              legend_breaks = c(-3, -1, 1, 3), legend_labels = c("-3","-1", "1","3"),
              
              # 
              #annotation_col = rev(annotation_col),
              #annotation_colors = annotation_colors,
              
              fontsize = 7, fontsize_row = 2.5, fontsize_col = 2.5, treeheight_row = 10, 
              annotation_names_row = T, annotation_names_col = F,
              
              # 
              clustering_distance_rows = "euclidean", clustering_method = "average", 
              cluster_cols = F, cluster_rows = T,
              
              # title
              main = '')
