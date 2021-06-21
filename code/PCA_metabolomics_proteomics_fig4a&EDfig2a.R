library(xlsx)
library(xlsxjars)
library(rJava)
library(ggplot2)
library(ggrepel)

#### 
df.PCA <-read.xlsx("E:/R/PHH/PCA_metabolomics_fig4a.xlsx", 
              sheetIndex = 1, startRow = 1, as.data.frame = T, header = T)
# 
rownames(df.PCA) <- df.PCA$Metabolite
df.PCA <- data.frame(t(df.PCA[,-1]), check.names = F)

#### 
pca_res <- prcomp(t(df.PCA), scale.=F)

## 
apply(pca_res$x, 2, sd)

## 
pca.group <-data.frame(pca_res$x, group = c(rep("Ctrl", 8), rep("ZIKV", 8)))

## 
pca_res.var <- apply(pca_res$x, 2, var)
# pca_res.var <- pca_res$sdev^2
percentage <- round(pca_res.var/sum(pca_res.var) * 100, 2)
percentage <- paste(colnames(pca.group)," (", as.character(percentage), "%", ")", sep = "")

pca.group$text <- rownames(pca.group)

# 
ggplot(pca.group, aes(x = PC1, y = PC2, color = group))+ 
  geom_point(size = 2)+  
  
  xlab(percentage[1]) +
  ylab(percentage[2])+
  
  stat_ellipse(level = 0.95, show.legend = T)+
  theme_bw() +
  # 
  # coord_fixed(ratio = 1) +
  
  #
  # theme(panel.border=element_rect(linetype = "solid",colour = "black"),
  #       # panel.grid.major=element_blank(),
  #       panel.grid.minor=element_blank(),
  #       axis.line= element_line(colour = "black")) +
  # 
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  
  # 
  coord_fixed(ratio = 1) +
  
  # 
  theme(axis.title.x =element_text(vjust = -2, hjust = 0.55, size=9, face="plain"),
        axis.title.y=element_text(vjust = 5, size=9, face="plain")) +
  # 
  theme(axis.text.y = element_text(size = 9,face="plain", colour = "black"),
        axis.text.x = element_text(size = 9,face="plain", colour = "black")) +

  # 
  geom_text_repel(data = pca.group, aes(x = PC1,
                                            y = PC2,
                                            label = text),
                  size = 3, fontface = 'plain',
                  # 
                  #xlim = c(-30, -10), ylim = c(5,15),

                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.4, "lines"),
                 # segment.color = "black",
                  show.legend = T)


