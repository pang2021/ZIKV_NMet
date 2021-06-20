library(ggplot2)
library(xlsx)
library(xlsxjars)

df.cp <- read.xlsx("E:/R/PHH/foldchange_metabolomics_fig4d.xlsx",
                   sheetIndex = 1, startRow = 1, as.data.frame = T, header = T)

ggplot(df.cp, aes(x = log2FC, y = reorder(Metabolite, log2FC))) +
  #
  geom_point(shape = 21, size = 1.8, fill = 'steelblue4') +
  
  theme_bw() +
  # 
  theme(plot.margin = margin(2, 1.5, 2, 1, "cm")) +
  
  #
  labs(x="Log2(Fold change)", y="", title="") +
  
  # title
  theme(plot.title = element_text(vjust = 6, hjust = 0.5,  size=9, face="bold"))  +
  # 
  theme(axis.title.x =element_text(vjust = -2, hjust = 0.55, size=10, face="plain"),
        axis.title.y=element_text(vjust = 5, size=10, face="bold")) +
  # 
  theme(axis.text.y = element_text(size = 8,face="plain", colour = "black"),
        axis.text.x = element_text(size = 9,face="plain", colour = "black"))

