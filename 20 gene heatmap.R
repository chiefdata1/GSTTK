rm(list = ls())
load('20-gene.rda')
#自定义标准化函数
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}


data <- as.data.frame(t(data_mut))
data <- merge(clin_edger[,1:2],data,by.x = 1, by.y = 0)
data <- data[order(data$Cluster),]
group <- data$Cluster
names(group) <- data$Sample
head(group)
#绘制热图
rownames(data) <- NULL
indata <- as.data.frame(t(tibble::column_to_rownames(data[,-2],'Sample')))
plotdata <- standarize.fun(indata[,names(group)],halfwidth = 2) # 绘图数据标准化（z-score并将绝对值超过2的数值截断）

annCol.tcga <- data.frame(Cluster = as.character(group), # 构建样本注释
                          row.names = names(group),
                          stringsAsFactors = F)
annColors <- list("Cluster" = c("C1" = '#FDBF6F',"C2" = '#1F78B4')) # 构建颜色注释
library(pheatmap)
heatmap.BLYW <- c("#3F047D","white","#B05900")
# 绘制热图
plotdata <- plotdata[c(3,4,8,11:14,19,20,1,2,5:7,9,10,15:18),]
pheatmap(plotdata,
         color = NMF:::ccRamp(x = heatmap.BLYW,n=64), # 原文颜色模版
         annotation_col = annCol.tcga[colnames(plotdata),,drop = F],
         annotation_colors = annColors,
         cluster_cols = FALSE, # 样本不聚类，按照亚型顺序排列
         cluster_rows = F, # 行聚类
         border_color = NA,
         treeheight_row = 20, # 修改行树高
         cellheight = 10, # 修改每个单元的高度
         cellwidth = 0.6, # 修改每个单元的宽度
         show_rownames = TRUE, # 显示行名
         show_colnames = FALSE, # 不显示列名
         filename = "Figure/heatmap of 20 genes in tcga cluster.pdf")






















