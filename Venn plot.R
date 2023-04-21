rm(list = ls())
load('Diff_results.rda')
rm(DEseq2,limma)
#####火山图
library(ggplot2)
#读取数据
dataset <- edgeR
# 设置pvalue和logFC的阈值
cut_off_pvalue = 0.05
cut_off_logFC = 1
dataset$change = ifelse(dataset$PValue < cut_off_pvalue & abs(dataset$logFC) >= cut_off_logFC,
                        ifelse(dataset$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
table(dataset$change)

ggplot(
  #设置数据
  dataset,
  aes(x = logFC,
      y = -log10(PValue),
      colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#1F78B4", "#d2dae2","#FDBF6F"))+

  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +

  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+

  # 图例
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank()
  )
library(ggsci)
library(export)
graph2pdf(file = 'Figure/Volcano Plot.pdf',width =5,height=5)

#单因素分析森林图-----------------
load('cluster.rda')
library(forestplot)

rs_forest <- outTab[outTab$id%in%rownames(data),]
# 读入数据的时候大家一定要把header设置成FALSE，保证第一行不被当作列名称〿
str(rs_forest)
rs_forest[,2:5] <- apply(rs_forest[,2:5],2,as.numeric)
rs_forest[,2:5] <- round(rs_forest[,2:5],digits = 3)

forestplot(labeltext = as.matrix(rs_forest[,c(1,2,5)]),
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           mean = rs_forest$HR, #设置均倿
           lower = rs_forest$HR.95L, #设置均值的lowlimits陿
           upper = rs_forest$HR.95H, #设置均值的uplimits陿
           #is.summary = c(T,T,T,T,T,T,T,T,T),
           #该参数接受一个逻辑向量，用于定义数据中的每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不昿0
           boxsize = 0.3, #设置点估计的方形大小
           lineheight = unit(8,'mm'),#设置图形中的行距
           colgap = unit(2,'mm'),#设置图形中的列间跿
           lwd.zero = 2,#设置参考线的粗绿
           lwd.ci = 2,#设置区间估计线的粗细
           col=fpColors(box='#458B00', summary= "#8B008B",lines = 'black',zero = '#7AC5CD'),
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           xlab="The estimates",#设置x轴标筿
           lwd.xaxis=2,#设置X轴线的粗绿
           lty.ci = "solid",
           graph.pos = 3)#设置森林图的位置，此处设置为4，则出现在第四列
graph2pdf(file = 'Figure/uni-cox.pdf',width=8,height = 14)

dev.off()


#韦恩图
library(grid)
library(futile.logger)
library(VennDiagram)
color <- c("#FDBF6F", "#1F78B4")
venn.plot <- draw.pairwise.venn(
  area1 = 102,
  area2 = 53,
  cross.area = 36,
  category = c('DEGs','Unicox Gene'),
  fill = color,
  lty = "blank",
  cex = c(1,1,1));

graph2pdf(file = 'Figure/Veen.pdf',width = 5, height =5)














