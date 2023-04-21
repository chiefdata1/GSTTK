rm(list = ls())
library(data.table)
load('20-gene.rda')

cna.region <- fread('CNV/HNSC-Gistic2.Level_4/all_lesions.conf_99.txt',data.table = F)
cna.gene <- fread('CNV/HNSC-Gistic2.Level_4/all_thresholded.by_genes.txt',data.table = F)
cna.region <- cna.region[,-c(3:9)]
#有个基因名字未匹配上，genecard找到别名
cna.gene$`Gene Symbol`[cna.gene$`Gene Symbol`=='MKI67IP'] <- 'NIFK'
cna.gene <- cna.gene[cna.gene$`Gene Symbol` %in% rownames(data_mut),]
cna.gene <- cna.gene[,-c(2,3)]
rownames(cna.gene) <- NULL
cna.gene <- tibble::column_to_rownames(cna.gene,'Gene Symbol')
cna.gene$Gain <- apply(cna.gene,1,function(x){sum(x>0)})
cna.gene$Loss <- apply(cna.gene,1,function(x){sum(x<0)})
cna.20 <- cna.gene[,109:110]
cna.20 <- tibble::rownames_to_column(cna.20,'Gene')

library(ggplot2)
library(reshape2)
darkblue <- "paleturquoise3"
lightblue <- "rosybrown3"
dd1 <- reshape(cna.20,
               varying = c("Gain","Loss"),
               v.names = "CNV",timevar = "GL",
               direction = "long")
dd1 <- dd1[,-4]
str(dd1)
dd1$GL <- ifelse(dd1$GL==1,"Gain","Loss")
dd1$GL <- factor(dd1$GL)
dd1$Gene <- factor(dd1$Gene,levels = gg)

ggplot(dd1,aes(CNV,Gene))+
  geom_segment(aes(yend=Gene),xend=0,color="black",
  size=0.6,lineend = "butt",alpha=1)+
  geom_point(aes(color=GL),size=4)+
  scale_color_brewer(palette = "Set2",direction = -1)+
  labs(x="CNA number",y="")+
  theme_classic()+
  # coord_flip()+
  scale_fill_manual(values = c("#FDBF6F", "#1F78B4")) +
  scale_color_manual(values = c("#FDBF6F", "#1F78B4")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1 ,colour = 'black'))+
  theme(legend.position = c("right"),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold",size = 12),
        axis.text = element_text(colour = "black",size = 12))
library(export)
graph2pdf(file='Figure/CNA-lollipop.pdf',width=5,height=7)



# 染色体圈图-------------------------------------------------------------------------
rm(list = ls())
library(data.table)
library(RCircos)
library(magrittr)
library(tidyverse)
library(rtracklayer)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
load('1-GRCh38_v39.gtf.rda')
load('20-gene.rda')
gtf_data <- gtf_data[gtf_data$gene_name%in%rownames(data_mut),]
gtf_data <- gtf_data[gtf_data$type=='gene',c(1,2,3,12)]

# 加载基因所在的位置和数值
gene_pos <- gtf_data
rm(gtf_data)
# 加载染色体Ideogram
(data("UCSC.HG38.Human.CytoBandIdeogram"))
pdf(file="circGene.pdf", height=5, width=5)

# 根据hg38构建染色体位置，只保留chr1-22,X,Y，在圈内部构建三圈轨道
RCircos.Set.Core.Components(UCSC.HG38.Human.CytoBandIdeogram,
                            chr.exclude = NULL,
                            tracks.inside = 3,
                            tracks.outside = 0)
RCircos.Set.Plot.Area()

# 绘制染色体
RCircos.Chromosome.Ideogram.Plot()

# 在第一圈用散点在基因所在的位置标注数值
# 调整配色
params <- RCircos.Get.Plot.Parameters()
params$track.background <- "grey" # 第三圈默认配色为wheat，模仿原文修改为灰色
RCircos.Reset.Plot.Parameters(params)

RCircos.Scatter.Plot(gene_pos,
                     # data.col = 5, # 用第5列的数值作为点的纵坐标
                     by.fold = 1, # 点的颜色cutoff，大于等于1的基因显示为红色点，小于等于-1的显示为蓝色点，-1到1之间为黑点
                     track.num = 1,
                     side = "in")

# 在第二圈绘制线段标注基因所在的位置
RCircos.Gene.Connector.Plot(genomic.data = gene_pos,
                            track.num = 2,
                            side = "in")

# 在第三圈标注基因名
RCircos.Gene.Name.Plot(gene_pos,
                       name.col = 4,
                       track.num = 3,
                       side = "in")
library(export)
graph2pdf(file ='Figure/chome-plot.pdf',width =10,height=10)

dev.off()

















