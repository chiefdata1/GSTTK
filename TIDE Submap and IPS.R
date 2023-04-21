# submap ------------------------------------------------------------------
rm(list = ls())
# 自定义函数用来产生submap需要的数据格式
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct,
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)

  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

# 创建submap需要的数据格式 (SKCM)

skcm.immunotherapy.logNC <- read.table("immune/submap/skcm.immunotherapy.47samples.log2CountsNorm.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) #原文提供的log2转化的标准化count值
skcm.immunotherapy.info <- read.table("immune/submap/skcm.immunotherapy.47sampleInfo.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R

load('TCGA-symbol.rda')

# 创建submap需要的数据格式
# submap不允许出现flat value, 因此最好选取过滤掉低表达的表达谱，这里使用的数据过滤了超过90%样本表达值均<1的基因
my_data <- my_data[apply(my_data,1,function(x)sum(abs(x)<1)<0.9*ncol(my_data)),]

tmp <- my_data
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) # 取交集后的基因列表

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# 产生输出数据的文件名
gct_file <- "immune/submap/skcm.immunotherapy.for.SubMap.gct"
cls_file <- "immune/submap/skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")
load('20-gene.rda')
# 提出亚型的样本，顺序排列
sample_C1 <- clin_edger$Sample[clin_edger$Cluster=='C1']
sample_C2 <- clin_edger$Sample[clin_edger$Cluster=='C2']

sam_info <- data.frame("ImmClust"=c(sample_C1,sample_C2),row.names = c(sample_C1,sample_C2))
sam_info$rank <- rep(c(1,2),times=c(length(sample_C1),length(sample_C2)))

# 产生输出数据的文件名
gct_file <- "immune/submap/Immune2.for.SubMap.gct"
cls_file <- "immune/submap/Immune2.for.SubMap.cls"
in_gct <- tmp[GENELIST,rownames(sam_info)]
# in_gct <- log2(tmp[GENELIST,rownames(sam_info)] + 1) # 产生和示例数据类似的形式，log2转化的标准化count值
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

library(scales)
library(RColorBrewer)
library(ggsci)
show_col(pal_igv()(8))




heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
heatmap.YlGnPe1 <- c('#802268FF','#466983FF','#BA6338FF','#5DB1DDFF','#F0E685FF')
library(scales)
show_col(heatmap.YlGnPe)
cherry    <- "#440259"
lightgrey <- "#dcddde"

# 绘图信息
# 把submap结果/130452/SubMap_SubMapResult.txt文件中的值填入相应的位置
# 输入文件中的名义p值和校正p值绘制热图
tmp <- matrix(c(1, 1, 1, 0.007992008,
                1, 1, 1, 1,
                0.8411588,0.2097902,0.8921079,0.000999001,
                0.5194805,0.2887113,0.4225774,0.944055944), # Bonferroni校正p值
              nrow = 4,byrow = T,dimnames = list(c("C1_p","C2_p","C1_b","C2_b"),c("CTLA4-noR","CTLA4-R","PD1-noR","PD1-R")))
library(pheatmap)
pheatmap(tmp, cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe1[5:1],
         gaps_row = 2,
         annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Bonferroni corrected","Bonferroni corrected"),row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
         filename = "Figure/immune/Submap.pdf")


# TIDE --------------------------------------------------------------------
# dir.create('immune/TIDE')
rm(list = ls())
load('TCGA-symbol.rda')
load('20-gene.rda')
TIDE <- my_data
# 这里为了得到比较好的结果，采用two direction median centered
TIDE <- sweep(TIDE,2, apply(TIDE,2,median,na.rm=T))
TIDE <- sweep(TIDE,1, apply(TIDE,1,median,na.rm=T))
write.table(TIDE,"immune/TIDE/TIDE_input.self_subtract",sep = "\t",row.names = T,col.names = NA,quote = F)

###在TIDE网页操作
result <- read.csv('immune/TIDE/TIDE_output.csv',h=T)
table(result$Responder)
result <- merge(result[,c(1,3)],clin_edger[,c(1,2)],by=1)


data <- data.frame('C1'=c(176,112),'C2'=c(151,52),row.names = c('False','True'))
fisher.test(data)
data2 <- data.frame('C1'=c(176/288,112/288),'C2'=c(151/203,52/203),row.names = c('False','True'))
data2 <- round(data2,2)
library(tidyr)
data3 <- gather(data2,'Risk','Value')
data3$Result <- rep(c('False','True'),2)
colnames(data3) <- c('Cluster','Frequency','Response')
data3$Value <- paste0(data3$Frequency*100,'%')

ggplot(data3,aes(x = Response,y = Frequency,fill = Cluster))+
  geom_col(position = 'dodge')+
  geom_text(aes(label=Value),
            colour='black',size=2,
            vjust=1.5,position = position_dodge(1))+
  theme_bw()+
  scale_fill_manual(values = c("#FDBF6F", "#1F78B4"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid = element_blank())+
  # ylab('Human leukocyte antigen expression')+
  xlab('')
library(export)
graph2pdf(file='Figure/immune/TIDE-percentage.pdf',height=4,width=3.5)


# IPS ---------------------------------------------------------------------
rm(list = ls())
library(dplyr)
library(tibble)
library(IOBR)
load('TCGA-symbol.rda')
load('20-gene.rda')
ips<-deconvo_tme(eset = my_data, method = "ips", plot= FALSE)
data_ips <- merge(clin_edger[,c(1,2)],ips,by.x = 1,by.y = 1)%>%
  column_to_rownames('Sample')

str(data_ips)
data_ips$Cluster <- factor(data_ips$Cluster,levels = c("C1","C2"))

library(ggsci)
ggplot(data_ips,aes(Cluster,AZ_IPS,fill=Cluster))+
  geom_boxplot(outlier.colour = 'grey30',
               outlier.size = 0.6,
               outlier.alpha = 0.6,
               outlier.shape = NA)+
  scale_fill_manual(values = c("#FDBF6F", "#1F78B4"))+
  stat_compare_means(symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.2, 1),
                                        symbols = c("***", "**", "*", "ns")),
                     label = "p.signif",label.y = 2.8)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle = 30,hjust = 1),
        panel.grid = element_blank())+
  labs(y='IPS z-score',x=NULL)

ggsave(filename = 'Figure/immune/IPS.pdf',width = 2,height = 3.5)


# 新抗原负荷-------------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(ggplot2)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
nal <- fread('Mutation/TCGA_PCA.mc3.v0.2.8.CONTROLLED.filtered.sample_neoantigens_10062017.tsv',data.table = F)
table(nal$cancer_type)
nal <- nal[nal$cancer_type=='HNSC',]
nal <- data.frame(ID=nal$sample,NAL=nal$neoantigen_num)
range(nal$NAL)
load('20-gene.rda')
tcga <-clin_edger[,c(1,2)]
tcga$Sample <- substr(tcga$Sample,1,12)
input <- merge(tcga,nal,by=1)
input$NAL <- log2(input$NAL+1)
input$NAL <- scale(input$NAL)
input$NAL <- as.numeric(input$NAL)
range(input$NAL)
ggplot(input,aes(Cluster,NAL)) +
  geom_boxplot(aes(color=Cluster),outlier.colour = NA,size=1,fill=NA)+
  geom_jitter(aes(fill=Cluster,color=Cluster),width = 0.2,shape=21,size=1,alpha=0.7)+
  stat_compare_means(method = "t.test",label = 'p.signif',label.x = 1.45,label.y = 3,size=6)+
  theme_bw(base_rect_size = 1.5)+
  labs(x=NULL,y='Neoantigen Load (NAL)',title = '')+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5,size = 15),
        axis.title.y = element_text(size = 14,colour = 'black'),
        axis.text.x = element_text(size = 14,colour = 'black'),
        axis.text.y = element_text(size = 12,colour = 'black'),
        axis.ticks = element_line(size = 1),
        panel.grid = element_blank())+
  scale_color_manual(values = c("#FDBF6F", "#1F78B4"))+
  scale_fill_manual(values = c("#FDBF6F", "#1F78B4"))+
  ylim(-2,3)

ggsave(filename = 'Figure/immune/Neoantigen-Boxplot.pdf',width = 2.5,height = 4.5)









































