# 差异分析 --------------------------------------------------------------------
rm(list = ls())
library(clusterProfiler)
library(rvcheck)
library(enrichplot)
library(edgeR)
library(DESeq2)
library(dplyr)
library(ggsci)
library(ggplot2)
library(org.Hs.eg.db)
load('HNSCC-tcga.rda')
load('20-gene.rda')

gene <- bitr(geneID = rownames(count),fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
count <- merge(gene,count,by.x = 1,by.y=0)[,-1]
count <- count[!duplicated(count$SYMBOL),]
rownames(count) <- NULL
count <- tibble::column_to_rownames(count,'SYMBOL')

count <-  as.data.frame(t(count))
head(count[1:3,1:3])
data <- merge(clin_edger[,1:2],count,by.x=1,by.y=0)
group <- data$Cluster
data <- tibble::column_to_rownames(data,'Sample')
data <- as.data.frame(t(data[,-1]))
keep <- rowSums(data>=1) >= ncol(data)*0.5
table(keep)
cc <- data[keep,]
colData <- data.frame(group)
dds <- DESeqDataSetFromMatrix(round(cc), colData, design= ~group)
dds <- DESeq(dds)
res<- results(dds,contrast=c("group","C1","C2"),independentFiltering=FALSE)
table(res$padj<0.05)
dediff <- as.data.frame(res)%>%na.omit()
DEseq2 <- dediff[order(dediff$log2FoldChange,decreasing = T),]

save(DEseq2,file = 'Cluster-risk-diff.rda')



df.id<- bitr(rownames(dediff), fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
easy.df<-merge(df.id,dediff,by.x=1,by.y=0)
sortdf<-easy.df[order(easy.df$log2FoldChange, decreasing = T),]
#sortdf$logFC <- 2^sortdf$logFC
gene.expr <- sortdf$log2FoldChange
names(gene.expr) <- sortdf$ENTREZID
head(gene.expr)

kk <- gseKEGG(gene.expr,minGSSize = 1,maxGSSize = 10000,
              pvalueCutoff = 1,pAdjustMethod = 'BH')
sortkk<-kk[order(kk$NES, decreasing = T),]

gg <- gseGO(geneList = gene.expr,ont = 'BP',OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',
            pvalueCutoff = 0.05,pAdjustMethod = 'BH',minGSSize = 1,maxGSSize = 10000)
sortgg<- gg[order(gg$NES, decreasing = T),]

sortkk<-kk[order(kk$NES, decreasing = T),]
sortgg<- gg[order(gg$NES, decreasing = T),]
save(sortkk,sortgg,file = 'Cluster-GSEA.rda')

library(ggsci)
library(scales)
show_col(mycol)
mycol <- c(pal_npg()(10),pal_jama()(7)[c(5:7,1:4)])
#C1
library(enrichplot)
i1 <- c('hsa04612','hsa04672','hsa04658','hsa04659','hsa04660')
d1 <- sortkk[sortkk$ID%in%i1,]
gseaplot2(kk,geneSetID = i1,
          color = mycol[c(1:5)],
          title = 'Specific KEGG pathways in C1 group',
          rel_heights = c(1.3, 0.3, 0.6))
# dir.create('Figure/GSEA')
ggsave(filename = 'Figure/GSEA/GSEA-KEGG-C1.pdf',width = 7,height = 5)
dev.off()

#KEGG-C1另一种作业里的画法-----------------------------------------------------
geneSetID <- c('hsa04612','hsa04672','hsa04658','hsa04659','hsa04660')
# 突出显示感兴趣的基因
selectedGeneID <- c(rownames(data_mut))
# 自定义足够多的颜色
mycol <- c(pal_npg()(10),pal_jama()(7)[c(5:7,1:4)])[1:5]
show_col(mycol)
# 字号
base_size = 12

# rank mode
rankmode <- "comb" #合起来画
#rankmode <- "sep" #如果通路多，分开画更好


## DIY多条通路
x <- kk
geneList <- position <- NULL ## to satisfy codetool

#合并多条通路的数据
gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
gsdata$gsym <- rep(sortdf$SYMBOL,5)#5是geneSetID数量

# 画running score
p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
  geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
  scale_color_manual(values = mycol) +

  scale_x_continuous(expand=c(0,0)) + #两侧不留空
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
  ylab("Enrichment\n Score") +

  theme_bw(base_size) +
  theme(panel.grid = element_blank()) + #不画网格

  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +

  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
p.res

# 画rank
if (rankmode == "comb") {
  rel_heights <- c(1.3, 0.3, 0.6) # 上中下三个部分的比例

  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_bar(position = "dodge",stat = "identity",
             aes(x, y = position*(-0.1), fill=Description))+
    xlab(NULL) + ylab(NULL) +
    scale_fill_manual(values = mycol) + #用自定义的颜色

    theme_bw(base_size) +
    theme(panel.grid = element_blank()) + #不画网格

    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))

} else if (rankmode == "sep") {
  rel_heights <- c(1.3, 0.3, 0.6) # 上中下三个部分的比例

  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }

  head(gsdata)

  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) +
    scale_color_manual(values = mycol) + #用自定义的颜色

    theme_bw(base_size) +
    theme(panel.grid = element_blank()) + #不画网格

    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
} else {stop("Unsupport mode")}

p2

# 画变化倍数
df2 <- p.res$data
df2$y <- p.res$data$geneList[df2$x]
df2$gsym <- p.res$data$gsym[df2$x]
head(df2)

# 提取感兴趣的基因的变化倍数
selectgenes <- data.frame(gsym = selectedGeneID)
selectgenes <- merge(selectgenes, df2, by = "gsym")
selectgenes <- selectgenes[selectgenes$position == 0,]
head(selectgenes)
selectgenes <- selectgenes[!duplicated(selectgenes$gsym),]
library(ggrepel)
p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) +
  geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
               color = "grey") +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = mycol, guide=FALSE) + #用自定义的颜色
  scale_color_manual(values = mycol, guide=FALSE) + #用自定义的颜色

  scale_x_continuous(expand=c(0,0)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
  ylab("Ranked list\n metric") +
  xlab("Rank in ordered dataset") +

  theme_bw(base_size) +
  theme(panel.grid = element_blank()) +

  # 显示感兴趣的基因的基因名
  geom_text_repel(data = selectgenes,
                  show.legend = FALSE, #不显示图例
                  direction = "x", #基因名横向排列在x轴方向
                  ylim = c(2, NA), #基因名画在-2下方
                  angle = 90, #基因名竖着写
                  size = 2.5, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))

p.pos

library(cowplot)
# 组图
plotlist <- list(p.res, p2, p.pos)
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text())
plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)

library(export)
dir.create('Figure/GSEA')
graph2pdf(file='Figure/GSEA/KEGG-C1.pdf',width=10,height=7)


#KEGG-C2另一种作业里的画法-----------------------------------------------------
geneSetID <- c('hsa03010','hsa00190','hsa00970','hsa01232','hsa03013')
# 突出显示感兴趣的基因
selectedGeneID <- c(rownames(data_mut))
# 自定义足够多的颜色
mycol <- c(pal_npg()(10),pal_jama()(7)[c(5:7,1:4)])[1:5]
show_col(mycol)
# 字号
base_size = 12

# rank mode
rankmode <- "comb" #合起来画
#rankmode <- "sep" #如果通路多，分开画更好


## DIY多条通路
x <- kk
geneList <- position <- NULL ## to satisfy codetool

#合并多条通路的数据
gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
gsdata$gsym <- rep(sortdf$SYMBOL,5)#5是geneSetID数量

# 画running score
p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
  geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
  scale_color_manual(values = mycol) +

  scale_x_continuous(expand=c(0,0)) + #两侧不留空
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
  ylab("Enrichment\n Score") +

  theme_bw(base_size) +
  theme(panel.grid = element_blank()) + #不画网格

  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +

  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
p.res

# 画rank
if (rankmode == "comb") {
  rel_heights <- c(1.3, 0.3, 0.6) # 上中下三个部分的比例

  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_bar(position = "dodge",stat = "identity",
             aes(x, y = position*(-0.1), fill=Description))+
    xlab(NULL) + ylab(NULL) +
    scale_fill_manual(values = mycol) + #用自定义的颜色

    theme_bw(base_size) +
    theme(panel.grid = element_blank()) + #不画网格

    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))

} else if (rankmode == "sep") {
  rel_heights <- c(1.3, 0.3, 0.6) # 上中下三个部分的比例

  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }

  head(gsdata)

  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) +
    scale_color_manual(values = mycol) + #用自定义的颜色

    theme_bw(base_size) +
    theme(panel.grid = element_blank()) + #不画网格

    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
} else {stop("Unsupport mode")}

p2

# 画变化倍数
df2 <- p.res$data
df2$y <- p.res$data$geneList[df2$x]
df2$gsym <- p.res$data$gsym[df2$x]
head(df2)

# 提取感兴趣的基因的变化倍数
selectgenes <- data.frame(gsym = selectedGeneID)
selectgenes <- merge(selectgenes, df2, by = "gsym")
selectgenes <- selectgenes[selectgenes$position == 0,]
head(selectgenes)
selectgenes <- selectgenes[!duplicated(selectgenes$gsym),]
library(ggrepel)
p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) +
  geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
               color = "grey") +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = mycol, guide=FALSE) + #用自定义的颜色
  scale_color_manual(values = mycol, guide=FALSE) + #用自定义的颜色

  scale_x_continuous(expand=c(0,0)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
  ylab("Ranked list\n metric") +
  xlab("Rank in ordered dataset") +

  theme_bw(base_size) +
  theme(panel.grid = element_blank()) +

  # 显示感兴趣的基因的基因名
  geom_text_repel(data = selectgenes,
                  show.legend = FALSE, #不显示图例
                  direction = "x", #基因名横向排列在x轴方向
                  ylim = c(2, NA), #基因名画在-2下方
                  angle = 90, #基因名竖着写
                  size = 2.5, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))

p.pos

library(cowplot)
# 组图
plotlist <- list(p.res, p2, p.pos)
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text())
plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)

library(export)
graph2pdf(file='Figure/GSEA/KEGG-C2.pdf',width=10,height=7)

gseaplot2(kk,geneSetID = geneSetID,
          color = mycol[c(1:5)],
          title = 'Specific KEGG pathways in C1 group',
          rel_heights = c(1.3, 0.3, 0.6))



#GO-C1另一种作业里的画法-----------------------------------------------------
geneSetID <- c('GO:0006958','GO:0002250','GO:0002377','GO:0050853','GO:0019724')
# 突出显示感兴趣的基因
selectedGeneID <- c(rownames(data_mut))
# 自定义足够多的颜色
mycol <- c(pal_npg()(10),pal_jama()(7)[c(5:7,1:4)])[1:5]
show_col(mycol)
# 字号
base_size = 12

# rank mode
rankmode <- "comb" #合起来画
#rankmode <- "sep" #如果通路多，分开画更好


## DIY多条通路
x <- gg
geneList <- position <- NULL ## to satisfy codetool

#合并多条通路的数据
gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
gsdata$gsym <- rep(sortdf$SYMBOL,5)#5是geneSetID数量

# 画running score
p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
  geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
  scale_color_manual(values = mycol) +

  scale_x_continuous(expand=c(0,0)) + #两侧不留空
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
  ylab("Enrichment\n Score") +

  theme_bw(base_size) +
  theme(panel.grid = element_blank()) + #不画网格

  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +

  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
p.res

# 画rank
if (rankmode == "comb") {
  rel_heights <- c(1.3, 0.3, 0.6) # 上中下三个部分的比例

  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_bar(position = "dodge",stat = "identity",
             aes(x, y = position*(-0.1), fill=Description))+
    xlab(NULL) + ylab(NULL) +
    scale_fill_manual(values = mycol) + #用自定义的颜色

    theme_bw(base_size) +
    theme(panel.grid = element_blank()) + #不画网格

    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))

} else if (rankmode == "sep") {
  rel_heights <- c(1.3, 0.3, 0.6) # 上中下三个部分的比例

  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }

  head(gsdata)

  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) +
    scale_color_manual(values = mycol) + #用自定义的颜色

    theme_bw(base_size) +
    theme(panel.grid = element_blank()) + #不画网格

    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
} else {stop("Unsupport mode")}

p2

# 画变化倍数
df2 <- p.res$data
df2$y <- p.res$data$geneList[df2$x]
df2$gsym <- p.res$data$gsym[df2$x]
head(df2)

# 提取感兴趣的基因的变化倍数
selectgenes <- data.frame(gsym = selectedGeneID)
selectgenes <- merge(selectgenes, df2, by = "gsym")
selectgenes <- selectgenes[selectgenes$position == 0,]
head(selectgenes)
selectgenes <- selectgenes[!duplicated(selectgenes$gsym),]
library(ggrepel)
p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) +
  geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
               color = "grey") +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = mycol, guide=FALSE) + #用自定义的颜色
  scale_color_manual(values = mycol, guide=FALSE) + #用自定义的颜色

  scale_x_continuous(expand=c(0,0)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
  ylab("Ranked list\n metric") +
  xlab("Rank in ordered dataset") +

  theme_bw(base_size) +
  theme(panel.grid = element_blank()) +

  # 显示感兴趣的基因的基因名
  geom_text_repel(data = selectgenes,
                  show.legend = FALSE, #不显示图例
                  direction = "x", #基因名横向排列在x轴方向
                  ylim = c(2, NA), #基因名画在-2下方
                  angle = 90, #基因名竖着写
                  size = 2.5, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))

p.pos

library(cowplot)
# 组图
plotlist <- list(p.res, p2, p.pos)
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text())
plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)

library(export)
graph2pdf(file='Figure/GSEA/GO-C1.pdf',width=10,height=7)

gseaplot2(gg,geneSetID = geneSetID,
          color = mycol[c(1:5)],
          title = 'Specific GO pathways in C1 group',
          rel_heights = c(1.3, 0.3, 0.6))


#GO-C2另一种作业里的画法-----------------------------------------------------
geneSetID <- c('GO:0042773','GO:0006119','GO:0019646','GO:0033108','GO:0042255')
# 突出显示感兴趣的基因
selectedGeneID <- c(rownames(data_mut))
# 自定义足够多的颜色
mycol <- c(pal_npg()(10),pal_jama()(7)[c(5:7,1:4)])[1:5]
show_col(mycol)
# 字号
base_size = 12

# rank mode
rankmode <- "comb" #合起来画
#rankmode <- "sep" #如果通路多，分开画更好


## DIY多条通路
x <- gg
geneList <- position <- NULL ## to satisfy codetool

#合并多条通路的数据
gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
gsdata$gsym <- rep(sortdf$SYMBOL,5)#5是geneSetID数量

# 画running score
p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
  geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
  scale_color_manual(values = mycol) +

  scale_x_continuous(expand=c(0,0)) + #两侧不留空
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
  ylab("Enrichment\n Score") +

  theme_bw(base_size) +
  theme(panel.grid = element_blank()) + #不画网格

  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +

  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
p.res

# 画rank
if (rankmode == "comb") {
  rel_heights <- c(1.3, 0.3, 0.6) # 上中下三个部分的比例

  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_bar(position = "dodge",stat = "identity",
             aes(x, y = position*(-0.1), fill=Description))+
    xlab(NULL) + ylab(NULL) +
    scale_fill_manual(values = mycol) + #用自定义的颜色

    theme_bw(base_size) +
    theme(panel.grid = element_blank()) + #不画网格

    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))

} else if (rankmode == "sep") {
  rel_heights <- c(1.3, 0.3, 0.6) # 上中下三个部分的比例

  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }

  head(gsdata)

  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) +
    scale_color_manual(values = mycol) + #用自定义的颜色

    theme_bw(base_size) +
    theme(panel.grid = element_blank()) + #不画网格

    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
} else {stop("Unsupport mode")}

p2

# 画变化倍数
df2 <- p.res$data
df2$y <- p.res$data$geneList[df2$x]
df2$gsym <- p.res$data$gsym[df2$x]
head(df2)

# 提取感兴趣的基因的变化倍数
selectgenes <- data.frame(gsym = selectedGeneID)
selectgenes <- merge(selectgenes, df2, by = "gsym")
selectgenes <- selectgenes[selectgenes$position == 0,]
head(selectgenes)
selectgenes <- selectgenes[!duplicated(selectgenes$gsym),]
library(ggrepel)
p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) +
  geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
               color = "grey") +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = mycol, guide=FALSE) + #用自定义的颜色
  scale_color_manual(values = mycol, guide=FALSE) + #用自定义的颜色

  scale_x_continuous(expand=c(0,0)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
  ylab("Ranked list\n metric") +
  xlab("Rank in ordered dataset") +

  theme_bw(base_size) +
  theme(panel.grid = element_blank()) +

  # 显示感兴趣的基因的基因名
  geom_text_repel(data = selectgenes,
                  show.legend = FALSE, #不显示图例
                  direction = "x", #基因名横向排列在x轴方向
                  ylim = c(2, NA), #基因名画在-2下方
                  angle = 90, #基因名竖着写
                  size = 2.5, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))

p.pos

library(cowplot)
# 组图
plotlist <- list(p.res, p2, p.pos)
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text())
plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)

library(export)
graph2pdf(file='Figure/GSEA/GO-C2.pdf',width=10,height=7)

gseaplot2(gg,geneSetID = geneSetID,
          color = mycol[c(1:5)],
          title = 'Specific GO pathways in C2 group',
          rel_heights = c(1.3, 0.3, 0.6))

dev.off()















































