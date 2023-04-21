rm(list = ls())
library(GSVA)
load('TCGA-symbol.rda')
load('hallmark.gs.RData')
load('20-gene.rda')
# 这一句就完成了GSVA分析
gsva_es <- gsva(as.matrix(my_data), gs)
group_list <- clin_edger[,1:2]
colnames(group_list) <- c('sample','group')
group_list$group <- ifelse(group_list$group=="High-risk",'High','Low')

# 设置对比
design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design
library(limma)
# 构建差异比较矩阵
contrast.matrix <- makeContrasts(C1-C2, levels = design)

# 差异分析，High vs. Low
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

#把通路的limma分析结果保存到文件
write.csv(x, "gsva_limma.csv", quote = F)
library(stringr)
#输出t值，用做画bar的输入数据
pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)
write.csv(df, "easy_input2_bar.csv", quote = F, row.names = F)



# 画图 ----------------------------------------------------------------------
#开始画图
head(df)
range(df$score)
#按照score的值分组
cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))
df$ID <- gsub(pattern = '_',' ',df$ID)
df$score <- ifelse(df$score>15,11,df$score)
#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)
library(ggplot2)
library(ggsci)
ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_manual(values = c('#1F78B4','snow3','#FDBF6F'), guide = FALSE) +

  #画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff),
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细

  #写label
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = "inward" ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=ID, color = group),
            size = 3, hjust = "outward") +
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +

  xlab("") +ylab("t value of GSVA score\n C1 vs C2 group")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴
ggsave("Figure/GSEA/gsva.pdf", width = 6.2, height = 8)
