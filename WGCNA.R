rm(list = ls())
library(stringr)
library(WGCNA)
library(export)
load('HNSCC-tcga.rda')
expr <- tpm[,str_sub(colnames(tpm),14,15)=='01']
##使用绝对中位差选择，推荐使用绝对中位差
WGCNA_matrix <- t(expr[order(apply(expr,1,mad), decreasing = T)[1:5000],])#mad代表绝对中位差
# 判断是否缺失较多的样本和基因
datExpr0 <- WGCNA_matrix
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
# 是否需要过滤，TRUE表示不需要，FALSE表示需要
gsg$allOK
### 通过样本聚类识别离群样本，去除离群样本 ###
sampleTree <-  hclust(dist(datExpr0), method = "average");#使用hclust函数进行均值聚类
# 绘制样本聚类图确定离群样本
sizeGrWindow(30,9)
# pdf(file = "figures/Step01-sampleClustering.pdf", width = 30, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# 根据上图判断，需要截取的高度参数h
abline(h = 180, col = "red")
dev.off()

# 去除离群得聚类样本，cutHeight参数要与上述得h参数值一致
clust <- cutreeStatic(sampleTree, cutHeight = 180, minSize = 10)
table(clust)
# clust 1聚类中包含着我们想要得样本，将其提取出来
keepSamples <- (clust==1)#转换为逻辑性向量
datExpr <- datExpr0[keepSamples, ]
# 记录基因和样本数，方便后续可视化
nGenes <- ncol(datExpr)#基因数
nSamples <- nrow(datExpr)#样本数
save(datExpr,nGenes,nSamples,file = "WGCNA_input.rda")


#WGCNA-------------------------------------------------------------------------
####数据准备--
# 清空所有变量
rm(list = ls())
# 加载包
library(WGCNA)
library(export)
# 允许多线程运行
enableWGCNAThreads()
# 加载表达矩阵
load("WGCNA_input.rda")
### 选择软阈值 ###
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# 进行网络拓扑分析
sft <-  pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)#β=power,就是软阈值
# 可视化结果
par(mfrow = c(1,2))#一个画板上，画两个图，一行两列
cex1 = 0.9;
softpower <- sft$powerEstimate
# 无尺度网络阈值得选择
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",#x轴
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",#y轴
     main = paste("Scale independence"));#标题
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     cex=cex1,
     col="red");
# 用红线标出R^2的参考值
abline(h=0.87,col="red")
# 平均连接度
plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels=powers,
     cex=cex1,
     col="red")
graph2pdf(file='Figure/WGCNA/软阈值.pdf',width=10,height=5)
dev.off()

# 无尺度网络检验，验证构建的网络是否是无尺度网络
softpower <- sft$powerEstimate

ADJ <-  abs(cor(datExpr,use="p"))^softpower#相关性取绝对值再幂次
k <-  as.vector(apply(ADJ,2,sum,na.rm=T))#对ADJ的每一列取和，也就是频次
par(mfrow = c(1,2))

hist(k)#直方图
scaleFreePlot(k,main="Check scale free topology")
graph2pdf(file='Figure/WGCNA/直方图.pdf',width=7,height=5)
dev.off()

###分步法构建网络--
## 计算邻接矩阵
adjacency <- adjacency(datExpr,power=softpower)
## 计算TOM拓扑矩阵
TOM <- TOMsimilarity(adjacency)
## 计算相异度
dissTOM <- 1- TOM
#模块初步聚类分析
library(flashClust)
geneTree <- flashClust(as.dist(dissTOM),method="average")
#绘制层次聚类树
plot(geneTree,xlab="",sub="",main="Gene clustering on TOM-based",labels=FALSE,hang=0.04)
dev.off()

#构建初步基因模块
#设定基因模块中至少30个基因
minModuleSize=30
# 动态剪切树识别网络模块
dynamicMods  <-  cutreeDynamic(dendro = geneTree,#hclust函数的聚类结果
                               distM = dissTOM,#
                               deepSplit = 2,
                               pamRespectsDendro = FALSE,
                               minClusterSize = minModuleSize)#设定基因模块中至少30个基因
# 将标签转换为颜色
dynamicColors <-  labels2colors(dynamicMods)
table(dynamicColors)#看聚类到哪些模块，哪些颜色
plotDendroAndColors(dendro = geneTree,
                    colors = dynamicColors,
                    groupLabels = "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE,
                    main = "Gene dendrogram and module colors")
graph2pdf(file = 'Figure/WGCNA/合并前模块.pdf',width=10,height=4)
dev.off()

#前面用动态剪切树聚类了一些模块，现在要对这步结果进一步合并，合并相似度大于0.75的模块，降低网络的复杂度
#计算基因模块特征向量
MEList <-  moduleEigengenes(datExpr,colors=dynamicColors)#计算特征向量
MEs <-  MEList$eigengenes;#提取特征向量
MEDiss1 <-  1-cor(MEs);#计算相异度
METree1 <-  flashClust(as.dist(MEDiss1),method="average")#对相异度进行flashClust聚类
#设置特征向量相关系数大于0.75
MEDissThres <-  0.25;#相异度在0.25以下，也就是相似度大于0.75，对这些模块合并
#合并模块
merge <-  mergeCloseModules(datExpr, #合并相似度大于0.75的模块
                            dynamicColors,
                            cutHeight = MEDissThres,
                            verbose=3)
mergedColors <-  merge$colors

table(dynamicColors)#动态剪切树的模块颜色
table(mergedColors)#合并后的模块颜色，可以看到从18个模块变成了14个模块
mergedMEs <-  merge$newMEs#合并后的14个模块

#重新命名合并后的模块
moduleColors = mergedColors;
colorOrder = c("grey",standardColors(50));
moduleLabels = match(moduleColors,colorOrder)-1;
MEs = mergedMEs;
MEDiss2 = 1-cor(MEs);#计算相异度
METree2 = flashClust(as.dist(MEDiss2),method="average");#对合并后的模块进行聚类

#绘制聚类结果图
# pdf(file="figures/Step03-MECombined.pdf",width=12,height=5)
par(mfrow=c(1,2))
plot(METree1,xlab="",sub="",main="Clustering of ME before combined")# METree1是动态剪切树形成的模块
abline(h=MEDissThres,col="red")#相异度为0.25
plot(METree2,xlab="",sub="",main="Clustering of ME after combined")# METree2是合并后的模块
dev.off()

# pdf(file="figures/Step03-MergedDynamics.pdf",width=10,height=4)
#合并后的图
plotDendroAndColors(dendro = geneTree,#剪切树
                    colors = cbind(dynamicColors,mergedColors),#将两种方法形成的模块颜色合并在一起
                    groupLabels = c("Dynamic Tree Cut","Merged Dynamics"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide=T,
                    guideHang=0.05,
                    main="Cluster Dendrogram")
graph2pdf(file = 'Figure/WGCNA/MergedDynamics.pdf',width=10,height=4)
dev.off()
#合并后的模块图
plotDendroAndColors(dendro = geneTree,
                    colors = mergedColors,
                    groupLabels = "Merged Dynamics",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE,
                    main = "Cluster Dendrogram")
graph2pdf(file = 'Figure/WGCNA/hebingqian.pdf',width=8,height=4)
dev.off()


# 模块中基因数
write.csv(table(moduleColors),"MEgeneCount.csv",quote = F,row.names = F)

# 保存构建的网络信息
moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file="Step_by_step_buildnetwork.rda")





#### 绘制样本间的相关性 ----
load("WGCNA_input.Rda")
load(file="Step_by_step_buildnetwork.rda")


MEs = orderMEs(MEs)

sizeGrWindow(5,7.5);
# pdf(file = "figures/Step03-moduleCor.pdf", width = 5, height = 7.5);

par(cex = 0.9)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2),
                      cex.lab = 0.8, xLabelsAngl = 90)
graph2pdf(file = 'Figure/WGCNA/moduleCor.pdf',width = 5, height = 7.5)
dev.off()

#### TOMplot ----
dissTOM <-  1-TOMsimilarityFromExpr(datExpr, power = softpower);
nSelect <-  400
# 随机选取400个基因进行可视化，s设置seed值，保证结果的可重复性
set.seed(10);
select <-  sample(nGenes, size = nSelect);
selectTOM <-  dissTOM[select, select];
# 对选取的基因进行重新聚类
selectTree <-  hclust(as.dist(selectTOM), method = "average")
selectColors <-  mergedColors[select];
# 打开一个绘图窗口
sizeGrWindow(9,9)
# 美化图形的设置
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors,
        main = "Network heatmap plot, selected genes")
graph2pdf(file = "Figure/WGCNA/TOMplot.pdf", width = 9, height = 9)
dev.off()
#Cluster 基因-------------------------------------------------------------------------
####模块与性状之间的关系
# 清空环境变量
rm(list = ls())
# 加载包
library(WGCNA)
library(export)
# 加载表达矩阵
load("WGCNA_input.Rda")
load('clean-clin.rda')
load('20-gene.rda')
rm(tpm,data_mut)
dim(datExpr)
clin <- merge(clin,clin_edger[,1:2],by = 1)
table(clin$Cluster)
clin$Cluster <- ifelse(clin$Cluster=='C1',1,0)

# 读入临床信息
clinical <- clin[,c(1,10:19)]
rm(clin_edger,clin)
# 查看临床信息
head(clinical)
# 对表达矩阵进行预处理
clinical1 <- clinical[clinical$ID%in%rownames(datExpr),]#去掉三个离群值
rownames(clinical1) <- NULL
clinical1 <- tibble::column_to_rownames(as.data.frame(clinical1),'ID')
datTraits <- as.data.frame(do.call(cbind,lapply(clinical1, as.numeric)))
rownames(datTraits) <- rownames(clinical1)
rm(clinical,clinical1)
# 对样本进行聚类
sampleTree2 <- hclust(dist(datExpr), method = "average")

# 将临床信息转换为颜色，白色表示低，红色表示高，灰色表示缺失
traitColors <- numbers2colors(datTraits, signed = FALSE)

# 样本聚类图与样本性状热图
plotDendroAndColors(sampleTree2,
                    traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
####网络的分析----
#基因模块与临床信息的关系
# 加载构建的网络
load(file = "Step_by_step_buildnetwork.rda")
# 对模块特征矩阵进行排序
MEs <- orderMEs(MEs)
#计算模型特征矩阵和样本信息矩阵的相关度。
moduleTraitCor <- cor(MEs, datTraits, use="p")
# write.csv(moduleTraitCor,file="data/Step04-modPhysiological-cor.csv",quote=F)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
# write.csv(moduleTraitPvalue,file="data/Step04-modPhysiological.p.xls",quote=F)

#使用labeledHeatmap()将上述相关矩阵和p值可视化
textMatrix <- paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix) <- dim(moduleTraitCor)
moduleTraitCor <- edit(moduleTraitCor)
textMatrix <- edit(textMatrix)
# 基因模块与临床信息相关性图
labeledHeatmap(Matrix=moduleTraitCor,#模块和表型的相关性矩阵，这个参数最重要，其他可以不变
               xLabels=colnames(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=F,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=F,
               cex.text=0.6,
               cex.lab=0.6,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))
graph2pdf(file = "Figure/WGCNA/性状相关性.pdf", width = 9, height = 9)
dev.off()



# 提取模块中的基因 ----------------------------------------------------------------
gene_name <- rownames(as.data.frame(t(datExpr)))
head(gene_name)
gene_C1 <- gene_name[which(moduleColors=="brown")]
gene_C2 <- gene_name[which(moduleColors=="red")]

#单一模块与性状之间的相关性
C1 <- as.data.frame(datTraits$Cluster)
names(C1) <- "C1"
modNames <- substring(names(MEs), 3)#模块对应的颜色
# 计算基因模块特征
nSamples <- 494#样本数
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
# 对结果进行命名
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
# 计算C2基因特征显著性
geneTraitSignificance <- as.data.frame(cor(datExpr, C2, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
# 对结果进行命名
names(geneTraitSignificance) <- paste("GS.", names(C2), sep="")
names(GSPvalue) <- paste("p.GS.", names(C2), sep="")
# 设置需要分析的模块名称，此处为brown模块
module <- "pink"
# 提取yellow模块数据
column <- match(module, modNames)
moduleGenes <- moduleColors==module
par(mfrow = c(1,1))
col <- c("#FDBF6F", "#1F78B4")
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for C1",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = col[1])
graph2pdf(file='Figure/WGCNA/C2散点图.pdf',width=7,height=7)
dev.off()


save(gene_C2,gene_C1,file = 'NTP-gene.rda')














































