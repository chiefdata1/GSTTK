# -------------------------------------------------------------------------
rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(survival)
library(survminer)
library(dplyr)
load('HNSCC-tcga.rda')
gene <- readxl::read_xlsx('ENSEMBL.xlsx')
gene <- na.omit(gene)
gene_expr <- tpm[rownames(tpm)%in%gene$ENSEMBL,]

gene_expr <- as.data.frame(t(gene_expr))
tmp <- merge(clin[,1:3],gene_expr,by.x=1,by.y=0)
library(stringr)
tmp_t <- tmp[str_sub(tmp$ID,14,16)=='01A',]
outTab=data.frame()
for(i in colnames(tmp_t[,4:ncol(tmp_t)])){
  cox <- coxph(Surv(OS.time,OS) ~ get(i), data = tmp_t)
  coxSummary = summary(cox)
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}

UNI <- outTab$id[outTab$pvalue<0.05]
rm(cox,coxSummary,gene_expr,tmp,i)
samples <- c(colnames(tpm)[str_sub(colnames(tpm),14,16)=='01A'],
             colnames(tpm)[str_sub(colnames(tpm),14,16)=='11A'])
table(str_sub(colnames(tpm),14,16))#494肿瘤,44正常
# 差异分析 --------------------------------------------------------------------

count <- count[,samples]
tpm <- tpm[,samples]
group <- rep(c('Tumor','Normal'),c(494,44))
count_GSK <- count[rownames(count)%in%gene$ENSEMBL,]
tpm_GSK <- count[rownames(tpm)%in%gene$ENSEMBL,]


#表达值差异分析（三种方法）
if(T){
  library(limma)
  library(DESeq2)
  library(edgeR)

  group <- factor(group, levels = c("Normal","Tumor"),labels =  c("Normal","Tumor")) #limma中normal在前
  table(group)
  #group
  #2   1
  #130 226
  keep <- rowSums(count_GSK>=1) >= ncol(count_GSK)*0.5
  table(keep)
  cc <- count_GSK[keep,]
  #-count--DEseq2----------------------------------------------------------
  colData <- data.frame(group)
  dds <- DESeqDataSetFromMatrix(round(cc), colData, design= ~group)
  dds <- DESeq(dds)
  res<- results(dds,contrast=c("group","Tumor","Normal"),independentFiltering=FALSE)
  table(res$padj<0.05)
  dediff <- as.data.frame(res)%>%na.omit()
  DEseq2 <- dediff[order(dediff$log2FoldChange,decreasing = T),]

  #--count--edgeR------------------------------------
  design <- model.matrix(~ group)  #注意这里与limma包不一样
  colnames(design) <- levels(group)
  rownames(design) <- colnames(count_GSK)

  y <- DGEList(counts=cc,group= group)
  y <- calcNormFactors(y)

  y <- estimateDisp(y, design, robust=TRUE)  #是下面三个函数的组合
  # y <- estimateGLMCommonDisp(y, design)
  # y <- estimateGLMTrendedDisp(y, design)
  # y <- estimateGLMTagwiseDisp(y, design)


  fit <- glmQLFit(y, design, robust=TRUE)
  lrt <- glmQLFTest(fit)

  ordered_tags <- topTags(lrt, n=100000)
  allDiff=ordered_tags$table
  allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
  edgeR <- allDiff
  #--tpm--limma-------------------------------------------------------------
  design <- model.matrix(~ group)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(tpm_GSK)

  fit <- lmFit(tpm_GSK,design = design)
  fit2 <- eBayes(fit)
  limma <- topTable(fit2,adjust='fdr',coef=2,number=Inf)
  #----------------------------
  save(outTab,limma,DEseq2,edgeR,file = 'Diff_results.rda')

}


#聚类分析-------------------------------------------------------------------------
rm(list = ls())
library(survminer)
library(survival)
library(tidyverse)
library(stringr)
library(ConsensusClusterPlus)
load('Diff_results.rda')
load('HNSCC-tcga.rda')
rm(limma)
uni_name <- outTab$id[outTab$pvalue<0.05]
edger_name <- rownames(edgeR)[edgeR$PValue<0.05]
deseq2_name <- rownames(DEseq2)[DEseq2$pvalue<0.05]
tpm_t <- tpm[,colnames(tpm)[str_sub(colnames(tpm),14,15)=='01']]
data <- tpm_t[intersect(uni_name,edger_name),]
data1 <- tpm_t[intersect(uni_name,deseq2_name),]
data2 <- data[6:36,]

dir.create('ConsensusCluster/')
results = ConsensusClusterPlus(d = as.matrix(data2),
                               maxK=10,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               title='ConsensusCluster/',
                               innerLinkage="complete",
                               finalLinkage="complete",
                               clusterAlg="pam",
                               distance="euclidean",
                               seed=123456,
                               plot="png")
icl <- calcICL(results,title = 'ConsensusCluster/',plot = 'png')

## 保存分型信息
clusterNum=2
cluster=results[[clusterNum]][["consensusClass"]]

sub <- data.frame(Sample=names(cluster),Cluster=cluster)
sub$Cluster <- paste0('C',sub$Cluster)
table(sub$Cluster)

clin_edger <- merge(sub,clin,by=1)

fit <- survfit(Surv(OS.time,OS)~Cluster,clin_edger)
ggsurvplot(fit,pval = T,pval.method = T)

# 计算轮廓系数 ------------------------------------------------------------------
group <- substr(sub$Cluster,2,2)
nmf.input <- data2
library(factoextra)
library(cluster)
sil <-  silhouette(as.numeric(group), dist(t(nmf.input)))
plot(sil)
sildata <- sil[,1:3]%>%as.data.frame()
Cluster <- cbind(sub,sildata)
Cluster <- Cluster[Cluster$sil_width>0,]  ###最终434个样本
table(Cluster$Cluster) ##C1  C2---120  314
library(ggplot2)
x1 <- Cluster[Cluster$Cluster=='C1',]
x1 <- x1[order(x1$sil_width),]
x2 <- Cluster[Cluster$Cluster=='C2',]
x2 <- x2[order(x2$sil_width),]
c <- mean(x2$sil_width) ##0.21   0.163
mydata <- x1
mydata$ID <- 1:120

mydata <- x2
mydata$ID <- 1:314
#mydata$Name <- as.factor('C1')
#library(ggthemr)
#ggthemr('solarized')
ggplot(mydata,aes(ID,sil_width,fill=Cluster))+
  scale_fill_manual(values = c("#377EB8"))+
  geom_bar(stat = 'identity',width = 0.5)+
  geom_hline(aes(yintercept=mean(sil_width)),linetype=2,col="grey")+
  coord_flip() + theme_classic()+
  facet_grid(.~Cluster)+
  scale_x_discrete(expand = c(0.001,1.2))+
  #scale_y_discrete(expand = c(0.001,1))+
  xlab('')+
  ylab('Silhouette Width')+
  theme(#axis.line.y= element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(colour = 'black',#linetype="solid",
                                    size = 0.05, fill="seagreen"),
    strip.text = element_text(size = 12,face = "bold"),
    legend.position="none",
    #legend.title = element_text(size = 12,colour = 'black'),
    #legend.background =  element_rect(fill ='#FFCCCC',colour = "gray"),
    #legend.key.size = unit(0.15, "inches"),
    #legend.text = element_text(size = 11,colour = 'black'),
    axis.text = element_text(size = 12,colour = 'black',face = "bold"),
    axis.title = element_text(size = 14,colour = 'black'),
    plot.title = element_text(size=14,hjust=0.5)) +
    annotate("text", x = 8, y = 0.155, label = "Ave SW = 0.16",
             color="black",size = 5, fontface="bold",family="serif" )# 6*4


#绘制生存曲线
clin_sil <- merge(clin,Cluster[,1:2],by=1)
colnames(clin_sil)[20] <- 'group'
fit <- survfit(Surv(OS.time,OS)~group,clin_sil)
ggsurvplot(fit,pval = T,pval.method = T)

save(data,data2,uni_name,edger_name,clin_sil,clin_edger,file = 'cluster.rda')








