#ssGSEA-----------------------------------------------------------------------
rm(list = ls())
library(survival)
library(ggkm)
library(ggplot2)
library(survminer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tibble)
library(pheatmap)
library(stringr)
library(dplyr)
library(GSVA)
immune <- read.table("pan_cancer_immune.txt", sep="\t",
                     header=T, check.names=F)

CellType <- as.vector(unique(immune$`Cell type`))

# 生成免疫细胞对应的基因标签的列表，用于ssGSEA分析
geneset <- lapply(CellType, function(x){
  x <- as.vector(immune[immune$`Cell type`==x,1])
  return(x)
})
names(geneset)<- CellType
#读入ssGSEA分析用的数据
load('HNSCC-tcga.rda')
load('20-gene.rda')
rm(count,clin,data_mut)
tpm <- tpm[,colnames(tpm)%in%clin_edger$Sample]
my_data <- tpm

class(my_data[1,1])
library(clusterProfiler)
library(org.Hs.eg.db)
Symbol<- bitr(rownames(my_data), fromType = "ENSEMBL", toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
my_data <- merge(Symbol,my_data,by.x=1,by.y=0)[,-1]

Symbol <- unique(my_data$SYMBOL)
my_data <- my_data%>%distinct(SYMBOL,.keep_all = T)

my_data <- tibble::column_to_rownames(my_data,'SYMBOL')
save(my_data,file = 'TCGA-symbol.rda')

# 进行ssGSEA分析
gsva_matrix<- gsva(as.matrix(my_data),
                   geneset,
                   method='ssgsea',#设置成ssgsea
                   kcdf='Gaussian',#高斯分布
                   abs.ranking=TRUE)

ssgseaScore <- t(scale(t(gsva_matrix)))
save(ssgseaScore,file = "data_ssgsea.rda")

# 绘制箱线图 -------------------------------------------------------------------
rm(list = ls())
library(dplyr)
library(export)
load("20-gene.rda")
load('data_ssgsea.rda')
ssGSEA <- as.data.frame(ssgseaScore)
ssGSEA <- as.data.frame(t(ssGSEA))
C1 <- ssGSEA[rownames(ssGSEA)%in%clin_edger$Sample[clin_edger$Cluster=='C1'],]
C2 <- ssGSEA[rownames(ssGSEA)%in%clin_edger$Sample[clin_edger$Cluster=='C2'],]
C1 <- as.data.frame(t(C1))
C1$sum <- apply(C1,1,mean)
C1$Group <- rep('C1',28)
C1 <- C1%>%
  dplyr::select(Group,sum)%>%
  tibble::rownames_to_column('Cell')
C2 <- as.data.frame(t(C2))
C2$sum <- apply(C2,1,mean)
C2$Group <- rep('C2',28)
C2 <- C2%>%
  dplyr::select(Group,sum)%>%
  tibble::rownames_to_column('Cell')

data <- rbind(C1,C2)
data$sum <- round(data$sum,2)
ggplot(data,aes(x = Cell,y = sum,fill = Group))+
  geom_col(position = 'dodge')+
  geom_text(aes(label=sum),
            colour='black',size=2,
            vjust=1.5,position = position_dodge(0.9))+
  theme_bw()+
  scale_fill_manual(values = c("#FDBF6F", "#1F78B4"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid = element_blank())+
  ylab('Mean expression level of immune cells')+
  xlab('')
dir.create('Figure/immune')
graph2pdf(file = 'Figure/immune/ssgsea.pdf',width = 10, height = 7)


# 免疫检查点 -------------------------------------------------------------------
rm(list = ls())
library(data.table)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
load('TCGA-symbol.rda')
icp <- fread(file='Immune checkpoints.txt',data.table = F)
co_sti <- icp[icp$Type=='Co-stimulatory immune checkpoint targets',][[2]]
co_inh <- icp[icp$Type=='Co-inhibitory immune checkpoint targets',][[2]]
hla <- c('HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQA2','HLA-DQB1','HLA-DQB2',
         'HLA-DRA','HLA-DRB1','HLA-DRB3','HLA-DRB4','HLA-DRB5')

co_sti <- intersect(rownames(my_data),co_sti) #27
co_inh <- intersect(rownames(my_data),co_inh) #15
hla <- intersect(rownames(my_data),hla)  # 9

dd1 <- my_data[co_sti,]
dd2 <- my_data[co_inh,]
dd3 <- my_data[hla,]

# 绘制共刺激免疫检查点箱图 ---
load('20-gene.rda')
dd1 <- merge(clin_edger[,c(1,2)],t(dd1),by.x=1,by.y=0)
dd1 <- tibble::column_to_rownames(dd1,'Sample')
C1 <- dd1[dd1$Cluster=='C1',]
C1 <- as.data.frame(t(C1[,-1]))
C1$sum <- apply(C1,1,mean)
C1$Group <- rep('C1',27)
C1 <- C1%>%
  dplyr::select(Group,sum)%>%
  tibble::rownames_to_column('Gene')

C2 <- dd1[dd1$Cluster=='C2',]
C2 <- as.data.frame(t(C2[,-1]))
C2$sum <- apply(C2,1,mean)
C2$Group <- rep('C2',27)
C2 <- C2%>%
  dplyr::select(Group,sum)%>%
  tibble::rownames_to_column('Gene')

sti_data <- rbind(C1,C2)

sti_data$sum <- round(sti_data$sum,2)
ggplot(sti_data,aes(x = Gene,y = sum,fill = Group))+
  geom_col(position = 'dodge')+
  geom_text(aes(label=sum),
            colour='black',size=2,
            vjust=1.5,position = position_dodge(0.9))+
  theme_bw()+
  scale_fill_manual(values = c("#FDBF6F", "#1F78B4"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid = element_blank())+
  ylab('Co-stimulatory immune checkpoint expression')+
  xlab('')
graph2pdf(file = 'Figure/immune/co-sti.pdf',width = 10, height = 7)

# 绘制hla家族点箱图 ---

dd3 <- merge(clin_edger[,c(1,2)],t(dd3),by.x=1,by.y=0)
dd3 <- tibble::column_to_rownames(dd3,'Sample')
C1 <- dd3[dd3$Cluster=='C1',]
C1 <- as.data.frame(t(C1[,-1]))
C1$sum <- apply(C1,1,mean)
C1$Group <- rep('C1',9)
C1 <- C1%>%
  dplyr::select(Group,sum)%>%
  tibble::rownames_to_column('Gene')

C2 <- dd3[dd3$Cluster=='C2',]
C2 <- as.data.frame(t(C2[,-1]))
C2$sum <- apply(C2,1,mean)
C2$Group <- rep('C2',9)
C2 <- C2%>%
  dplyr::select(Group,sum)%>%
  tibble::rownames_to_column('Gene')

hla_data <- rbind(C1,C2)

hla_data$sum <- round(hla_data$sum,2)
hla_data <- hla_data[!hla_data$Gene%in%c('CD58','TNFSF9'),]

ggplot(hla_data,aes(x = Gene,y = sum,fill = Group))+
  geom_col(position = 'dodge')+
  geom_text(aes(label=sum),
            colour='black',size=2,
            vjust=1.5,position = position_dodge(0.9))+
  theme_bw()+
  scale_fill_manual(values = c("#FDBF6F", "#1F78B4"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid = element_blank())+
  ylab('Human leukocyte antigen expression')+
  xlab('')
graph2pdf(file = 'Figure/immune/HLA.pdf',width = 7, height = 5)













