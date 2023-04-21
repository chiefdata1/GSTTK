rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(maftools)
library(ggsci)
library(readxl)
library(MethylMix)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(clusterProfiler)
library(org.Hs.eg.db)
# -------------------------------------------------------------------------
load("methy/TCGA-HNSC_methy450.Rdata")
ee <- methy450
rm(methy450)
ee <- as.data.frame(ee)
data("Other")
xx <- as.data.frame(Other)
xx <- xx[,4:5]
xx$gene <- sapply(xx$UCSC_RefGene_Name, function(x)unlist(strsplit(x,';'))[1])
methyann <- data.frame(probe=rownames(xx),gene=xx$gene)
methyann <- na.omit(methyann)
rm(xx,Other)
ee <- merge(methyann,ee,by=1)
ee <- na.omit(ee)
ee <- ee%>%
  dplyr::select(-probe)%>%
  group_by(gene)%>%
  summarise_all(mean)%>%
  column_to_rownames('gene')
rm(methyann)
colnames(ee) <- substr(colnames(ee),1,16)
#读入原来整理的TCGA的tpm数据,scale之前的
load('TCGA-logtpm.rda')#scale之前的logtpm数据
load('HNSCC-tcga.rda')
rm(count,tpm)
logtpm <- as.data.frame(t(logtpm))
logtpm <- logtpm[rownames(logtpm)%in%clin$ID,]
rm(clin)
library(clusterProfiler)
library(org.Hs.eg.db)
x <- bitr(colnames(logtpm),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
logtpm <- as.data.frame(t(logtpm))
tpm <- merge(x,logtpm,by.x=1,by.y=0)
tpm <- distinct(tpm,SYMBOL,.keep_all = T)[,-1]
tpm <- column_to_rownames(tpm,"SYMBOL")
rm(logtpm)
range(tpm)

METnormal <- ee[,substr(colnames(ee),14,15)=="11"]
METcancer <- ee[,substr(colnames(ee),14,15)=="01"]
x <- intersect(colnames(tpm),colnames(METcancer))
METcancer <- METcancer[,x]
GEcancer <- log2(tpm[,x]+1)
y <- intersect(rownames(GEcancer),rownames(METcancer))
GEcancer <- GEcancer[y,]
METcancer <- METcancer[y,]
METnormal <- METnormal[y,]

save(GEcancer,METcancer,METnormal,file = 'methy/HNSC-MEdata.rda')

# -------------------------------------------------------------------------

MethylMixResults=MethylMix(as.matrix(METcancer),
                           as.matrix(GEcancer),
                           as.matrix(METnormal))
save(MethylMixResults,file="methy/MethylMixResults.Rda")
# -------------------------------------------------------------------------

sort(MethylMixResults$MethylationDrivers)#600个基因
tt <- data.frame(ss=substr(colnames(GEcancer),1,12),ID=colnames(GEcancer))
#提取分组信息
load('20-gene.rda')
Sinfo <- clin_edger[,c(1,2)]
rm(data_mut,clin_edger)
tt <- merge(Sinfo,tt,by.x=1,by.y=2)
colnames(tt) <- c('ID','Group','ss')
tt <- distinct(tt,ID,.keep_all = T)
G1 <- GEcancer[,colnames(GEcancer)%in%tt$ID[tt$Group=='C1']]
M1 <- METcancer[,colnames(METcancer)%in%tt$ID[tt$Group=='C1']]
G2 <- GEcancer[,colnames(GEcancer)%in%tt$ID[tt$Group=='C2']]
M2 <- METcancer[,colnames(METcancer)%in%tt$ID[tt$Group=='C2']]


# -------------------------------------------------------------------------

MethylMixResults1=MethylMix(as.matrix(M1),
                            as.matrix(G1),
                            as.matrix(METnormal))
MethylMixResults2=MethylMix(as.matrix(M2),
                            as.matrix(G2),
                            as.matrix(METnormal))


sort(MethylMixResults1$MethylationDrivers)
sort(MethylMixResults2$MethylationDrivers)

save(MethylMixResults,MethylMixResults1,MethylMixResults2,file = 'methy/MethylMix.Rda')

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(maftools)
library(ggsci)
library(readxl)
library(MethylMix)
library(impute)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# -------------------------------------------------------------------------

ann <- as.data.frame(Other)
ann <- na.omit(data.frame(Probe=rownames(ann),Gene=sapply(ann$UCSC_RefGene_Name,
                                                          function(x)unlist(strsplit(x,';'))[1])))

# -------------------------------------------------------------------------

load("methy/TCGA-HNSC_methy450.Rdata")
mm <- methy450
rm(methy450)
mm <- as.data.frame(mm)
mm <- na.omit(mm)

mm <- mm[rowMeans(mm[,substr(colnames(mm),14,15)=='11'])<=0.2,]
rownames(mm) <- NULL
mm <- column_to_rownames(mm,"V1")
t_mm <- mm[,substr(colnames(mm),14,15)=='01']

t_mm2 <- apply(t_mm, 2, function(x)ifelse(x>0.3,'M','U'))%>%as.data.frame()
t_mm2 <- t_mm2[apply(t_mm2, 1, function(x){sum(x=='M')>=0.05*ncol(t_mm2)}),]
t_mm3 <- merge(ann,t_mm2,by.x=1,by.y=0)

#读入原来整理的TCGA的tpm数据,scale之前的
load('Easy-Data/TCGA_m_l.rda')
tpm <- rbind(lnc_TCGA,m_TCGA)
rm(lnc_TCGA,m_TCGA)
logtpm <- as.data.frame(t(tpm))
logtpm <- logtpm[rownames(logtpm)%in%TCGA_clin$ID,]
rm(tpm)
library(clusterProfiler)
library(org.Hs.eg.db)
x <- bitr(colnames(logtpm),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
logtpm <- as.data.frame(t(logtpm))
tpm <- merge(x,logtpm,by.x=1,by.y=0)
tpm <- distinct(tpm,SYMBOL,.keep_all = T)[,-1]
tpm <- column_to_rownames(tpm,"SYMBOL")
rm(logtpm,TCGA_clin)
range(tpm)
HNSC_tpm <- log2(tpm[,substr(colnames(tpm),14,15)=='01']+1)

HNSC_tpm <- HNSC_tpm[apply(HNSC_tpm, 1, function(x){sum(x==0)<=0.1*ncol(HNSC_tpm)}),]
t_mm4 <- t_mm3[t_mm3$Gene%in%rownames(HNSC_tpm),]
colnames(t_mm4) <- substr(colnames(t_mm4),1,16)
t_mm4 <- t_mm4[,c('Probe','Gene',intersect(colnames(t_mm4),colnames(HNSC_tpm)))]
HNSC_tpm <- HNSC_tpm[,intersect(colnames(t_mm4),colnames(HNSC_tpm))]
rm(t_mm,t_mm2,t_mm3)
# -------------------------------------------------------------------------
################
################
################
probe <- data.frame()

system.time(
  for (i in t_mm4$Probe) {
    M_id <- colnames(t_mm4[,-c(1:2)])[which(t_mm4[t_mm4$Probe==i,-c(1:2)]=='M')]
    U_id <- colnames(t_mm4[,-c(1:2)])[which(t_mm4[t_mm4$Probe==i,-c(1:2)]!='M')]
    gene <- as.character(t_mm4$Gene)[t_mm4$Probe==i]
    M_expr_mean <- mean(as.numeric(HNSC_tpm[gene,M_id]))
    U_expr_mean <- mean(as.numeric(HNSC_tpm[gene,U_id]))
    M_expr_sd <- sd(as.numeric(HNSC_tpm[gene,M_id]))
    if((U_expr_mean - M_expr_mean) > 1.64*M_expr_sd){
      probe <- rbind(probe,cbind(P=i,Gene=gene,sig='y'))
    }
    else{
      probe <- rbind(probe,cbind(P=i,Gene=gene,sig='n'))
    }
  }
)

save(probe,file = "methy/probe.Rda")
load('methy/MethylMix.Rda')

gg <- c()
for (i in unique(probe$Gene)) {
  dd <- probe[probe$Gene==i,]
  if(sum(dd$sig=='y')>=length(dd$sig)/2){
    gg <- c(gg,i)
  }
}

drivers=sort(intersect(gg,MethylMixResults$MethylationDrivers))
save(drivers,t_mm,t_mm4,file = 'methy/Driver-methy.rda')

load('methy/Driver-methy.rda')

x <- MethylMixResults$MethylationDrivers
save(x,gg,file = 'methy/mm-divers.rda')
load('methy/mm-divers.rda')
y <- intersect(x,gg)
xx <- cbind(gg,x,y)
colnames(xx) <- c('The Wheeler criterion','The MethylMix algorithm','Intersection')
xx <- as.data.frame(xx)
xx$`The Wheeler criterion` <- c(sort(unique(xx$`The Wheeler criterion`)),rep(NA,600-153))
xx$`The MethylMix algorithm` <- sort(unique(xx$`The MethylMix algorithm`))
xx$Intersection <- c(sort(unique(xx$Intersection)),rep(NA,600-80))
write.csv(xx,'methy/Divers.csv',row.names = F)




# load('Driver-methy.rda')
# my <- read.csv('Diver-methy.csv')
# my2 <- merge(xx,my,by.x = 3,by.y = 1,all = T)
# my2 <- my2[,c(2,3,1,4:12,14,15,16)]
# write.csv(my2,'dd.csv',row.names = F)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
rm(list = ls())
library(dplyr)
library(tibble)
library(ggplot2)
library(ggsci)
library(ggthemes)
library(ggpubr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
load('methy/Driver-methy.rda')
tt <- t_mm4[t_mm4$Gene%in%drivers,]
colnames(t_mm) <- substr(colnames(t_mm),1,16)
t_mm <- t_mm[,colnames(t_mm4)[-c(1,2)]]
mm <- merge(tt[,1:2],t_mm,by.x=1,by.y=0)
rm(t_mm,t_mm4)
#读入原来整理的TCGA的tpm数据,scale之前的
load('Easy-Data/TCGA_m_l.rda')
tpm <- rbind(lnc_TCGA,m_TCGA)
rm(lnc_TCGA,m_TCGA)
logtpm <- as.data.frame(t(tpm))
logtpm <- logtpm[rownames(logtpm)%in%TCGA_clin$ID,]
rm(tpm)
library(clusterProfiler)
library(org.Hs.eg.db)
x <- bitr(colnames(logtpm),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
logtpm <- as.data.frame(t(logtpm))
tpm <- merge(x,logtpm,by.x=1,by.y=0)
tpm <- distinct(tpm,SYMBOL,.keep_all = T)[,-1]
tpm <- column_to_rownames(tpm,"SYMBOL")
rm(logtpm,TCGA_clin)
range(tpm)

ee <- log2(tpm[drivers,colnames(mm)[-c(1,2)]]+1)
rm(tpm)
mm2 <- mm[,-1]%>%
  group_by(Gene)%>%
  summarise_all(mean)%>%
  column_to_rownames('Gene')
colnames(ee) <- substr(colnames(ee),1,12)
colnames(mm2) <- substr(colnames(mm2),1,12)

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
#提取RS，计算分组信息
#提取分组信息
load('OS data/risk_score_list.rda')
Sinfo <- rsl$TCGA[,c(1,5)]
rm(rsl)
rownames(Sinfo) <- Sinfo$ID

Sinfo$ID <- substr(Sinfo$ID,1,12)
Sinfo <- distinct(Sinfo,ID,.keep_all = T)
rownames(Sinfo) <- substr(rownames(Sinfo),1,12)
x <- Reduce(intersect,list(rownames(Sinfo),colnames(mm2),colnames(ee)))


ee <- ee[,x]
mm2 <- mm2[,x]
Sinfo <- Sinfo[x,]
sub <- Sinfo
# -------------------------------------------------------------------------

cc <- data.frame()
for (i in drivers) {
  d <- cor.test(as.numeric(ee[i,]),as.numeric(mm2[i,]),method = 'spearman')
  d1 <- cor.test(as.numeric(ee[i,Sinfo$ID[Sinfo$Group=='High-risk']]),as.numeric(mm2[i,Sinfo$ID[Sinfo$Group=='High-risk']]),method = 'spearman')
  d2 <- cor.test(as.numeric(ee[i,Sinfo$ID[Sinfo$Group=='Low-risk']]),as.numeric(mm2[i,Sinfo$ID[Sinfo$Group=='Low-risk']]),method = 'spearman')
  cc <- rbind(cc,cbind(data.frame(ID=i,Corr=d$estimate,Pval=d$p.value,
                                  Corr1=d1$estimate,Pval1=d1$p.value,
                                  Corr2=d2$estimate,Pval2=d2$p.value)))
}

# -------------------------------------------------------------------------

x <- merge(Sinfo,t(mm2),by=0)[,-c(1)]

sig <- data.frame()
for (i in drivers) {
  sig <- rbind(sig,data.frame(ID=i,
                              mean1=mean(x[x$Group=='High-risk',i]),
                              mean2=mean(x[x$Group=='Low-risk',i]),
                              Pval=wilcox.test(x[,i]~x[,2])$p.value))
}
sig$sig <- ifelse(sig$Pval<0.001,'***',ifelse(sig$Pval<0.01,'**',ifelse(sig$Pval<0.05,'*','')))
colnames(sig) <- c('ID','Beta-mean1','Beta-mean2','Beta.KW.test.Pval','Beta.Sig')
my <- merge(cc,sig,by=1)

ee2 <- t(scale(t(ee)))
x <- merge(Sinfo,t(ee2),by=0)[,-c(1)]

sig <- data.frame()
for (i in drivers) {
  sig <- rbind(sig,data.frame(ID=i,
                              mean1=mean(x[x$Group=='High-risk',i]),
                              mean2=mean(x[x$Group=='Low-risk',i]),
                              Pval=wilcox.test(x[,i]~x[,2])$p.value))
}
sig$sig <- ifelse(sig$Pval<0.001,'***',ifelse(sig$Pval<0.01,'**',ifelse(sig$Pval<0.05,'*','')))
colnames(sig) <- c('ID','Zscore-mean1','Zscore-mean2','Zscore.KW.test.Pval','Zscore.Sig')
my <- merge(my,sig,by=1)

# -------------------------------------------------------------------------

x <- merge(Sinfo,t(ee2),by.x=0,by.y=0)[,-c(1)]
h1 <- gather(x,'Gene','Expr',-1)
x <- merge(Sinfo,t(mm2),by.x=0,by.y=0)[,-c(1)]
h2 <- gather(x,'Gene','Methy',-1)
h <- h1
h$Methy <- h2$Methy
rm(h1,h2)
ID <- my$ID[my$Beta.KW.test.Pval<0.00001&my$Zscore.KW.test.Pval<0.00001]#调整P值基因的数量也会调整

# PHYHD1画箱线图-------------------------------------------------------------------------
col <- pal_nejm()(2)
hh <- mm2['PHYHD1',]#分别画出两个基因的图
Sinfo <- Sinfo[order(Sinfo$Group),]
hh <- merge(Sinfo,t(hh),by.x=0,by.y=0)[,-c(1,2)]
ggdata <- gather(hh,'Gene','Value',-1)
ggdata$Group <- factor(ggdata$Group,levels = c('High-risk','Low-risk'))
options(digits = 2)
p1 <- ggplot(ggdata,aes(Gene,Value,fill=Group))+
  geom_boxplot(size=0.25,outlier.colour = 'grey30',outlier.size = 0.1,
               outlier.alpha = 0.6,lwd=1,fatten=1)+
  stat_compare_means(label = 'p.signif',method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")))+
  # scale_fill_npg(alpha = 0.8)+
  scale_fill_manual(values = col)+
  theme_classic()+
  xlab('')+
  ylab('Methylation Level (Beta)')+
  theme(axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=9),
        legend.position = 'top',
        legend.title = element_blank());p1

ggsave(filename = 'Figure/methy/PHYHD1-Diver-beta-boxplot.pdf',width = 3.5,height = 4.1)
dev.off()
# --
hh <- ee['PHYHD1',]
hh <- merge(Sinfo,t(hh),by.x=0,by.y=0)[,-c(1,2)]

ggdata1 <- gather(hh,'Gene','Value',-1)
ggdata1$Group <- factor(ggdata$Group,levels = c('High-risk','Low-risk'))


options(digits = 2)
p2 <- ggplot(ggdata1,aes(Gene,Value,fill=Group))+
  geom_boxplot(size=0.25,outlier.colour = 'grey30',outlier.size = 0.1,
               outlier.alpha = 0.6,lwd=1,fatten=1)+
  stat_compare_means(label = 'p.signif',method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")))+
  # scale_fill_npg(alpha = 0.8)+
  scale_fill_manual(values = col)+
  theme_classic()+
  xlab('')+
  ylab('mRNA normalized')+
  theme(axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=9),
        legend.position = 'top',
        legend.title = element_blank());p2

ggsave(filename = 'Figure/methy/PHYHD1-Diver-expr-boxplot.pdf',width = 3.5,height = 4.1)
library(patchwork)
p1/p2


# --
my2 <- my[my$ID=='PHYHD1',]
colnames(ggdata)[3] <- 'Beta'
colnames(ggdata1)[3] <- 'Expr'
ggdata$Expr <- ggdata1$Expr
ggdata1 <- ggdata
ggdata1$Group <- 'All'
ggdata <- rbind(ggdata,ggdata1)
ggdata$Group <- factor(ggdata$Group,levels = c('All','High-risk','Low-risk'))
options(digits = 1)
col1 <- pal_nejm()(8)[c(3,1,2)]
ggplot(ggdata,aes(Expr,Beta,color=Group))+
  geom_point(size=1.5,alpha=0.6)+
  stat_smooth(method = 'loess',color='grey40',size=0.5,se = F)+
  facet_grid(Gene~Group,scales = 'free')+
  theme_bw()+
  theme(legend.position = 'none',
        strip.text = element_text(size=11))+
  scale_color_manual(values = col1)+
  # scale_color_d3()+
  ylab('Methylation Level (Beta)')+
  xlab('mRNA Normalized')+
  labs(fill='Group')

ggsave(filename = 'Figure/methy/PHYHD1-E-M-Mehty-drivers.pdf',width = 4,height = 4)

write.csv(my,file = 'methy/Diver-methy.csv',row.names = F)

# --

msi <- t(mm2['PHYHD1',])
msi <- merge(sub,msi,by=0)
msi$msi <- ifelse(msi$PHYHD1>0.3,'M','U')
table(msi$Group,msi$msi)
options(digits = 3)
fisher.test(data.frame(c1=c(107,33),c2=c(217,133)))
ggdata <- data.frame(MSI=c('M','U'),msc1=c(107,33),msc2=c(217,133))
ggdata$msc1 <- ggdata$msc1/140
ggdata$msc2 <- ggdata$msc2/350
ggdata2 <- gather(ggdata,'Cluster','value',-1)
ggdata2$pp <- paste0(round(ggdata2$value,2)*100,'%')
ggdata2$MSI <- factor(ggdata2$MSI,levels = c('M','U'))
fix(ggdata2)
ggdata2$Cluster <- factor(ggdata2$Cluster,levels = c('msc1','msc2'))
ggplot(ggdata2,aes(Cluster,value,fill=MSI))+
  geom_bar(stat = 'identity',width = 0.85)+
  xlab('')+ylab('Fraction')+
  geom_text(aes(label=pp),size=3.8,
            position = position_stack(vjust = 0.5),
            color='white')+
  scale_fill_manual(values = col)+
  theme_bw()+
  scale_x_discrete(labels=c('High-risk','Low-risk'))+
  ggtitle('*')+
  theme(legend.position = 'right',
        plot.title = element_text(hjust=0.5,face = 'plain'),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12))+
  labs(fill='PHYHD1')
ggsave(filename = 'Figure/methy/Cluster-percent-PHYHD1-methy.pdf',height=7,width = 5)
# SVIP画箱线图-------------------------------------------------------------------------
#SVIP
col <- pal_nejm()(2)
hh <- mm2['SVIP',]#分别画出两个基因的图
Sinfo <- Sinfo[order(Sinfo$Group),]
hh <- merge(Sinfo,t(hh),by.x=0,by.y=0)[,-c(1,2)]
ggdata <- gather(hh,'Gene','Value',-1)
ggdata$Group <- factor(ggdata$Group,levels = c('High-risk','Low-risk'))
options(digits = 2)
p1 <- ggplot(ggdata,aes(Gene,Value,fill=Group))+
  geom_boxplot(size=0.25,outlier.colour = 'grey30',outlier.size = 0.1,
               outlier.alpha = 0.6,lwd=1,fatten=1)+
  stat_compare_means(label = 'p.signif',method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")))+
  # scale_fill_npg(alpha = 0.8)+
  scale_fill_manual(values = col)+
  theme_classic()+
  xlab('')+
  ylab('Methylation Level (Beta)')+
  theme(axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=9),
        legend.position = 'top',
        legend.title = element_blank());p1

ggsave(filename = 'Figure/methy/SVIP-Diver-beta-boxplot.pdf',width = 3.5,height = 4.1)
dev.off()
# --

hh <- ee['SVIP',]
hh <- merge(Sinfo,t(hh),by.x=0,by.y=0)[,-c(1,2)]

ggdata1 <- gather(hh,'Gene','Value',-1)
ggdata1$Group <- factor(ggdata$Group,levels = c('High-risk','Low-risk'))


options(digits = 2)
p2 <- ggplot(ggdata1,aes(Gene,Value,fill=Group))+
  geom_boxplot(size=0.25,outlier.colour = 'grey30',outlier.size = 0.1,
               outlier.alpha = 0.6,lwd=1,fatten=1)+
  stat_compare_means(label = 'p.signif',method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")))+
  # scale_fill_npg(alpha = 0.8)+
  scale_fill_manual(values = col)+
  theme_classic()+
  xlab('')+
  ylab('mRNA normalized')+
  theme(axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=9),
        legend.position = 'top',
        legend.title = element_blank());p2

ggsave(filename = 'Figure/methy/SVIP-Diver-expr-boxplot.pdf',width = 3.5,height = 4.1)
library(patchwork)
p1/p2


# --
my2 <- my[my$ID=='SVIP',]
colnames(ggdata)[3] <- 'Beta'
colnames(ggdata1)[3] <- 'Expr'
ggdata$Expr <- ggdata1$Expr
ggdata1 <- ggdata
ggdata1$Group <- 'All'
ggdata <- rbind(ggdata,ggdata1)
ggdata$Group <- factor(ggdata$Group,levels = c('All','High-risk','Low-risk'))
options(digits = 1)
col1 <- pal_nejm()(8)[c(3,1,2)]
ggplot(ggdata,aes(Expr,Beta,color=Group))+
  geom_point(size=1.5,alpha=0.6)+
  stat_smooth(method = 'loess',color='grey40',size=0.5,se = F)+
  facet_grid(Gene~Group,scales = 'free')+
  theme_bw()+
  theme(legend.position = 'none',
        strip.text = element_text(size=11))+
  scale_color_manual(values = col1)+
  # scale_color_d3()+
  ylab('Methylation Level (Beta)')+
  xlab('mRNA Normalized')+
  labs(fill='Group')

ggsave(filename = 'Figure/methy/SVIP-E-M-Mehty-drivers.pdf',width = 4,height = 4)

# --

msi <- t(mm2['SVIP',])
msi <- merge(sub,msi,by=0)
msi$msi <- ifelse(msi$SVIP>0.3,'M','U')
table(msi$Group,msi$msi)
options(digits = 3)
fisher.test(data.frame(c1=c(87,53),c2=c(138,212)))
ggdata <- data.frame(MSI=c('M','U'),msc1=c(87,53),msc2=c(138,212))
ggdata$msc1 <- ggdata$msc1/140
ggdata$msc2 <- ggdata$msc2/350
ggdata2 <- gather(ggdata,'Cluster','value',-1)
ggdata2$pp <- paste0(round(ggdata2$value,2)*100,'%')
ggdata2$MSI <- factor(ggdata2$MSI,levels = c('M','U'))
ggdata2$Cluster <- factor(ggdata2$Cluster,levels = c("msc1","msc2"))
fix(ggdata2)

ggplot(ggdata2,aes(Cluster,value,fill=MSI))+
  geom_bar(stat = 'identity',width = 0.85)+
  xlab('')+ylab('Fraction')+
  geom_text(aes(label=pp),size=3.8,
            position = position_stack(vjust = 0.5),
            color='white')+
  scale_fill_manual(values = col)+
  theme_bw()+
  scale_x_discrete(labels=c('High-risk','Low-risk'))+
  ggtitle('****')+
  theme(legend.position = 'right',
        plot.title = element_text(hjust=0.5,face = 'plain'),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12))+
  labs(fill='SVIP')
ggsave(filename = 'Figure/methy/Cluster-percent-SVIP-methy.pdf',height=7,width = 5)
