rm(list = ls())
library(magrittr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(Rmisc)
library(dplyr)
Segment <- read.delim("Mutation/TCGA-HNSC.cnv.tsv.gz") ##from UCSC
Segment$bases=Segment$End-Segment$Start
load('clean-clin.rda')
load('20-gene.rda')
rm(data_mut)
tcga <- clin_edger[,c(1,2)]
tcga$Sample <- substr(tcga$Sample,1,12)
Segment$sample <- substr(Segment$sample,1,12)

data=data.frame()
for (i in 1:length(table(Segment$sample))) {
  tmp=Segment[Segment$sample==names(table(Segment$sample))[i],]
  FGA=sum(tmp[abs(tmp$value)>0.2,"bases"])/ sum(tmp[,"bases"])
  FGG=sum(tmp[tmp$value>0.2,"bases"])/sum(tmp[,"bases"])
  FGL=sum(tmp[tmp$value< -0.2,"bases"])/sum(tmp[,"bases"])
  tmp=data.frame(Patient=names(table(Segment$sample))[i],FGA=FGA,FGG=FGG,FGL=FGL)
  data=rbind(data,tmp)
}

data2 <- merge(tcga,data,by.x=1,by=1)
colnames(data2) <- c("ID",'Cluster',"FGA","FGG","FGL")
clin$ID <- str_sub(clin$ID,1,12)
data3 <- merge(clin[,c(1,10:16)],data2,by = 1)
data3 <- data3[!duplicated(data3$ID),]


#自定义函数计算mean、se和p值
if(T){
cal_ms <- function(x,y){
  table1=data.frame()
  for (i in y) {
    name=summarySE(x,measurevar=c(colnames(x)[ncol(x)]), groupvars=c(colnames(x)[i]))[[1]]
    mean=summarySE(x, measurevar=c(colnames(x)[ncol(x)]), groupvars=c(colnames(x)[i]))[[3]] #mean
    se=summarySE(x, measurevar=c(colnames(x)[ncol(x)]), groupvars=c(colnames(x)[i]))[[5]] #se

    table=data.frame(group=rep(colnames(x)[i],length(mean)),subgroup=name,mean=mean,se=se)
    table1=rbind(table1,table)

  }
  table1$subgroup=factor(table1$subgroup,levels=unique(table1$subgroup))
  table1$group=factor(table1$group,levels=unique(table1$group))
  return(table1)
}#计算mean和se,排序
cal_p=function(x,y){
  table2=data.frame()
  for (i in y) {
    m<-aov(x[,ncol(x)]~as.factor(x[,i]))
    summary(m)
    table=cbind(rep(colnames(x)[i],nrow((TukeyHSD(m))[[1]])),(TukeyHSD(m))[[1]])
    table2=rbind(table2,table)
  }
  return(table2)
}#计算p值
signif=function(table1,table2){
  signif1=c()
  for (i in 1:length(table(table2[,1]))) {
    signif1=c(signif1,"Ref",table2[table2[,1]==levels(table1[,1])[i],5][1:length(table1[table1[,1]==names(table(table2[,1])[i]),"subgroup"])-1])
  }#提取p值
  for (i in 1:length(signif1)) {
    if (!is.na(as.numeric(signif1[i]))) {
      signif1[i]=ifelse(signif1[i]>0.05,"",ifelse(signif1[i]<=0.0001,"****",ifelse(signif1[i]<=0.001,"***",ifelse(signif1[i]<=0.01,"**",ifelse(signif1[i]<=0.05,'*')))))
    }
  }#转换标记
  return(signif1)
}#准备显著性标记的数据
mar=function(tableGL){
  mar=round(max(tableGL[,"mean"]+tableGL[,"se"]),1)+0.05
  return(mar)
}#ggplot坐标长度
cal_ms <- function(x,y){
  table1=data.frame()
  for (i in y) {
    name=summarySE(x,measurevar=c(colnames(x)[ncol(x)]), groupvars=c(colnames(x)[i]))[[1]]
    mean=summarySE(x, measurevar=c(colnames(x)[ncol(x)]), groupvars=c(colnames(x)[i]))[[3]] #mean
    se=summarySE(x, measurevar=c(colnames(x)[ncol(x)]), groupvars=c(colnames(x)[i]))[[5]] #se

    table=data.frame(group=rep(colnames(x)[i],length(mean)),subgroup=name,mean=mean,se=se)
    table1=rbind(table1,table)

  }
  table1$subgroup=factor(table1$subgroup,levels=unique(table1$subgroup))
  table1$group=factor(table1$group,levels=unique(table1$group))
  return(table1)
}#计算mean和se,排序
cal_p=function(x,y){
  table2=data.frame()
  for (i in y) {
    m<-aov(x[,ncol(x)]~as.factor(x[,i]))
    summary(m)
    table=cbind(rep(colnames(x)[i],nrow((TukeyHSD(m))[[1]])),(TukeyHSD(m))[[1]])
    table2=rbind(table2,table)
  }
  return(table2)
}#计算p值
signif=function(table1,table2){
  signif1=c()
  for (i in 1:length(table(table2[,1]))) {
    signif1=c(signif1,"Ref",table2[table2[,1]==levels(table1[,1])[i],5][1:length(table1[table1[,1]==names(table(table2[,1])[i]),"subgroup"])-1])
  }#提取p值
  for (i in 1:length(signif1)) {
    if (!is.na(as.numeric(signif1[i]))) {
      signif1[i]=ifelse(signif1[i]>0.05,"",ifelse(signif1[i]<=0.0001,"****",ifelse(signif1[i]<=0.001,"***",ifelse(signif1[i]<=0.01,"**",ifelse(signif1[i]<=0.05,'*')))))
    }
  }#转换标记
  return(signif1)
}#准备显著性标记的数据
mar=function(tableGL){
  mar=round(max(tableGL[,"mean"]+tableGL[,"se"]),1)+0.05
  return(mar)
}#ggplot坐标长度
}
data4 <- na.omit(data3)
data4[,2:9] <- apply(data4[,2:9],2,as.factor)
str(data4)
colnames(data4)[7] <- 'pN'

data4$Age <- ifelse(data4$Age=='0','≤65','>65')
data4$Gender <- ifelse(data4$Gender=='0','Female','Male')
data4$Stage <- ifelse(data4$Stage=='0','I-II','III-IV')
data4$Grade <- ifelse(data4$Grade=='0','G0-1','G1-2')
data4$T <- ifelse(data4$T=='0','T1-2','T3-4')
data4$pN <- ifelse(data4$pN=='0','N0-1','N2-3')
data4$M <- ifelse(data4$M=='0','M0','M1')

# 计算mean和se和p值并画图
table1 <- cal_ms(data4[,c(1:9,10)],2:9) #(保持最后一列是目标变量)
table2=cal_p(data4[,c(1:9,10)],2:9)#(保持最后一列是目标变量)
signif1 <- rep('',16)
table2$`p adj` <- as.numeric(table2$`p adj`)
p_adj <- ifelse(table2$`p adj`<0.0001,'****',
                ifelse(table2$`p adj`<0.001,'***',
                       ifelse(table2$`p adj`<0.01,'**',
                              ifelse(table2$`p adj`<0.05,'*','ns'))))
signif1[seq(1,16,2)] <- p_adj

table3=cal_ms(data4[,c(1:9,12)],2:9) #FGL
table4=cal_p(data4[,c(1:9,12)],2:9)  #FGL

signif2 <- rep('',16)
table4$`p adj` <- as.numeric(table4$`p adj`)
p_adj <- ifelse(table4$`p adj`<0.0001,'****',
                ifelse(table4$`p adj`<0.001,'***',
                       ifelse(table4$`p adj`<0.01,'**',
                              ifelse(table4$`p adj`<0.05,'*','ns'))))
signif2[seq(1,16,2)] <- p_adj

table5=cal_ms(data4[,c(1:9,11)],2:9) #FGG
table6=cal_p(data4[,c(1:9,11)],2:9)  #FGG
signif3 <- rep('',16)
table6$`p adj` <- as.numeric(table6$`p adj`)
p_adj <- ifelse(table6$`p adj`<0.0001,'****',
                ifelse(table6$`p adj`<0.001,'***',
                       ifelse(table6$`p adj`<0.01,'**',
                              ifelse(table6$`p adj`<0.05,'*','ns'))))
signif3[seq(1,16,2)] <- p_adj
tableGL=rbind(cbind(table3,class=rep('FGL',nrow(table3))),cbind(table5,class=rep('FGG',nrow(table5)))) #合并FGL、FGG




gap=0.05 #图纵坐标间隔

# 左侧部分(FGA)
if (T) {
  p1=ggplot(table1, aes(x = subgroup,
                        y = -mean,fill=rep("0",nrow(table1))))+
    geom_bar(stat = 'identity') +
    geom_errorbar(aes(ymax = -mean -se, ymin = -mean+se),position = position_dodge(0.9), width = 0.15)+
    scale_x_discrete(name="",position = "top")+
    geom_text(aes(y = -mean-se-0.007, label = signif1),position = position_dodge(0.9), size = ifelse(signif1=='ns',1.8,3), fontface = "bold",angle=90)+
    theme_bw()+
    theme(axis.line.y =element_line(size=0.8),
          axis.ticks.y =element_line(size=0.2),
          axis.text.y = element_blank(),
          axis.title.x = element_text(vjust = -2),
          plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
          legend.title=element_blank())+
    coord_flip()+
    scale_fill_manual(values  = c("#FFD034"),breaks=c("0"),labels=c("Copy number-altered genome"))+
    scale_y_continuous(limits = c(-mar(table1),0),
                       breaks =seq(-0.5,0,gap),labels=abs(seq(-0.5,0,gap)),expand = c(0.01,0),
                       name  = "FGA (Fraction of Genome Altered)",position = "left")

  p3=p1 +
    facet_grid(rows = vars(group),scales = "free_y",space = "free_y")+
    theme(strip.text.y = element_text(size = 10, angle= -90,vjust =1))
}#画图
p3

# 右侧部分(FGG/FGL)
if (T) {

  p2 <- ggplot(tableGL, aes(x = subgroup,
                            y = ifelse(class=='FGG',-mean,mean),fill=class))+
    geom_bar(stat = 'identity') +
    geom_errorbar(data=tableGL[tableGL$class=='FGL',],aes(ymax = mean+se, ymin =mean-se),position = position_dodge(0.9), width = 0.15)+
    geom_errorbar(data=tableGL[tableGL$class=='FGG',],aes(ymax = -mean-se, ymin =-mean+se),position = position_dodge(0.9), width = 0.15)+
    geom_text(data=tableGL[tableGL$class=='FGG',],aes(y = -mean -se-0.01, label = signif2),position = position_dodge(0.9), size = ifelse(signif2=='ns',1.8,3), fontface = "bold",angle=90)+
    geom_text(data=tableGL[tableGL$class=='FGL',],aes(y = mean+se+0.01, label = signif3),position = position_dodge(0.9), size = ifelse(signif3=='ns',1.8,3), fontface = "bold",angle=-90)+
    scale_x_discrete(name="")+
    theme_bw()+
    theme(axis.line.y =element_line(size=0.8),
          axis.ticks.y =element_line(size=0.2),
          axis.text.y = element_blank(),
          axis.title.x = element_text(vjust = -2),
          plot.margin = unit(c(0.3, 0.3, 0.3, -1), "lines"),
          legend.title=element_blank())+
    coord_flip()+
    scale_fill_manual(values  = c("#f17d80","#79BEDB"),breaks=c("FGG","FGL"),
                      labels=c("Copy number-Gained genome","Copy number-lost genome"))+
    scale_y_continuous(limits = c(-mar(tableGL[tableGL$class=='FGG',]),mar(tableGL[tableGL$class=='FGL',])),breaks =seq(-0.5,0.5,gap),
                       labels=abs(seq(-0.5,0.5,gap)),expand = c(0.01,0),
                       name  = "FGG or FGL (Fraction of Genome Gained or lost)")


  p4 <- p2+
    facet_grid(rows = vars(group),scales = "free_y",space = "free_y",switch = "y")+
    theme(strip.text.y = element_text(size = 10, angle= -90,vjust = 1))

}
p4

# 中间部分
pp=ggplot()+
  geom_text(data = tableGL,
            aes(label = subgroup, x=subgroup), y = 0.5,
            size = 0.8*11/.pt, # match font size to theme
            hjust = 0.5, vjust = 0.5)+
  theme_minimal()+
  theme(axis.line.y =element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.y =element_blank(),
        axis.title.y =element_blank(),
        axis.title.x =element_blank(),
        plot.margin = unit(c(0.3, 0, 0.3, 0), "lines")
  )+
  coord_flip()+
  scale_y_reverse()
ppa=pp+
  facet_grid(rows = vars(group),scales = "free_y",space = "free_y",switch = "y")+
  theme(strip.text.y = element_text(size = 0, angle= 45,vjust = 200))
library(patchwork)
#输出
pal <- p3 + ppa + p4 +
  plot_layout(widths = c(7,1,7),guides = 'collect')& theme(legend.position = 'bottom')
pal

#ggsave(file='FGA.png',pal,height = 6,width=12,units = c("in"),dpi = 1200)
ggsave(file='Figure/FGA_plus.pdf',pal,height = 6,width=12,units = c("in"))




















