# TCGA --------------------------------------------------------------------
rm(list = ls())
library(ggplot2)
library(dplyr)
load('20-gene.rda')
load('clean-clin.rda')
clin <- merge(clin_edger[,1:2],clin,by = 1)
clin <- clin[,c(1:3,11:19)]
clin <- tibble::column_to_rownames(clin,'Sample')
# 按Cluster分成C1和C2，计算各列数值。
gname <- "Cluster"
dat <- clin
vname <- setdiff(colnames(dat), gname)
pie.C1 <- pie.C2 <- list()
fisher.p <- c()
for (i in vname) {
  tmp <- table(dat[,gname], dat[,i])
  p <- format(fisher.test(tmp)$p.value,digits = 2)
  names(p) <- i
  fisher.p <- c(fisher.p, p)

  pie.dat <-
    tmp %>% as.data.frame() %>% group_by(Var1) %>% mutate(Pct = Freq/sum(Freq)) %>% as.data.frame()

  # 表格内的两行对应Risk的两类：Risk high和Risk low
  pie.C1[[i]] <- pie.dat[which(pie.dat$Var1 == "C1"),]
  pie.C2[[i]] <- pie.dat[which(pie.dat$Var1 == "C2"),]
}
#设置颜色
library(scales)
library(RColorBrewer)
library(ggsci)
show_col(pal_d3()(10))
display.brewer.all()

black  <- "#1E1E1B"
blue   <- "#3C4E98"
yellow <- "#E4DB36"
orange <- "#E19143"
green  <- "#57A12B"
cherry <- "#8D3A86"

# 创建颜色
status.col <- c("grey80",black)
Age.col <- alpha(orange, c(0.5,  1))
Gender.col <- alpha('#17BECFFF', c(0.5,  1))
stage.col <- alpha(blue, c(0.5,  1))
Grade.col <- alpha('#008280FF', c(0.5,  1))
M.col <- alpha('#8C564BFF', c(0.5,  1))
N.col <- alpha(green, c(0.5, 1))
T.col <- alpha(cherry, c(0.5, 1))
# 硬核base plot一块一块画，当然也可以把其中的pie chart提取出来后期AI或者PPT拼接也是比较方便的
# pdf("pieTable.pdf",width = 7, height = 5)
showLayout <- T # 默认不在最终pdf的首页显示layout结构，不过建议初次绘制的时候改为TRUE看一下，方便理解
if(T){
# 设置画面布局，相同数字代表同一区块，数字越多代表该区块所占面积越大（一共25个区域）
layout(matrix(c( 1, 1, 1,  2, 2, 2,  3, 3, 3,  4, 4, 4,  5, 5, 5,  6, 6, 6,  7, 7, 7,  8, 8, 8,  9, 9, 9,
                 10,10,10, 11,11,11, 12,12,12, 13,13,13, 14,14,14, 15,15,15, 16,16,16, 17,17,17, 18,18,18,
                 10,10,10, 11,11,11, 12,12,12, 13,13,13, 14,14,14, 15,15,15, 16,16,16, 17,17,17, 18,18,18,
                 19,19,19, 20,20,20, 21,21,21, 22,22,22, 23,23,23, 24,24,24, 25,25,25, 26,26,26, 27,27,27,
                 19,19,19, 20,20,20, 21,21,21, 22,22,22, 23,23,23, 24,24,24, 25,25,25, 26,26,26, 27,27,27,
                 28,28,28, 29,29,29, 30,30,30, 31,31,31, 32,32,32, 33,33,33, 34,34,34, 35,35,35, 36,36,36,
                 37,37,37, 37,37,37, 37,37,37, 37,37,37, 37,37,37, 37,37,37, 37,37,37, 37,37,37, 37,37,37),
              byrow = T,nrow = 7))

if(showLayout) {
  layout.show(n = 37) # 直观展示画布分布
}

# 画布区域1-9：绘制图抬头 ##############
#-------------------------#

par(bty="n", mgp = c(0,0,0), mar = c(0,0,0,0), lwd = 2) # 基础参数，各边界距离为0
plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "HNSC",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Status",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Age",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Gender",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Stage",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Grade",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "T",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "N",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "M",cex = 2, col = "white") # 显示图标题

# 画布区域10-18：绘制C1组抬头和扇形图 ##############
#--------------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "C1\n(n = 291)",cex = 2, col = "white") # 显示图标题

# High group
pie(pie.C1$OS$Pct,
    col = status.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C1$Age$Pct,
    col = Age.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C1$Gender$Pct,
    col = Gender.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C1$Stage$Pct,
    col = stage.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C1$Grade$Pct,
    col = Grade.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C1$T$Pct,
    col = T.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C1$N$Pct,
    col = N.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C1$M$Pct,
    col = M.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

abline(v = par("usr")[2], col = "black") # 右侧封上黑线

# 画布区域19-27：绘制C2组抬头和扇形图 ##############
plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "C2\n(n = 203)",cex = 2, col = "white") # 显示图标题

# High group
pie(pie.C2$OS$Pct,
    col = status.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C2$Age$Pct,
    col = Age.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C2$Gender$Pct,
    col = Gender.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C2$Stage$Pct,
    col = stage.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C2$Grade$Pct,
    col = Grade.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C2$T$Pct,
    col = T.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C2$N$Pct,
    col = N.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.C2$M$Pct,
    col = M.col,
    border = "white",
    radius = 1,
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
abline(v = par("usr")[2], col = "black") # 右侧封上黑线


# 画布区域28-36：绘制空抬头和p值------------------------
plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p[1]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black") # 底部封上黑线

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p[2]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n",# 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p[3]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p[4]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p[5]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black") # 底部封上黑线


plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p[6]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black") # 底部封上黑线

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p[7]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black") # 底部封上黑线

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p[8]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black") # 底部封上黑线
abline(v = par("usr")[2], col = "black") # 右侧封上黑线

# 画布区域37：绘制图例 #-------------------------------------------

plot(0,0,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
legend("topleft",
       legend = c('','',
                  "Alive","Dead",
                  "≤65",">65",
                  "Female","Male",
                  "I-II","III-IV",
                  "G0-1","G2-3",
                  "T1-2","T3-4",
                  "N0-1","N2-3",
                  "M0","M1"),
       fill = c(c('white','white'),
                status.col,
                Age.col,
                Gender.col,
                stage.col,
                Grade.col,
                T.col,
                N.col,
                M.col),
       border = NA, # 图例颜色没有边框
       bty = "n", # 图例没有边框
       cex = 1.2,
       #box.lwd = 3,
       x.intersp = 0.05,
       y.intersp = 1,
       text.width = 0.075, # 图例的间隔
       horiz = T) # 图例水平放置

}

library(export)
graph2pdf(file='Figure/clin-cir/TCGA.pdf',width = 10,height =4)
dev.off()















































