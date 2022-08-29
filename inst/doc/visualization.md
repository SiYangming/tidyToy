---
title: "visualization"
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{visualization}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(tidyToy)
```

# 散点图

## plot函数绘制散点图示例

``` r
# 生成测试数据
x <- runif(100,0,2)
y <- runif(100,0,2)
## 第1种：使用plot一步绘制散点图
plot(x, y, main="散点图", sub="子标题", 
     xlab="横坐标", ylab="纵坐标", pch=16)
# 在图中标出指定的点
text(0.6,0.6,"(0.6,0.6)")
```

``` r
## 第1种：使用plot一步绘制散点图
plot(x, y, main="散点图", sub="子标题", 
     xlab="横坐标", ylab="纵坐标", pch=16)

# 在图中标出指定的点
text(0.6,0.6,"(0.6,0.6)")

# 在横纵轴指定处绘制红色线条
abline(h=.6,v=.6, col='red')
```

``` r
## 第2种：使用plot分步绘制散点图

# 生成空白绘图面板
plot(x, y, type="n", axes=F)
```

``` r
## 第2种：使用plot分步绘制散点图

# 生成空白绘图面板
plot(x, y, type="n", axes=F) 

#添加坐标点
points(x,y) 
```

``` r
## 第2种：使用plot分步绘制散点图

# 生成空白绘图面板
plot(x, y, type="n", axes=F) 

#添加坐标点
points(x,y) 
#添加横轴
axis(side=1) 
```

``` r
## 第2种：使用plot分步绘制散点图

# 生成空白绘图面板
plot(x, y, type="n", axes=F) 

#添加坐标点
points(x,y) 
#添加横轴
axis(side=1) 
#添加纵轴
axis(at=seq(0,2,0.5), side=2) 
```

``` r
## 第2种：使用plot分步绘制散点图

# 生成空白绘图面板
plot(x, y, type="n", axes=F) 

#添加坐标点
points(x,y) 
#添加横轴
axis(side=1) 
#添加纵轴
axis(at=seq(0,2,0.5), side=2) 
#补齐散点图的边框
box() 
```

``` r
## 第2种：使用plot分步绘制散点图

# 生成空白绘图面板
plot(x, y, type="n", axes=F) 

#添加坐标点
points(x,y) 
#添加横轴
axis(side=1) 
#添加纵轴
axis(at=seq(0,2,0.5), side=2) 
#补齐散点图的边框
box() 
# 添加图标题
title(main="散点图", sub="子标题", 
      xlab="横坐标", ylab="纵坐标")
```

``` r
## 第2种：使用plot分步绘制散点图

# 生成空白绘图面板
plot(x, y, type="n", axes=F) 

#添加坐标点
points(x,y) 
#添加横轴
axis(side=1) 
#添加纵轴
axis(at=seq(0,2,0.5), side=2) 
#补齐散点图的边框
box() 
# 添加图标题
title(main="散点图", sub="子标题", 
      xlab="横坐标", ylab="纵坐标")

# 在图中标出指定的点
text(0.6,0.6,"(0.6,0.6)")
```

``` r
## 第2种：使用plot分步绘制散点图

# 生成空白绘图面板
plot(x, y, type="n", axes=F) 

#添加坐标点
points(x,y) 
#添加横轴
axis(side=1) 
#添加纵轴
axis(at=seq(0,2,0.5), side=2) 
#补齐散点图的边框
box() 
# 添加图标题
title(main="散点图", sub="子标题", 
      xlab="横坐标", ylab="纵坐标")

# 在图中标出指定的点
text(0.6,0.6,"(0.6,0.6)")

# 在横纵轴指定处绘制红色线条
abline(h=.6,v=.6, col='red')
```

## nature文章实例

<https://www.nature.com/articles/nature10098> <https://static-content.springer.com/esm/art%3A10.1038%2Fnature10098/MediaObjects/41586_2011_BFnature10098_MOESM298_ESM.ppt>

b, In the replicate experiment mRNA levels explained 37% of protein levels in NIH3T3 cells (middle bar in a).

c, The model explains 85% of variance in protein levels from measured mRNA levels (middle bar in a).

Nature 473,337-342(19 May 2011)

``` r
# R基础绘图版本
a <- read.delim("data/nature10098-s5.txt")
#png("figs/mRNA-Protein_Scatter.png")
par(mgp=c(1.6,0.6,0),mar=c(3,3,1,1))
x<-a$mRNA.copy.number.replicate..molecules.cell.
y<-a$Protein.copy.number.replicate..molecules.cell.
plot(x,y,log="xy",pch=15,bty="l",
     xlab="mRNA copies per cell,replicate",
     ylab="Protein copies per cell,replicate",
     cex=0.7,font.lab=2,font.axis=2,
     xaxt="n",yaxt="n",#不绘制X，Y轴
     ylim=c(10^2,10^7))
axis(1,at=c(1,10^0.5,10,10^1.5,100,10^2.5,1000),
     labels=c("1","","10","","100","","1,000"))
axis(2,at=c(10^2,10^3,10^4,10^5,10^6,10^7),
     labels=c(expression(10*phantom()^2), expression(10*phantom()^3),
              expression(10*phantom()^4), expression(10*phantom()^5),
              expression(10*phantom()^6), expression(10*phantom()^7)))

#线性拟合
lm.obj<-lm(y~x)
summary(lm.obj)
sprintf("%.2f",0.2844)
text(x=1,y=10^7,pos=4,labels=expression(paste(R*phantom()^2,"=",0.31,sep="")))
```

``` r
# ggplot2版本
library(ggplot2)
a <- read.delim("data/nature10098-s5.txt")

daf <- data.frame(
  x<-a$mRNA.copy.number.replicate..molecules.cell.,
  y<-a$Protein.copy.number.replicate..molecules.cell.)

ggplot(data = daf) +
  geom_point(mapping = aes(x=x, y=y), alpha=0.2) +
  #geom_point(mapping = aes(x=x, y=y), alpha=0.2,colour="blue") +
  scale_y_log10() +
  scale_x_log10() +
  xlab("mRNA copies per cell,replicate") +  
  ylab("Protein copies per cell,replicate") + 
  annotate(geom = "text",x = 5,y = 10^7,
           label = "R^(2)==0.31", 
           parse = TRUE, vjust = 1)
```

## 火山图

``` r
### 火山图 ###
# 第一种：plot基础版
#读取数据
gene_list<-read.delim("data/volcano_plot_GSE20986.txt")
#查看数据情况
head(gene_list)
```

``` r
#par(mar(3,3,2,1),mgp=c(1.6,0.6,0),tck=.01)
plot(gene_list$logFC,-log10(gene_list$P.Value),
     pch=16,#设置点的类型
     cex=0.5,#设置点的大小
     cex.lab=1.2,#设置坐标轴名称的字体大小
     font.lab=2,#设置坐标轴名称的字体类型，粗体
     col=ifelse(gene_list$P.Value<=0.01#设置点的颜色
                &abs(gene_list$logFC)>1,"blue","gray"),
     xlim=c(-10,10),ylim=c(0,15),#Set limits
     xlab="log2 fold chage",
     ylab="-log10p-value")#设置坐标轴标签
```

``` r
# 第二种：
diff0 = read.table("data/volcano_plot_GSE20986.txt",sep="\t",header=T)
P.value = diff0$adj.P.Val
FC = diff0$logFC
df <- data.frame(P.value,FC)
df$threshold = as.factor(abs(df$FC) > log(1.5,2) & df$P.value < 0.01)
levels(df$threshold) = c('#f1f1f3','red')
df$logp = -log10(df$P.value)

plot(x = df[,2],y = df[,4], pch=16, col=df$threshold, 
     xlab="Fold change", ylab="-Log10(Pvalue)",
     cex=0.2, main="Volcano Plot")
abline(v=c(-log(1.5,2),log(1.5,2)), h=-log10(0.01),col="green")
```

``` r
# ggpolt2版本
library(ggplot2)
ggplot(data=df)+
  geom_point(aes(x=FC,y=-log10(P.value)),colour=df$threshold)
```

## 相关性图

``` r
## 相关性散点图
rt=read.table("data/corplot.txt",sep="\t",header=T,check.names=F,row.names=1)
i=1  #gene 1
j=2  #gene 2 
x=as.numeric(rt[i,])
y=as.numeric(rt[j,])
xName=row.names(rt[i,])
yName=row.names(rt[j,])
R=cor(x,y)
R=round(R,3)

#pdf(file="figs/cor.pdf")
z=lm(y~x)
plot(x,y, type="p",pch=16,cex=0.8, col="blue",main=paste("Pearson's correlation=",R,sep=""),
     xlab=paste(xName,"mRNA expression"),ylab=paste(yName,"mRNA expression"))
lines(x,fitted(z),col="red")
#dev.off()

## corplot包
rt <- read.table("data/corplot2.txt",sep="\t",header=T,row.names=1)       #读取文件

library(corrplot)
M = cor(rt)

#pdf(file="figs/corpot2.pdf",width=10,height=10)
corrplot(M,type = "upper")         #type: "full", "upper" or "lower"
#dev.off()

#pdf(file="figs/corpotMixed.pdf",width=13.5,height=13.5)
corrplot.mixed(M)
#dev.off()
```

## 散点图_火山图

``` r
# ==========================================================
#
#      散点图
#      •   描述两个变量之间的关联
#      •   火山图
#      •   nature文章实例
#
# ==========================================================

##### 散点图 #####

### nature文章实例
# https://www.nature.com/articles/nature10098
# https://static-content.springer.com/esm/art%3A10.1038%2Fnature10098/MediaObjects/41586_2011_BFnature10098_MOESM298_ESM.ppt
# b, In the replicate experiment mRNA levels explained 37% of protein levels in NIH3T3 cells (middle bar in a).
# c, The model explains 85% of variance in protein levels from measured mRNA levels (middle bar in a).
### Nature 473,337-342(19 May 2011)
rm(list = ls())
# 1.R基础绘图版本
# 不等长
load("data/ScatterData.rda")
a = ScatterData
#png("figs/mRNA-Protein_Scatter.png")
#par(mgp=c(1.6,0.6,0),mar=c(3,3,1,1))
x<-a$mRNA.copy.number.experiment..molecules.cell.
y<-a$Protein.copy.number.replicate..molecules.cell.
plot(x, y,
     log="xy",# 对xy轴值取log
     pch=1, # 点的形状
     bty="l", # 边框的形状
     xlab="mRNA copies per cell,replicate",
     ylab="Protein copies per cell,replicate",
     cex=0.7,font.lab=2,font.axis=2,#坐标轴的文字大小
     xaxt="n",yaxt="n",##不绘制X，Y轴
     ylim=c(10^2,10^7))

# 加上x轴标签
axis(1,at=c(1,10^0.5,10,10^1.5,100,10^2.5,1000),
     labels=c("1","","10","","100","","1,000"))

# 加上y轴标签
axis(2,at=c(10^2,10^3,10^4,10^5,10^6,10^7),
     labels=c(expression(10*phantom()^2), expression(10*phantom()^3),
              expression(10*phantom()^4), expression(10*phantom()^5),
              expression(10*phantom()^6), expression(10*phantom()^7)))

##线性拟合
lm.obj<-lm(y~x)
summary(lm.obj)
sprintf("%.2f",0.3371)
text(x=1,y=10^7,pos=4,
     labels=expression(paste(R*phantom()^2,"=",0.33,sep="")))
# dev.off()

# 2.ggplot2版本
library(ggplot2)

## 提取用于绘图的数据
daf <- data.frame(
  x = a$mRNA.copy.number.replicate..molecules.cell.,
  y = a$Protein.copy.number.replicate..molecules.cell.)

## 调用ggplot包的geom_point绘制散点图
scatter0 <- ggplot(data = daf)
scatter0
scatter0 +
  geom_point(mapping = aes(x=x, y=y), alpha=0.2)

# 修改点的颜色
scatter0 +
  geom_point(mapping = aes(x=x, y=y),
             alpha=0.2,colour="blue")

scatter1 <- scatter0 +
  geom_point(mapping = aes(x=x, y=y), alpha=0.2)
scatter1

# 对x轴y轴的值取log
scatter2 <- scatter1 +
  scale_y_log10() +
  scale_x_log10() #+
scatter2


# 添加x轴y轴的标题
scatter3 <- scatter2 +
  xlab("mRNA copies per cell,replicate") +
  ylab("Protein copies per cell,replicate")

scatter3

# 添加注释
scatter4 <- scatter3 + annotate(geom = "text",
                                x = 3, y = 10^7,
           label = "R^2==0.33",
           parse = TRUE, #数学表达式
           vjust = 1) # 注释的位置
scatter4 +
  theme_bw()

ggsave("figs/01-scatter_ggplot.tiff",
       plot = scatter4,
         width = 6, height = 6)

##### 火山图 #####
### 第1种：plot基础版

##读取数据
load(file = "data/limma_fit2.rda")

##查看数据情况
head(allDiff)

# 坐标轴上限
xMax <- max(abs(allDiff$logFC))
yMax <- max(-log10(allDiff$adj.P.Val))
#y1Max <- max(allDiff$adj.P.Val)

## 提取绘图所需数据
x <- allDiff$logFC
y <- -log10(allDiff$adj.P.Val)
#y1 <- allDiff$adj.P.Val

## 通过颜色来区分上下调基因
colour1 <- ifelse(allDiff$P.Value <= 0.01
                 &abs(allDiff$logFC) > 2,
                 "green","gray")
colour2 <- ifelse(allDiff$P.Value>0.01,"gray",
                 ifelse( allDiff$logFC > 2,"red",
                         ifelse( allDiff$logFC < -2, "green", "gray") )
)
## 绘图火山图

# tiff(file="vol.tiff",
#      width = 12,            #图片的宽度
#      height =12,            #图片的高度
#      units ="cm",
#      compression="lzw",
#      bg="white",
#      res=600)
plot(x,y,
  #y, x,
     pch=16,#设置点的类型
     cex=0.5,#设置点的大小
     cex.lab=1.2,#设置坐标轴名称的字体大小
     font.lab=2,#设置坐标轴名称的字体类型，粗体
     #col=colour2,#设置点的颜色
     xlim=c(-xMax,xMax),ylim=c(0,yMax),#设置横纵坐标
  #ylim=c(-xMax,xMax),xlim=c(0,yMax),
     xlab="log2 fold chage",
     ylab="-log10p-value",#设置坐标轴标签
     main="Volcano")

# 标记符合条件的点
diffSub=subset(allDiff, P.Value<0.01 & logFC>2)
points(-log10(diffSub$P.Value),
       diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=subset(allDiff, P.Value<0.01 & logFC<(-2))
points(-log10(diffSub$P.Value),
       diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,v = 0, lty=2,lwd=3)

# dev.off()
# 在图上标记基因
grep("CD36",DiffGene$symbol)
x[grep("CD36",DiffGene$symbol)]
y[grep("CD36",DiffGene$symbol)]
text(5.780170,11.456527,"CD36")

text(-4.212683,10.842402,"DUSP6")

# 第2种：ggpolt2版本
library(ggplot2)

## 提取用于绘图的数据
allDiff$log10PValue <- -log10(allDiff$FDR)

## 设置差异倍数阈值
allDiff$change <- ifelse(allDiff$FDR>0.01,'stable',
                          ifelse( allDiff$logFC > 1,'up',
                                  ifelse( allDiff$logFC < -1,
                                          'down', 'stable') )
)

## 绘制火山图
volcano0 <- ggplot(data = DiffGene)

volcano0 + geom_point(aes(x = logFC,y = log10PValue),
                      colour = colour1)# 指定每个点的颜色

volcano0 + geom_point(aes(x = logFC,y = log10PValue),
             colour = colour2)# 指定每个点的颜色

# 根据分组自动生成颜色
DiffGene$change <- factor(DiffGene$change)
volcano1 <- volcano0 +
  geom_point(aes(x = logFC,y = log10PValue,
                 colour = change))

volcano1
# 添加x轴和y轴标题
volcano2 <- volcano1 +
  xlab("FoldChange") +
  ylab("P Value")

# 添加点的标签
DiffGene$label <- ifelse(DiffGene$P.Value < 0.001
                         & abs(DiffGene$logFC) >= 3,
                         DiffGene$symbol,"")
install.packages("ggrepel")
library(ggrepel)
volcano2 +
  geom_text_repel(mapping = aes(x = logFC,
                                y = log10PValue,
                                label = label),
                  data = DiffGene,
    size = 3,# 字体大小
    segment.color = "black") #标签颜色

ggsave("figs/01-volcano_ggplot.tiff",
       width = 6, height = 6)

### ggplot2绘制散点图
a=ggplot(allDiff,aes(logFC,-log10(FDR)))
b=a+ geom_point(aes(color =change))+
  xlim(-xMax,xMax) + ylim(0,yMax) +
  labs(title="Volcano",x="log2FC", y="-log10(FDR)")+
  theme(plot.title=element_text(hjust=0.5))+
  scale_color_manual(values =c("green","black", "red"))
b+coord_flip()+geom_vline(xintercept = 0,lty=2,lwd=1)

# 第3种：ggpubr版本
library(ggpubr)
ggscatter(data = DiffGene,
          x = "logFC", y = "log10PValue",
          color = "change",size = 0.5,
          #shape = 18,
          #label = "symbol", repel = T,
          #label.select = DiffGene$symbol[1:30] ,
          #label.select = c('CD36','DUSP6'), #挑选一些基因在图中显示出来
          #palette = c("#00AFBB", "#999999", "#FC4E07"),
          ylab = "-log10p.value")
#palette里面是RGB颜色十六进制编号，
#参考https://www.114la.com/other/rgb.htm
ggsave("figs/01-scatter_valcano_ggpubr.tiff",
       width = 6, height = 6)
```

# hist直方图

``` r
#-------------------------------------
# name: hist.R
# date: 2014-07-08
# copyright: genedenovo
#-------------------------------------

# hist
data = rnorm(5000)
hist(data)

# prob/freq, breaks
layout(matrix(c(1,2,3,4),ncol=2,byrow=T))
hist(data)
hist(data,col="wheat")
hist(data,col="wheat",prob=T)
hist(data,col="wheat",freq=F,breaks=seq(min(data)-0.3,max(data)+0.3,0.2))
dev.off()


hist(data,col="wheat",breaks=seq(min(data)-0.3,max(data)+0.3,0.2),freq=F)
lines(density(data),col="red",lwd=2)

#----------------------------------------
# name: hist.practise.R
# date: 2014/7/17
# copyright: genedenovo
# desc: used for hist practice
#----------------------------------------

# read data
load(file = "data/limma_fit2.rda")
values = allDiff$logFC
len = length(values)

# 差异倍数设定在10到-10
for (i in 1 : len)
{
    if (values[i] > 10 || values[i] == "Inf")
    {
        values[i] = 10;
    }
    else if (values[i] < -10 || values[i] == "-Inf")
    {
        values[i] = -10;
    }
}

# draw histogram
hist(values,col="lightgreen",breaks=seq(-10,10,0.5),freq=F)
lines(density(values),lwd=2,col="red")
```

# 柱状图

``` r
### barplot版本 ###
chrom_distribution <- read.table("data/chrom_distribution.txt",
                                 sep="\t", row.names=1, header=T)              

height <- t(as.matrix(chrom_distribution))

#pdf(file="figs/barplot_chrom_distribution.pdf",width =7,height = 7)

# 绘制柱状图
barplot(height)
barplot(height, beside = TRUE)

# 添加颜色
barplot(height, 
        col = c("green","brown","blue") )
barplot(height,
        beside = TRUE,
        col = c("green","brown","blue") )

# 添加图例
legend("topright", colnames(chrom_distribution), 
       cex = 1.2, fill = c("green","brown","blue") )
#dev.off()
```

``` r
### ggplot版本 ###
library(ggplot2)                                                   #引用包
#读取输入文件
GOResult <- read.table("data/GOResult.txt", 
                       header = T,sep="\t")                
#定义输出文件
#pdf(file="figs/barplot-GOResult.pdf",width=8,height=4) 
barplot0 <- ggplot(data = GOResult) + 
  geom_bar(aes(x=Term, y=Count, fill=PValue), 
           stat='identity')
barplot0  

# 坐标轴旋转
barplot0 + coord_flip()

# 修改柱的颜色
barplot0 + coord_flip() + 
  scale_fill_gradient(low="red", high = "blue")

# 去除横纵坐标标题
barplot1 <- barplot0 + coord_flip() +
  scale_fill_gradient(low="red", high = "blue") +
  xlab("") + ylab("")

barplot1

# 修改x轴和y轴标题
barplot2 <- barplot1 + 
  theme(axis.text.x=element_text(color="black", size=12),
        axis.text.y=element_text(color="black", size=12)) 
barplot2

# 修改柱与x轴和y轴的距离
barplot2 + 
  scale_y_continuous(expand=c(0, 0)) 

barplot2 + 
  scale_y_continuous(expand=c(0, 0)) + 
  scale_x_discrete(expand=c(0,0))
```

## 柱状图barplot

``` r
# ==========================================================
#
#      柱状图
#      •   某一类别的计数展示
#      •   barplot
#      •   ggplot2
#
# ==========================================================

# 清空变量环境
rm(list = ls())
# 关闭绘图窗口
dev.off()
load("data/barplotData.rda")
##### barplot绘制示例 #####
#读入数据
data = barplotData1
#绘制条形图
barplot(data[,2])
barplot(data[,2],names.arg = data[,1])
barplot(data[,2],names.arg = data[,1],main="条形图",xlab="分组",ylab="统计量")
barplot(data[,2],names.arg = data[,1],main="条形图",xlab="分组",ylab="统计量",col="blue")
barplot(data[,2],names.arg = data[,1],main="条形图",xlab="分组",ylab="统计量",col=c("grey","red","blue","orange","green"))

###----多种分组的柱状图
#堆积柱状图
#转换数据
data2 = t(data[,c(2,3)])
#绘制柱状图
barplot(as.matrix(data2))
#调整图形
barplot(as.matrix(data2),
        names.arg = data[,1],main="条形图",xlab="分组",ylab="统计量",
        col=c("blue","red"),
        legend=c("Low Dose","High Dose"),
        ylim=c(0,50))

box() #边框

#非堆积柱状图
barplot(as.matrix(data2),
        names.arg = data[,1],main="条形图",xlab="分组",ylab="统计量",
        col=c("blue","red"),
        legend=c("Low Dose","High Dose"),
        ylim=c(0,30),
        beside=TRUE)
box() #边框

###------水平柱状图
par(las=2)
barplot(as.matrix(data2),
        names.arg = data[,1],main="条形图",xlab="统计量",
        col=c("blue","red"),
        legend=c("Low Dose","High Dose"),
        xlim=c(0,50),
        horiz=TRUE,
        space=0.5)

box() #边框

##### barplot1绘制柱状图 #####
barplot1 <- barplotData2
n=as.numeric(barplot1[,1])
names(n)=row.names(barplot1)

par(mar=c(3,9,3,3),xpd=T)
bar=barplot(n,
            horiz=TRUE,
            col="skyblue",
            names=FALSE,
            xlab = "")
text(x=n*0.9,y=bar,n)
text(x=-0.2,y=bar,
     label=names(n),
     xpd=T,pos=2)

##### barplot2版本 #####
height <- t(as.matrix(chrom_distribution))

#pdf(file="figs/barplot_chrom_distribution.pdf",width =7,height = 7)

# 绘制柱状图
barplot(height)
barplot(height, beside = TRUE)

# 添加颜色
barplot(height,
        col = c("green","brown","blue") )
barplot(height,
        beside = TRUE,
        col = c("green","brown","blue") )

# 添加图例
legend("topright", colnames(chrom_distribution),
       cex = 1.2, fill = c("green","brown","blue") )
#dev.off()

##### 免疫细胞比例柱状图 #####
#outpdf="figures/barplot.pdf"
data=t(CibersortData)
col=rainbow(nrow(data),s=0.7,v=0.7)

#pdf(outpdf,height=6,width=10)
par(las=1,mar=c(8,5,4,13))
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.02,colnames(data),adj=1,cex=1);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
#text(par('usr')[2],(ytick1+ytick2)/2,rownames(data),cex=0.6,adj=0)
legend(par('usr')[2]*0.97,par('usr')[4],
       legend=rownames(data),col=col,
       pch=15,bty="n",cex=0.8)
#dev.off()

##### STRING基因统计 #####
tb=table(c(as.vector(string_interactions[,1]),as.vector(string_interactions[,2])))
tb=sort(tb,decreasing =T)
# write.table(tb,file="results/count.xls",
#             sep="\t",quote=F,col.names=F)

#定义可视化的基因数目
n=as.matrix(tb)[1:30,]

#pdf(file="barplot.pdf",width=8,height=6)
par(mar=c(3,10,3,3),xpd=T)
bar=barplot(n,horiz=TRUE,col="skyblue",names=FALSE)
text(x=n-0.5,y=bar,n)
text(x=-0.2,y=bar,label=names(n),xpd=T,pos=2)
#dev.off()

##### barplot绘制GO分组可视化 #####
x1 <- nrow(bp); x2=nrow(cc); x3=nrow(mf); x=x1+x2+x3
m <- c(bp[,2],cc[,2],mf[,2])
l <- c(bp[,1],cc[,1],mf[,1])
y <- ceiling(max(m)/10)*15
#tiff(file="pics/GOBarplot.tiff",width = 30,height = 20,units ="cm",compression="lzw",bg="white",res=300)
par(mar=c(8,4,3,3),mgp=c(0.8,0.3,0),cex.axis=.7)
barplot(m,
        beside=T,
        col=c(rep(rgb(153/255,216/255,201/255),x1),
              rep(rgb(44/255,127/255,184/255),x2),
              rep(rgb(201/255,148/255,199/255),x3)),
        space=0,xaxs='i',yaxs='i',yaxt='n',las=2,
        names.arg=l,ylab="target genes")
box()
abline(h=0)
par(xpd=T)
lx <- max(nchar(l))
y1 <- 5/100*y
y2 <- 10/100*y
segments(0,-y1,0,-y2); segments(0,-y2,x,-y2); segments(x1,-y1,x1,-y2); segments(x1+x2,-y1,x1+x2,-y2); segments(x,-y1,x,-y2)
text(x1/2,-(y2-2.5/10*y2),labels='biological_process',pos=1,cex=1,col=rgb(153/255,216/255,201/255))
text(x1+x2/2,-(y2-2.5/10*y2),labels='cellular_component',pos=1,cex=1,col=rgb(44/255,127/255,184/255))
text(x1+x2+x3/2,-(y2-2.5/10*y2),labels='molecular_function',pos=1,cex=1,col=rgb(201/255,148/255,199/255))
axis(2)
#dev.off()

##### ggplot版本：绘制GO柱状图 #####
library(ggplot2)                                                   #引用包
#读取输入文件
#绘制基础版柱状图
barplot3 <- ggplot(data = GOResult) +
  geom_bar(aes(x=Term, y=Count,
               fill=-log10(PValue)),
           stat='identity') +
  # 坐标轴旋转
  coord_flip() +
  # 修改柱的颜色
  scale_fill_gradient(low="red",
                      high = "blue") +
  # 去除横纵坐标标题
  xlab("") + ylab("") +
  # 修改x轴和y轴的字体颜色
  theme(axis.text.x=element_text(color="black", size=12),
        axis.text.y=element_text(color="black", size=12)) +
  # 修改柱与x轴和y轴的距离
  scale_y_continuous(expand=c(0, 0)) +
  scale_x_discrete(expand=c(0,0))

barplot3

# ggsave("figures/barplotGOResult.tiff", width = 8, height = 4)

```

## GO_barplot

``` r
#David网站分析后，Go可视化富集分析结果，BP CC MF分别提取Counts数前5的Term
setwd("E:\\clinical shengxin\\Keratitis\\GSE58291\\NetworkAnalyst")
up = read.table(file = 'go_all.txt',sep = '\t',header = T,quote = '')
up_rt = up[up$PValue < 0.05,]
library(tidyr)
up_rt = separate(up_rt, Term, sep ="~",
                 into = c("ID", "Term"))

bp_df = up_rt[up_rt$Category == 'GOTERM_BP_DIRECT',]
bp_df = bp_df[order(bp_df$Count,decreasing = T),]
bp = bp_df[1:10,]

cc_df = up_rt[up_rt$Category == 'GOTERM_CC_DIRECT',]
cc_df = cc_df[order(cc_df$Count,decreasing = T),]
cc = cc_df[1:10,]

mf_df = up_rt[up_rt$Category == 'GOTERM_MF_DIRECT',]
mf_df = mf_df[order(mf_df$Count,decreasing = T),]
mf = mf_df[1:10,]

allGo = rbind(bp,cc,mf)
library(stringr)
table(allGo$Category)
allGo$Category = substr(allGo$Category,8,9)

#画出条形图
library(ggpubr)
colnames(allGo)
p = ggbarplot(data = allGo,x = "ID",y = "Count",
              fill = "Category",
              palette = c("cadetblue3","mediumslateblue","mediumorchid3"),
              sort.by.groups = T,xlab = '',ylab = "Target genes")
p <- ggpar(p,x.text.angle = 45)
ggsave(plot = p,'barplot.pdf',width = 10,height =5)
```

# PCA图

``` r
## 1.各组织基因表达量的PCA分析
a <- read.table("data/mcp.M113.035600-2.txt",
                head = TRUE, sep = "\t",
                check.names = FALSE)
#第一列为基因ID，去掉第一列，
#把剩余的列保存到x里
x <- a[,c(-1,-ncol(a))]
#转置
x <- t(x)
x <- log2(x+1)
#加载包，利用包里的pca函数进行PCA分析
library(pcaMethods)
#进行PCA分析
pres <- pca(x,scale = "pareto",center = TRUE)
#提取主成分的贡献率
xlab <- paste("PC1",sprintf("(%.2f%%)",100*pres@R2[1]))
ylab <- paste("PC2",sprintf("(%.2f%%)",100*pres@R2[2]))
#绘制PCA的得分图
#pdf("figs/pca-mutil_tissue.pdf",width = 4,height = 4)
#par(mar=c(3,3,1,1),mgp=c(1.6,0.6,0))
plot(pres@scores[,1],pres@scores[,2],
     xlim=c(-150,150),
     xlab=xlab,ylab=ylab,pch=16,col="green")
abline(v=0,h=0,lty=2,col="gray",lwd=2)
text(pres@scores[,1],pres@scores[,2],labels = names(a)[-1],pos = 4,cex=0.5)
dev.off()

## 2.不同样本基因表达量的PCA分析
data <- read.table("data/exp.txt",header=T,sep="\t",row.names=1)   
#矩阵转置
data <- t(as.matrix(data))   
data.class <- rownames(data)
#PCA分析
data.pca <- prcomp(data, scale. = TRUE)   
#输出特征向量
write.table(data.pca$rotation,file="data/pca-PC.xls",quote=F,sep="\t")   
#输出新表
write.table(predict(data.pca),file="data/pca-newTab.xls",quote=F,sep="\t")   
pca.sum <- summary(data.pca)
#输出PC比重
write.table(pca.sum$importance,file="data/pca-importance.xls",quote=F,sep="\t")   

#柱状图
pdf(file="figs/pcaBarplot.pdf",width=15)   
barplot(pca.sum$importance[2,]*100,xlab="PC",ylab="percent",col="skyblue")
dev.off()

#碎石图
pdf(file="figs/pcaPlot.pdf",width=15)   
plot(pca.sum$importance[2,]*100,type="o",col="red",xlab="PC",ylab="percent")
dev.off()

# #pca 2d plot
# library(ggplot2)
# group <- c(rep("con",5),rep("A",5),rep("B",3))  
# pcaPredict <- predict(data.pca)
# PCA <- data.frame(PCA1 = pcaPredict[,1], PCA2 = pcaPredict[,2],group=group)
# PCA.mean <- aggregate(PCA[,1:2],list(group=PCA$group),mean)
# #定义椭圆
# veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
# {
#   theta <- (0:npoints) * 2 * pi/npoints
#   Circle <- cbind(cos(theta), sin(theta))
#   t(center + scale * t(Circle %*% chol(cov)))
# }
# df_ell <- data.frame()
# for(g in levels(PCA$group)){
#   df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$group==g,],
#                                                    veganCovEllipse(cov.wt(cbind(PCA1,PCA2),
#                                                                           wt=rep(1/length(PCA1),length(PCA1)))$cov,
#                                                                    center=c(mean(PCA1),mean(PCA2))))),group=g))
# }
# 
# pdf(file="figs/PCA2d.pdf")
# ggplot(data = PCA, aes(x=PCA1, y=PCA2)) + geom_point(aes(color = group)) +
#   geom_path(data=df_ell, aes(x=PCA1, y=PCA2,colour=group), size=1, linetype=2)+
#   annotate("text",x=PCA.mean$PCA1,y=PCA.mean$PCA2,label=PCA.mean$group)+
#   theme_bw()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# dev.off()
```

# 箱线图

``` r
# Nature Genetics 43,811-814(2011)  SPSS绘制

# boxplot版本
a<-read.delim("data/ng.864-S2.txt")
#数据预处理
d1<-data.frame(sample="Sperm",
               meth_ratio=a$Sperm.Rep1..Cs)
d2<-data.frame(sample="GV",
               meth_ratio=a$GV.Rep1..Cs)
d3<-data.frame(sample="MII",
               meth_ratio=a$MII.Rep1..Cs)
d4<-data.frame(sample="X3aWT",
               meth_ratio=a$X3aWT.Rep1..Cs)
md<-rbind(d1,d2)
md<-rbind(md,d3)
md<-rbind(md,d4)
#绘图
pdf("figs/boxplot.pdf",width=3.1,height=3.4)
par(mgp=c(1.6,0.6,0),mar=c(3,3,1,1))
boxplot(md[md[,2]!=0,2]~md[md[,2]!=0,1],
        boxwex=0.3,
        col="red",cex=0.6,
        xlab="Sample",ylab="Percent(%)")
dev.off()

#ggplot2版本
a<-read.delim("data/ng.864-S2.txt")
#数据预处理
d1<-data.frame(sample="Sperm",
               meth_ratio=a$Sperm.Rep1..Cs)
d2<-data.frame(sample="GV",
               meth_ratio=a$GV.Rep1..Cs)
d3<-data.frame(sample="MII",
               meth_ratio=a$MII.Rep1..Cs)
d4<-data.frame(sample="X3aWT",
               meth_ratio=a$X3aWT.Rep1..Cs)
md<-rbind(d1,d2)
md<-rbind(md,d3)
md<-rbind(md,d4)
md<-md[md[,2]!=0,]
library(ggplot2)
ggplot(data=md,mapping=aes(x=sample,y=meth_ratio))+
  geom_boxplot(fill="red",width=0.5) +  xlab("Sample") + ylab("Percent(%)")

# 不同样本间的基因表达
rt <- read.table("data/expression_Type.txt",sep="\t",header=T)                               #¶ÁÈ¡ÊäÈëÎÄ¼þ
#pdf(file = "boxplot-expression_Type.pdf",width=7,height=7)                                 #±£´æÊä³öÎÄ¼þ
boxplot(expression ~ Type, rt, col = c("green", "red"),outline = FALSE)
#dev.off()

## 小提琴图
library(ggplot2)
library(gplots)
library(RColorBrewer)

#read in the data file
data = read.table('data/violin_plot.txt', sep="\t", header=T)
#take a glance at the data
head(data)
dim(data)
data$EXP = log(data$ACIN1,2)
#start to draw the violin plot
p <- ggplot(data,aes(CANCER,EXP))+
  geom_violin(adjust=1,trim=T)+
  ggtitle("Expression of ACIN1 in different cancers")
print(p)

#changeing params
p <- ggplot(data,aes(CANCER,EXP))+
  geom_violin(adjust=1,trim=T)+
  geom_boxplot(width=0.3,fill="black",alpha=1,outlier.colour=NA)+
  stat_summary(fun.y=mean,geom="point",fill='white',shape=21,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60,hjust=1,size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab('Expression level')+
  xlab(NULL)+
  ggtitle("Expression of ACIN1 in different cancers")
print(p)

#changing params again, mainly change the color and sort the color by mean value
#library(RColorBrewer)
color1 = brewer.pal(12,"Set3")
color2 = brewer.pal(9,"Set1")
color = c(color1,color2)
p <- ggplot(data,aes(CANCER,EXP))+
  geom_violin(adjust=1,trim=T)+
  geom_boxplot(width=0.3,fill=color,alpha=1,outlier.colour=NA)+
  stat_summary(fun.y=mean,geom="point",fill='white',shape=21,size=2)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60,hjust=1,size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab('Expression level')+
  xlab(NULL)+
  ggtitle("Expression of ACIN1 in different cancers")
print(p)
```

## 箱线图系列boxplot_小提琴图vioplot

``` r
# ==========================================================
#
#      箱线图
#      •   查看数据分布
#      •   小提琴图
#      •   花生图
#
# ==========================================================
load("data/boxplotData.rda")
##### boxplot版本 #####

### 示例1                             #读取输入文件
#pdf(file = "boxplot.pdf",width=7,height=7)                                 #保存输出文件
boxplot(expression ~ Type,
        data=boxplotData1, col = c("green", "red"),outline = FALSE)
#dev.off()

### 示例2
boxplot(Value~Group,data = boxplotData2)

#调整
boxplot(Value~Group,data = boxplotData2,
        varwidth=TRUE,
        col="#6DC9F2",
        outline=FALSE,
        names=c("Control","Drug1","Drug2","Drug3","Drug4","Drug5"),
        horizontal=FALSE )
#高级调整
points(jitter(as.numeric(boxplotData2$Group)),boxplotData2$Value, pch = 16, col = "red")

##using ggplot2
ggplot(boxplotData2, aes(Group,Value))+
  geom_boxplot()


### 示例3
# Nature Genetics 43,811-814(2011)  SPSS绘制
a <- boxplotData3
##数据预处理
d1<-data.frame(sample="Sperm",
               meth_ratio=a$Sperm.Rep1..Cs)
d2<-data.frame(sample="GV",
               meth_ratio=a$GV.Rep1..Cs)
d3<-data.frame(sample="MII",
               meth_ratio=a$MII.Rep1..Cs)
d4<-data.frame(sample="X3aWT",
               meth_ratio=a$X3aWT.Rep1..Cs)
md<-rbind(d1,d2)
md<-rbind(md,d3)
md<-rbind(md,d4)
##绘图
#pdf("figs/boxplot.pdf",width=3.1,height=3.4)

boxplot(md[md[,2]!=0,2]~md[md[,2]!=0,1],
        boxwex=0.3,
        col="red",cex=0.6,
        xlab="Sample",ylab="Percent(%)")
#dev.off()

##ggplot2版本

library(ggplot2)
ggplot(data=md) +
  geom_boxplot(mapping=aes(x=sample,y=meth_ratio),
               fill="red", width=0.5) +
  xlab("Sample") + ylab("Percent(%)")

#ggsave("figures/04-boxplot.tiff")

##### ggpubr箱线图系列boxplot #####
library(ggpubr, quietly = TRUE)
#读取输入文件
#绘制基因在所有样品的boxplot
p=ggboxplot(boxplotData4, x="organ", y="value", color = "organ",
            palette = rainbow(length(levels(boxplotData4$organ))),
            ylab="Gene expression",
            xlab="",
            #add = "jitter",                                            #绘制每个样品的散点
            legend = "right")
# pdf(file="figures/boxplot.pdf",width=12,height=6)                          #输出图片文件
p+rotate_x_text(60)
# dev.off()

#根据男女分别绘制
p=ggboxplot(boxplotData4, x="organ", y="value", color = "Gender",
            ylab="Gene expression",
            xlab="",
            #add = "jitter",                                              #绘制每个样品的散点
            palette = c("green","red") )
p=p+stat_compare_means(aes(group=Gender),
                       label = "p.signif",
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")))
p=p+rotate_x_text(60)
# pdf(file="figures/boxplot.group.pdf",width=8,height=5)                       #输出图片文件
p
# dev.off()

##### 小提琴图vioplot #####
#引用包
library(vioplot, quietly = TRUE)
#正常样品数目
normal=6
tumor=6
#读取输入文件
load("data/PairedPlot.rda")
#保存图片的文件名称
# pdf("figures/小提琴图vioplot.pdf",height=8,width=15)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(df))
y=c(1:ncol(df))
plot(x,y,
     #xlim=c(0,37),ylim=c(min(df),max(df)*1.1),
     xlim=c(0,63),ylim=c(min(df),max(df)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     #cex.lab=1.5,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，正常用蓝色表示，肿瘤用红色表示
for(i in 1:ncol(df)){
  normalData=df[1:normal,i]
  tumorData=df[(normal+1):(normal+tumor),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.8)
  text(seq(1,64,3),-0.03,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
}
#text(seq(1,37,3),-0.03,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
#dev.off()

##### ggplot绘制小提琴图 #####
library(ggpubr)
Type=Type[order(Type[,2]),]

score=score[row.names(Type),]
colnames(Type)=c("cluster","Subtype")
cluster=cbind(Type,score)
cluster=cluster[,-1]

cluster$Subtype=factor(cluster$Subtype,
                       levels=c("Immunity_L","Immunity_M","Immunity_H"))
my_comparisons=list(c("Immunity_L","Immunity_M"),
                    c("Immunity_M","Immunity_H"),
                    c("Immunity_H","Immunity_L"))

#pdf(file="vioplot.pdf",width=6,height=5)
ggviolin(cluster, x="Subtype", y="TumorPurity",
         fill = "Subtype",
         palette = c("green", "blue", "red"),
         add = "boxplot", add.params = list(fill="white")) +
  stat_compare_means(comparisons =
                       my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                       symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")


##### 带统计的箱线图 #####
geneName="SMAD4"
cnv=paste(geneName,"|cnv",sep="")
exp=paste(geneName,"|exp",sep="")

rt=boxplotData5
group=sapply(strsplit(colnames(rt),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
rt=rt[,group==0]

data=rbind(cnv=rt[cnv,],exp=log2(rt[exp,]+1))
data=t(data)
ksTest<-kruskal.test(exp ~ cnv, data=data)
ksPval=ksTest$p.value
f=factor(data[,"cnv"])
labels=levels(f)
labels=gsub("-2","Double deletion",labels)
labels=gsub("-1","Single deletion",labels)
labels=gsub("0","Normal",labels)
labels=gsub("1","Single gain",labels)
labels=gsub("2","Amplification",labels)
cols=rainbow(length(labels))

pval=0
if(ksPval<0.001){
  pval=signif(ksPval,4)
  pval=format(pval, scientific = TRUE)
}else{
  pval=round(ksPval,3)
}

# tiffFile=paste("figures/",geneName,".tiff",sep="")
# tiff(file=tiffFile,width = 20,height = 17 ,units ="cm",compression="lzw",bg="white",res=600)
boxplot(exp ~ cnv,
        data = data,
        col=cols,
        names=labels,
        main=paste("p=",pval,sep=""),
        ylab = paste(geneName," expression",sep=""),
        xlab = paste(geneName," copy number",sep=""),
        cex.main=1.5,
        cex.lab=1.3,
        cex.axis=1.2,
        outline = FALSE)
#dev.off()
```

# 热图

``` r
##  heatmap.2
data0 = read.table("data/heatmap.txt", header=T)
data = data0[,-1]
rownames(data) = data0[,1]

library(gplots)
heatmap.2(as.matrix(data))

#调整参数
#library(RColorBrewer)
color = c(rep("orange",3),rep("blue",3))
heatmap.2(as.matrix(data),col = greenred(75), 
          scale = "row",dendrogram = 'both',
          key = TRUE, symkey = FALSE, density.info = "none", 
          trace = "none", cexRow = 0.5,
          main = "Heatmap",
          ColSideColors = color  #ColSideColors就按照样本读入的顺序分配颜色
) 


## pheatmap
# 第一个例子：
rt=read.table("data/heatmap.txt",sep="\t",header=T,row.names=1)         #读取文件
library(pheatmap)                                                #引用pheatmap包
annotation=read.table("data/group.txt",sep="\t",header=T,row.names=1) #读取样品属性文件
#pdf(file="figs/pheatmap1.pdf",width=7,height=7)                        #定义输出文件
pheatmap(rt, annotation = annotation)
pheatmap(rt, annotation=annotation, color = colorRampPalette(c("green", "black", "red"))(50))
pheatmap(rt, annotation=annotation, display_numbers = TRUE)
pheatmap(rt, annotation=annotation, cluster_cols = FALSE)
#dev.off()

# 第二个例子：相关性热图
a <- read.table("data/mcp.M113.035600-2.txt",
                #数据文件名
                head = TRUE,#第一行为表头
                sep = "\t",#制表符（"\t"）分割
                #为TRUE时，进行列名检查，如果有空格
                #会用" "代替，为FALSE时，不进行列名检查
                check.names = FALSE)
#第一列为基因ID，去掉第一列，
#把剩余的列保存到x里
x <- a[,c(-1,-ncol(a))]
#对x进行log转换
x <- log2(x+1)
#计算相关系数
corX <- cor(x)

#border_color=NA 表示去掉边框颜色
pheatmap(corX,border_color=NA)

## ComplexHeatmap ：
# https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
# https://jokergoo.github.io/ComplexHeatmap-reference/book/
# BiocManager::install("ComplexHeatmap")


library(ComplexHeatmap)


# 单细胞拷贝数的热图

# 读入数据
filename <- "data/copy_number_data.txt"

my_data <- read.table(filename, sep="\t", quote="", stringsAsFactors=FALSE,header=TRUE)

head(my_data)

dim(my_data) # (rows columns)

nrow(my_data) # rows: locations (bins) in genome
ncol(my_data) # columns: cells

# 转换为矩阵
my_matrix <- as.matrix(my_data[  ,c(4:100)]) # [all rows, columns 4-100]
# 前三行，绘制热图时不需要

# 检查数据的类型
class(my_data)
class(my_matrix)

head(my_matrix)

# 提取注释信息
chromosome_info <- data.frame(chrom = my_data$CHR)
chromosome_info



# 绘制热图

# 默认参数
Heatmap(my_matrix)

# 转置
my_matrix <- t(my_matrix)  # "transpose"
Heatmap(my_matrix)

# 保持基因组顺序，不要聚类
Heatmap(my_matrix, cluster_columns=FALSE)

fontsize <- 0.6

# 将标签放左边
Heatmap(my_matrix, cluster_columns=FALSE,
        row_names_side = "left", 
        row_dend_side = "left",
        row_names_gp=gpar(cex=fontsize))

# 不同的聚类计算方法
# "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
# 默认euclidean

# 不同聚类方法
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)

# 改变计算和聚类的方法
Heatmap(my_matrix, 
        cluster_columns=FALSE,
        row_names_side = "left", 
        row_dend_side = "left",
        row_names_gp=gpar(cex=fontsize),
        row_dend_width = unit(3, "cm"),
        clustering_distance_rows ="maximum",
        clustering_method_rows = "ward.D")


# 改变聚类的颜色

# 安装dendextend
library(dendextend)
# color_brances()
# 1. 计算距离 (method="maximum")
# 2. 聚类 (method="ward.D")
dend = hclust(dist(my_matrix,method="maximum"),method="ward.D")

Heatmap(my_matrix, 
        cluster_columns=FALSE,
        row_names_side = "left", 
        row_dend_side = "left",
        row_names_gp=gpar(cex=fontsize),
        cluster_rows = color_branches(dend, k = 3))


# 按聚类切割热图

Heatmap(my_matrix, 
        cluster_columns=FALSE,
        row_names_side = "left", 
        row_dend_side = "left",
        row_names_gp=gpar(cex=fontsize),
        row_dend_width = unit(3, "cm"),
        clustering_distance_rows ="maximum",
        clustering_method_rows = "ward.D",
        km=2) # number of clusters you want



# 加载额外的注释信息
chromosome_info

chromosome.colors <- c(rep(c("black","white"),11),"red")
chromosome.colors

names(chromosome.colors) <- paste("chr",c(seq(1,22),"X"),sep="")
chromosome.colors


Heatmap(my_matrix, 
        cluster_columns=FALSE,
        row_names_side = "left", 
        row_dend_side = "left",
        row_names_gp=gpar(cex=fontsize),
        clustering_distance_rows ="maximum",
        clustering_method_rows = "ward.D",
        km=2, # 想要的聚类数目
        bottom_annotation = HeatmapAnnotation(df = chromosome_info,col = list(chrom=chromosome.colors),show_legend=FALSE)
        )
```

# 韦恩图

``` r
# 韦恩图
listA <- read.csv("data/genes_list_A.txt",header=FALSE)
A <- listA$V1
A

listB <- read.csv("data/genes_list_B.txt",header=FALSE)
B <- listB$V1
B

listC <- read.csv("data/genes_list_C.txt",header=FALSE)
C <- listC$V1
C

listD <- read.csv("data/genes_list_D.txt",header=FALSE)
D <- listD$V1
D

length(A)
length(B)
length(C)
length(D)

# VennDiagram绘制韦恩图
library(VennDiagram)

# 直接保存图片到本地

venn.diagram(list("list C"=C, "list D"=D), fill = c("yellow","cyan"), cex = 1.5, filename="Venn_diagram_genes_2.png")

venn.diagram(list(A = A, C = C, D = D), fill = c("yellow","red","cyan"), cex = 1.5,filename="Venn_diagram_genes_3.png")

venn.diagram(list(A = A, B = B, C = C, D = D), fill = c("yellow","red","cyan","forestgreen"), cex = 1.5,filename="Venn_diagram_genes_4.png")

# 在线工具
# BioVenn: http://www.cmbi.ru.nl/cdd/biovenn/index.php
# VENNY: http://bioinfogp.cnb.csic.es/tools/venny/index.html
```

## 韦恩图venn_upset

``` r
rm(list = ls())

# Cai, H., Chen, H., Yi, T., Daimon, C. M., Boyle, J. P. et al PLoS ONE, 8(1), e53388. http://doi.org/10.1371/journal.pone.0053388
# E53388, Alicia A Midland et al Cell Research(2012) 22(4)：620–623.
# http://doi.org/10.1038/cr.2012.25
# BioVenn: http://www.cmbi.ru.nl/cdd/biovenn/index.php
# VENNY: http://bioinfogp.cnb.csic.es/tools/venny/index.html
#
# 最多四个集合，不可以根据大小调整面积
#
# Venn: http://bioinformatics.psb.ugent.be/webtools/Venn/
#
#   1. 上传绘图数据：文件仅有一列，为元素名称，
# 如基因（或蛋白）名称等
# 2. 如有更多集合，可以点击添加文件，可以绘制大于
# 5个集合的文氏图
# 3. 绘图数据设好后，点击Submit
# Venn Diagram Plotter:
#   https://omics.pnl.gov/software/venn-diagram-plotter
# 最多三个集合，可以根据大小调整面积
# VennDiagram: R packages,
# https://cran.r-project.org/web/packages/VennDiagram/index.html
# Vennerabl: R packages, http://r-forge.r-project.org/R/?group_id=474

##### 韦恩图venn_upset #####

##### plotrix #####
library(plotrix)
x=0:10
y=seq(0,10,length=11)

plot(x,y,type="n",axes=F,xlab="",ylab="")
draw.circle(2,5,2,col=rgb(154/255,0/255,205/255,0.6))
draw.circle(4,5,2,col=rgb(21/255,3/255,252/255,0.6))
text(1,5,labels="10.12%",col="white",font=2,cex=0.3)
text(5,5,labels="40.38%",col="white",font=2,cex=0.3)
text(2,5,labels="49.50%",col="white",font=2,cex=0.3)
legend(6.2,5,pch=15,xjust=0,yjust=0.5,bty="n",cex=0.3,
       col=c(rgb(154/255,0/255,205/255),rgb(74/255,2/255,233/255),rgb(21/255,3/255,252/255)),
legend=c("sample 1 uniq","sample 1 & sample 2","sample 2 uniq"))
text(3.5,6.5,labels="venn chart for uniq_sRNAs",font=2,cex=0.5)


##### VennDiagram #####
files=c("inst/Venn/G1.txt","inst/Venn/G2.txt",
        "inst/Venn/G3.txt","inst/Venn/G4.txt",
        "inst/Venn/G5.txt")          #文件列表
#files=dir()
#files=grep("txt",files,value=T)
targetList=list()
#循环读取文件
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile,header=F)
  header=unlist(strsplit(inputFile,"\\.|\\-"))
  targetList[[header[1]]]=as.vector(rt[,1])
  uniqLength=length(unique(as.vector(rt[,1])))
  print(paste(header[1],uniqLength,sep=" "))
}


library(VennDiagram, quietly = TRUE)
library(grid)#引用包
grob.list <- venn.diagram(targetList,
             filename="inst/Venn/venny.tiff",
             imagetype = "tiff",
             height = 3000,
             width = 3000,
             main="", main.cex = 2,
             fill=rainbow(length(targetList)),
             cat.cex=0.8)
#grid.draw(grob.list)

intersectGenes=Reduce(intersect,targetList)

#输出交集基因
#输出交集基因
# upGenes=cbind(upGenes,"up")
# downGenes=cbind(downGenes,"down")
# intersectGenes=rbind(upGenes,downGenes)
# colnames(intersectGenes)=c("Gene","Regulation")
write.table(file="inst/Venn/target.xls",
            intersectGenes,sep="\t",
            quote=F,col.names=F,row.names=F)


####first generate the test data
G1 = read.delim("inst/Venn/venn1.txt")
Group1 = G1$V1
G2 = read.delim("inst/Venn/venn2.txt")
Group2 = G2$V1
G3 = read.delim("inst/Venn/venn3.txt")
Group3 = G3$V1
G4 = read.delim("inst/Venn/venn3.txt")
Group4 = G4$V1

#four groups
# venn.diagram(list("Group1"=Group1,
#                   "Group2"=Group2,
#                   "Group3"=Group3,
#                   "Group4"=Group4) ,
#              height=5000,
#              width=5200,
#              resolution=500,
#              imagetype="tiff",
#              filename="venn.diagram.4groups2.tiff",
#              col="transparent",
#              fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),alpha = 0.50,
#              label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
#              cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
#              cex=1.5,
#              cat.cex=1.4
# )

##### 更多实例
####first generate the test data
#a function to generate gene names
generateGeneName <- function()
{
  paste(sample(LETTERS,5,replace=T), collapse='')
}

#generate 100 gene names
geneNames <- replicate(1000, generateGeneName())


####method 1
####
library(gplots)
#venn plot of two groups
GroupA <- sample(geneNames,400, replace=F)
GroupB <- sample(geneNames,400, replace=F)
input <- list(GroupA,GroupB)
venn(input)
#venn plot of three groups
GroupA <- sample(geneNames,400, replace=F)
GroupB <- sample(geneNames,400, replace=F)
GroupC <- sample(geneNames,500, replace=F)
input <- list(GroupA,GroupB,GroupC)
venn(input)
#venn plot of four groups
GroupA <- sample(geneNames,400, replace=F)
GroupB <- sample(geneNames,300, replace=F)
GroupC <- sample(geneNames,500, replace=F)
GroupD <- sample(geneNames,700, replace=F)
input <- list(GroupA,GroupB,GroupC,GroupD)
venn(input)
#venn plot of five groups
GroupA <- sample(geneNames,400, replace=F)
GroupB <- sample(geneNames,300, replace=F)
GroupC <- sample(geneNames,500, replace=F)
GroupD <- sample(geneNames,700, replace=F)
GroupE <- sample(geneNames,600, replace=F)
input <- list(GroupA,GroupB,GroupC,GroupD,GroupE)
venn(input)

####method2
####using package "venneuler", with weights. using the size of the area to indicate different sets and the intersection
library(venneuler)
m <- data.frame(c(GroupA, GroupB, GroupC),
                c(rep("GroupA", length(GroupA)),
                  rep("GroupB", length(GroupB)),
                  rep("GroupC", length(GroupC))
                )
)
v <- venneuler(m)
plot(v)

####method3
####
library(VennDiagram)
#three groups
venn.diagram(list(GroupA=sample(geneNames,400, replace=F),
                  GroupB=sample(geneNames,300, replace=F),
                  GroupC=sample(geneNames,500, replace=F) ),
             height=3000,
             width=3000,
             resolution=500,
             imagetype="tiff",
             filename="venn.diagram.3groups.tiff",
             fill = c("cornflowerblue", "green", "yellow"),
             label.col="black",
             cat.col = c("darkblue", "darkgreen", "orange")
)
#four groups
venn.diagram(list(GroupA=sample(geneNames,400, replace=F),
                  GroupB=sample(geneNames,300, replace=F),
                  GroupC=sample(geneNames,500, replace=F),
                  GroupD=sample(geneNames,700, replace=F) ),
             height=3000,
             width=3000,
             resolution=500,
             imagetype="tiff",
             filename="inst/Figures/venn.diagram.4groups.tiff",
             col="transparent",
             fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),alpha = 0.50,
             label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
             cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4")
)
#five groups
venn.diagram(list(GroupA=sample(geneNames,400, replace=F),
                  GroupB=sample(geneNames,300, replace=F),
                  GroupC=sample(geneNames,500, replace=F),
                  GroupD=sample(geneNames,700, replace=F),
                  GroupE=sample(geneNames,800, replace=F)),
             height=3000,
             width=3000,
             resolution=500,
             imagetype="tiff",
             filename="inst/Figures/venn.diagram.5groups.tiff",
             col = "black",
             fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
             alpha = 0.50,
             cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
                     1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
             cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
             cat.cex = 1.5,
             cat.fontface = "bold",
             margin = 0.1
)

##### overLapper #####
source("code/overLapper.R")
##下面这句命令获得每个交集的ID
OLlist <- overLapper(setlist=targetList,
                     sep="_",type="vennsets")
count <- list(sapply(OLlist$Venn_List,length))
vennPlot(counts=count,mysub="",ccol=c(rep(1,30),2),
         lcex=1.5,
         ccex=c(rep(1.5,5),rep(0.6,25),1.5))


##### 韦恩图UpSetR #####
library(UpSetR, quietly = TRUE)
#读取文件
load(file = "data/UpSetData.rda")
rt=UpSetData
gene=sapply(strsplit(rownames(rt),"\\|"),"[",1)
asType=sapply(strsplit(rownames(rt),"\\|"),"[",3)
upsetList=list(AA=unique(gene[asType=="AA"]),
               AD=unique(gene[asType=="AD"]),
               AP=unique(gene[asType=="AP"]),
               AT=unique(gene[asType=="AT"]),
               ES=unique(gene[asType=="ES"]),
               ME=unique(gene[asType=="ME"]),
               RI=unique(gene[asType=="RI"]) )
upsetData=fromList(upsetList)

# pdf(file="figures/upset.pdf",
#     onefile = FALSE,width=9,height=6)              #保存图片
upset(upsetData,
      nsets = 7,                                    #展示可变剪切类型个数
      nintersects = 50,                             #展示基因集数目
      order.by = "freq",                            #按照数目排序
      show.numbers = "yes",                         #柱状图上方是否显示数值
      number.angles = 20,                           #字体角度
      point.size = 1.5,                             #点的大小
      matrix.color="red",                           #交集点颜色
      line.size = 0.8,                                #线条粗线
      mainbar.y.label = "Gene Intersections",
      sets.x.label = "Set Size")
# dev.off()
```

# 气泡图

``` r
rm(list = ls())
rt <- read.table("data/bubble_plot.txt",header=T,sep="\t")           #读取文件

library(ggplot2)                                         #引用ggplot2这个包
# 画图
p <- ggplot(rt,aes(Ratio,Term))

# 四维数据的展示
pbubble = p + geom_point(aes(size=Count,color=-1*log10(FDR)))

# 绘制pathway富集散点图，可以绘制各种颜色，只要将green和red改成自己喜欢的颜色就行
pr = pbubble + 
  scale_colour_gradient(low="green",high="red") + 
  labs(color=expression(-log[10](FDR)),
       size="Gene number",
       x="Gene ratio",
       y="Term")

ggsave("figs/bubble.pdf",width=9,height=7)                  #保存图片
```

# 饼图

``` r
rm(list = ls())
#read in data
data = read.table("data/pieplot.txt", sep="\t", header=T)
data

#draw pie chart using function "pie"
pie(data$Reads, labels=data$Type)
pie(data$Reads, labels=data$Type, radius=0.4)
pie(data$Reads, labels=data$Type, radius=0.8, clockwise=T)
pie(data$Reads, labels=data$Type, radius=0.8, clockwise=T, init.angle=90)
pie(data$Reads, labels=data$Type, radius=0.8, clockwise=T, init.angle=90, density=20)
pie(data$Reads, labels=data$Type, radius=0.8, clockwise=T, init.angle=90, density=20, col=rainbow(8))
pie(data$Reads, labels=data$Type, radius=0.8, clockwise=T, init.angle=90, density=20, col=rainbow(8), border="black")
pie(data$Reads, labels=data$Type, radius=0.8, clockwise=T, init.angle=90, density=NULL, col=rainbow(8), border="black", lty=2)
pie(data$Reads, labels=data$Type, radius=0.8, clockwise=T, init.angle=90, density=NULL, col=rainbow(8), border="black", lty=2, main="My First Pie Chart using R")

#draw a perfect par chart with optimized parameters
par(mar=c(0,6,6,6))
pie(data$Reads, 
    labels=paste(data$Type,"(",substring(data$Reads,0,4), ")"), 
    radius=0.8, 
    clockwise=T, 
    init.angle=9, 
    density=NULL, 
    col=rainbow(8), 
    border="black", 
    lty=2, 
    main="Fig1. XXX Reads of different conditions"
)

#draw a round rainbow
par(mar=c(0,0,0,0))
#first start with 10
pie(rep(1,10), col=rainbow(10), lty=0, labels='', init.angle=90, border=NA)
#what will happen if set the number of the Readss to a big one, try 200 or even bigger one 2000
#Let's see the charm of R language
par(mar=c(0,0,0,0))
pie(rep(1,2000), col=rainbow(2000), lty=0, labels='', init.angle=90, border=NA)

## 3d pie
data = read.table("data/pieplot.txt", sep="\t", header=T)

#计算比例
ratio=sprintf("%.2f",100*data[,2]/sum(data[,2]))
ratio=paste(ratio,"%",sep="")
label=paste(data[,1],ratio,sep="\n")
#载入R包
library(plotrix)
#保存为pdf
#pdf("figs/pie3D.pdf",width=4,height=4)
#设置全局图形参数
#par(mgp=c(1.6,0.6,0),mar=c(3,3,1,1))
#绘制3维饼图
pie3D(data[,2],col=1:6,
      main="Pie chart for miRNAs uniq",
      border="purple",#设置边框颜色
      labels=label,#设置标签
      labelcex=0.8,#标签字体的大小
      explode=0,#各部分分开比例，默认为0
      radius=1)#半径
#dev.off()

## 3D饼图2
rt = read.table("data/pieplot.txt", sep="\t", header=T)

x=rt[,2]
labels=as.character(rt[,1])

library("plotrix")                                               #引用包
piepercent<- paste(round(100*x/sum(x), 2), "%")                  #求各部分百分率

#pdf(file = "figs/pie3D2.pdf",width=7,height=7)                         #保存图片
pie3D(x,labels = piepercent,explode =0.1,thete=0.8)             #explode:离中心距离    theta:观看角度          
legend("topright", labels, cex = 1.1,fill = rainbow(length(x)))
#dev.off()
```

# 曼哈顿图

``` r
rm(list = ls())
rt=read.table("data/manhattanplot.txt",header=T)                             

library(qqman)                                                  
#pdf(file = 'figs/manhattan.pdf', width=15,height=7)  
manhattan(rt, 
          chr="CHR", 
          bp ="BP", 
          p="P", 
          snp="SNP", 
          col=c("blue4", "orange3"), 
          ylim = c(0, 10),
          suggestiveline = F, 
          genomewideline = -log10(3e-06))
#dev.off()
```

## 曼哈顿图manhattan&QQPlot

``` r
##### 绘制曼哈顿manhattan图 #####
library(qqman, quietly = TRUE)
load(file = "data/manhattanData.rda")
rt=manhattanData

manhattan(rt,
          chr="CHR",
          bp ="BP",
          p="P",
          snp="SNP",
          col=c("blue4", "orange3"),
          ylim = c(0, 20),
          suggestiveline = F,
          genomewideline = -log10(0.01))

##### 曼哈顿图实例 #####
library(qqman)
#have a glance at the GWAS result
str(gwasResults)
head(gwasResults)
unique(gwasResults$CHR)
#see how many SNPs on each chromosome?
as.data.frame(table(gwasResults$CHR))

##Creating manhattan plots
manhattan(gwasResults)
####Changing params
# Add a title (main=), increase the y-axis limit (ylim=),
# reduce the point size to 60% (cex=), and reduce the font
# size of the axis labels to 90% (cex.axis=). While we're
# at it, let's change the colors (col=), remove the suggestive
# and genome-wide significance lines, and supply our own
# labels for the chromosomes
####
manhattan(gwasResults, main = "Manhattan Plot", ylim = c(0, 10), cex = 0.6,
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,
          chrlabs = paste("Chr", c(1:22), sep = '') )
#changing colors
my.man.color = c('#F05859','#5862AB','#C36DAA','#F8ADAF','#808080','#BE262B',
                 '#314399','#32B34D','#BDBF36','#A44297','#21BBBD','#3E3F3F',
                 '#EF4747','#4E5BA7','#BE62A5','#BFBFBF','#7F1B1D','#3F427F',
                 '#0B8143','#7E8037','#782D7C','#008081')
manhattan(gwasResults, main = "Manhattan Plot", ylim = c(0, 10), cex = 0.6,
          cex.axis = 0.9, col = my.man.color, suggestiveline = F, genomewideline = F,
          chrlabs = paste("Chr", c(1:22), sep = '') )
#see a single chromosome
manhattan(subset(gwasResults, CHR == 1), xlim = c(0,0.001))

#highlight SNPs of Interest
str(snpsOfInterest)
manhattan(gwasResults, highlight = snpsOfInterest)
manhattan(subset(gwasResults, CHR == 3),
          highlight = snpsOfInterest,
          xlim = c(0, 0.002),
          main = "Chr 3")

#the manhattan function can be used to plot any value, not just p-values
#simply call the function passing to the p= argument the name of the column
#we want to plot instead of the default “P” column.
gwasResults <- transform(gwasResults, zscore = qnorm(P/2, lower.tail = F))
head(gwasResults)
manhattan(gwasResults, p = "zscore", logp = FALSE,
          ylab = "Z-score", genomewideline = FALSE,
          suggestiveline = FALSE,
          main = "Manhattan plot of Z-scores")


####Creating Q-Q plots
qq(gwasResults$P)
#change other graphical paras
qq(gwasResults$P, main = "Q-Q plot of GWAS p-values",
   xlim = c(0, 7), ylim = c(0, 12), pch = 18, col = "blue4",
   cex = 1.5, las = 1)

##### Q-Q_plot #####
#生成100个正态分布的随机数
data = rnorm(100,0,1)

#生成柱状图，查看所生成随机数的分布情况
hist(data,10)
#利用经验累积分布函数（empirical cumulative distribution function, eCDF)
#来评价数据分布的情况
data.cdf = ecdf(data)
plot(data.cdf)


#QQ-plot 正态分布
#50 numbers
data50 = rnorm(100,0,1)
hist(data50,20)
qqnorm(data50, pch=16)
qqline(data50, pch=16, col="red")
#100 numbers
data100 = rnorm(100,0,1)
hist(data100,20)
qqnorm(data100, pch=16)
qqline(data100, pch=16, col="red")
#150 numbers
data150 = rnorm(150,0,1)
hist(data150,20)
qqnorm(data150, pch=16)
qqline(data150, pch=16, col="red")
#1000 numbers
data1000 = rnorm(1000,0,1)
hist(data1000,20)
qqnorm(data1000, pch=16)
qqline(data1000, pch=16, col="red")

#put the figures all together to see the change and the differents
nf <- layout(matrix(c(1,2,3,4),2,2,byrow=TRUE), c(2,2), c(2,2), TRUE)
layout.show(nf)

#fig1
qqnorm(data50, pch=16, main="QQ-plot for 100 numbers")
qqline(data50, pch=16, col="red")
#fig2
qqnorm(data100, pch=16, main="QQ-plot for 200 numbers")
qqline(data100, pch=16, col="red")
#fig3
qqnorm(data150, pch=16, main="QQ-plot for 500 numbers")
qqline(data150, pch=16, col="red")
#fig4
qqnorm(data1000, pch=16, main="QQ-plot for 1000 numbers")
qqline(data1000, pch=16, col="red")

##QQ plot for two data sets
data1 = runif(200,-3,3)
data2 = rnorm(1000,0,1)
qqplot(data1, data2, pch=16, main="QQ-plot")
qqline(data1, col="blue", lwd=2, lty=2)
qqline(data2, col="red", lwd=2)

#Shapiro-Wilk test of normality
shapiro.test(data)
```

# 进化树图

``` r
rm(list = ls())

# ggtree: a phylogenetic tree viewer for different types of tree annotations
# http://www.bioconductor.org/packages/3.1/bioc/html/ggtree.html
# http://www.bioconductor.org/packages/3.1/bioc/vignettes/ggtree/inst/doc/ggtree.html
# http://www.bioconductor.org/packages/release/bioc/html/ggtree.html
# http://www.bioconductor.org/packages/release/bioc/vignettes/ggtree/inst/doc/ggtree.html

# 第一个案例
a <- read.table("data/mcp.M113.035600-2.txt",
                head = TRUE,sep = "\t",
                check.names = FALSE)
#第一列为基因ID，去掉第一列，
#把剩余的列保存到x里
x <- a[,c(-1,-ncol(a))]
#转置
x <- t(x)
x <- log2(x+1)
cfit <- hclust(dist(x))
library(ape)
#png("figs/treemap.png",width = 700,height = 700,res = 120)
plot(as.phylo(cfit),type="fan",
     label.offset=0.1,
     show.node.label=TRUE,  
     no.margin=TRUE,cex=0.8)
#dev.off()
plot(as.phylo(cfit),type="un",
     label.offset=0.1,
     show.node.label=TRUE,  
     no.margin=TRUE,cex=0.8)
plot(as.phylo(cfit),type="phylogram",
     label.offset=0.1,
     show.node.label=TRUE,  
     no.margin=TRUE,cex=0.8)

# 第二个案例
library(ggplot2)                                             #引用ggplot2包
library(ggtree)                                              #引用ggtree包
library(colorspace)                                          #引用colorspace包

cls=list()
rt=read.table("data/ggtree_group.txt",sep="\t",header=T)                   #读取属性文件
for(i in 1:nrow(rt)){
  otu=as.character(rt[i,1])
  phylum=as.character(rt[i,2])
  cls[[phylum]]=c(cls[[phylum]], otu)
}
phylumNames=names(cls)
phylumNum=length(phylumNames)

tree <- read.tree("data/ggtree_input.tre")                                 #读取树文件
tree <- groupOTU(tree, cls)

#pdf(file="figs/treemap2.pdf",width=7,height=7)                          #保存图片
ggtree(tree, 
       layout="circular", 
       ladderize = FALSE, 
       branch.length="none", 
       aes(color=group)) + 
  scale_color_manual(values=c(rainbow_hcl(phylumNum+1)),
                     breaks=phylumNames, labels=phylumNames ) + 
  theme(legend.position="right") + 
  geom_text(aes(label=paste("                ",label,sep=""), 
                angle=angle+45), 
            size=2.2)
#dev.off()
```

# 多组图

``` r
# 第一个例子
# Chow et al. BMC Genomics 2011 12:449
# Assay performance on RNA from frozen 
# tissues and artificially degraded reference RNA.
a <- read.table("data/mcp.M113.035600-2.txt",#数据文件名
                head = TRUE,#第一行为表头
                sep = "\t",#制表符（"\t"）分割
                check.names = FALSE)
#第一列为基因ID，去掉第一列，
#把剩余的列保存到x里
x <- a[,-1]
#取前四个样品为例
x <- x[,1:4]
x <- log2(x)
#######绘制简单的pair图
pairs(x)

## 第二种
dat = read.table("data/multi_plots.txt",sep="\t",header=T,fill=T)
max(dat$Count)
max(dat$logpvalue)
pvaluemax = 5
countmax = 30
par(mar=c(5,15,5,2))
barplot(dat$logpvalue[20:1],horiz = T,xlim=c(0,pvaluemax),xlab="-Log(Pvalue)",main="Biological Process",col="#1FA5FF")
axis(2,at=seq(0.8,23.5,22.7/19),labels = dat$Term[20:1],las=2,cex=1)
points(dat$Count[20:1] / countmax * pvaluemax,seq(0.8,23.5,22.7/19),type="o",col="#FF3F8C",pch=16)
axis(3,at=seq(0,countmax,5) / countmax * pvaluemax,labels=seq(0,countmax,5))

## 第三个
library(ggplot2)
theme_set(theme_bw())
theme_update(strip.background=element_rect(colour="white"))

drugs <- c("Drug1", "Drug2", "Drug3", "Drug4")
measures <- paste(c(2,4,6,8,10,12,24,48),"hours",sep=" ")
df <- data.frame(expand.grid(drug=drugs, measure=measures))
df$value <- rexp(nrow(df))

#without facet
p <- (ggplot(df, aes(x=measure, y=value, color=drug, group=drug)) +
        scale_color_discrete() +
        geom_point() + 
        geom_line(alpha=0.5) +
        theme(axis.text.x=element_text(angle=60,hjust=1,size=10))+
        xlab(""))

p
#with facet
p <- (ggplot(df, aes(x=measure, y=value, color=drug, group=drug)) +
        scale_color_discrete() +
        geom_point() + 
        facet_wrap(~ drug) +
        geom_line(alpha=0.5) +
        theme(axis.text.x=element_text(angle=60,hjust=1,size=10))+
        xlab(""))

p
```

# 物种解剖图gganatogram

``` r
# https://github.com/jespermaag/gganatogram
# https://jespermaag.shinyapps.io/gganatogram/
# install.packages("ggpolypath")

# devtools::install_github("jespermaag/gganatogram")

# runGitHub( "gganatogram", "jespermaag",
#            subdir = "shiny")

library(gganatogram)
library(dplyr)
library(viridis)
library(gridExtra)
load(file = "data/gganatogramData.rda")
# pdf("figures/gganatogramMale.pdf",width=8,height=6)                                  #男性解剖图保存文件
#读取输入文件
gganatogram(data=gganatogramMaleData, fillOutline='white',
            organism='human', sex='male', fill="value")+
  theme_void() +
  scale_fill_gradient2(low = "green", mid="black",
                       high = "red",midpoint =3)
# dev.off()

# pdf("figures/gganatogramFemale.pdf",width=8,height=6)                                #女性解剖图保存文件
#读取输入文件
gganatogram(data=gganatogramFemaleData, fillOutline='white', organism='human', sex='female', fill="value")+
  theme_void() +
  scale_fill_gradient2(low = "green", mid="black",high = "red",midpoint =3)
# dev.off()
```

# 富集分析GO_KEGG及可视化

``` r
# ==========================================================
#
#      富集分析
#      •   GO/KEGG富集分析
#      •   GSEA富集分析
#      •   多种可视化图形
#
# ==========================================================

##### 在线网站 #####
### The Database for Annotation, Visualization and Integrated Discovery (DAVID ):
# https://david.ncifcrf.gov/

### KOBAS网站
# http://kobas.cbi.pku.edu.cn/kobas3

### WebGestalt一站式工具网站：
# http://www.webgestalt.org/

### Enrichr网页的链接为：
# http://amp.pharm.mssm.edu/Enrichr

### Enrichr网页的链接为：
# http://amp.pharm.mssm.edu/Enrichr

### TAM 2.0网站：
# http://www.lirmed.com/tam2/

# 清空环境变量
rm(list = ls())

##### R包 #####
### clusterProfiler包
# http://yulab-smu.top/clusterProfiler-book/

# 加载软件包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 读入数据
load(file = "data/GOgeneList.rda")
rt <- GOgeneList
rt=rt[is.na(rt[,"entrezID"])==F,]

cor=rt$cor
# 差异倍数
# geneFC=rt$logFC
# p值
# pvalue=rt$pvalue
# 数据库数目
# sum=rt$Sum
gene=rt$entrezID
names(cor)=gene

# GO富集分析
GOResultfile <- "inst/Results/ego_GO_Result.rda"
if(!file.exists(GOResultfile)){
  ego <- enrichGO(gene = gene,
                  OrgDb = org.Hs.eg.db,
                  pvalueCutoff =0.05,
                  qvalueCutoff = 0.05,
                  ont="all",
                  readable =TRUE)
  save(ego, file = GOResultfile)
}else{
  load(file = GOResultfile)
}

# 写出到文件
#write.table(ego, file="data/GO.txt", sep="\t", quote=F, row.names = F)

#柱状图
#tiff(file="figures/barplot.tiff",width = 26,height = 20,units ="cm",compression="lzw",bg="white",res=600)
barplot(ego, drop = TRUE, showCategory =10,
        split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
#dev.off()

#气泡图
#tiff(file="figures/dotplot.tiff",width = 26,height = 20,units ="cm",compression="lzw",bg="white",res=600)
dotplot(ego,showCategory = 10,split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
#dev.off()

#热图
#tiff(file="figures/heatplot.tiff",width = 40,height = 20,units ="cm",compression="lzw",bg="white",res=600)
heatplot(ego,showCategory =20, foldChange=cor)
#dev.off()


# kegg富集分析
KEGGResultfile <- "inst/Results/kk_KEGG_Result.rda"
if(!file.exists(KEGGResultfile)){
  kk <- enrichKEGG(gene = gene,
                   organism = "hsa",
                   pvalueCutoff =0.05,
                   qvalueCutoff =0.05)
  save(kk, file = KEGGResultfile)
}else{
  load(file = KEGGResultfile)
}

barplot(kk, drop = TRUE, showCategory = 20)
cnetplot(kk, categorySize = "geneNum", foldChange = cor)

#通路图
library(pathview)
load(file = "data/KEGG_for_pathview.rda")
keggxls=KEGG_for_pathview
# for(i in keggxls$ID){
#   pv.out <- pathview(gene.data = cor,
#                      pathway.id = i,
#                      species = "hsa",
#                      out.suffix = "pathview",
#                      kegg.dir = "results/")
# }
#
# for(i in keggxls$ID){
#   pv.out <- pathview(gene.data = -log10(pvalue),
#                      pathway.id = i,
#                      species = "hsa",
#                      out.suffix = "pathview",
#                      limit=list(gene=10, cpd=10),
#                      bins = list(gene = 10, cpd= 10),
#                      kegg.dir = "results/")
# }

### GOplot包
library(GOplot, quietly = TRUE)
# GOplot数据准备
# 从clusterProfiler结果中读取
# go=data.frame(Category = ego$ONTOLOGY,ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)
# kegg=data.frame(Category = "KEGG",ID = kk$ID,Term = kk$Description, Genes = gsub("/", ", ", kk$geneID), adj_pval = kk$p.adjust)

# genelist <- rt[,1:2]
# colnames(genelist) <- c("ID", "logFC")
# circ <- circle_dat(kegg, genelist)

# 从文件中读取
load(file = "data/GOplotData.rda")
david <- david_for_GOplot
genelist <- gene_for_GOplot
circ <- circle_dat(david, genelist)

GOplot1 <- GOBar(circ)
print(GOplot1)
# GOBar(subset(circ, category == 'BP'))
# ggsave(GOplot1, filename = "figures/GOplot1.png")

GOplot2 <- GOBar(circ, display = 'multiple')
print(GOplot2)
# ggsave(GOplot2, filename = "figures/GOplot2.png")

GOplot3 <- GOBar(circ, display = 'multiple',
                 title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))
print(GOplot3)
# ggsave(GOplot3, filename = "figures/GOplot3.png")

GOplot4 <- GOBubble(circ,labels = 1)
print(GOplot4)
GOBubble(circ,
         title = 'Bubble plot',
         colour = c('orange', 'darkred', 'gold'),
         display = 'multiple', labels = 1)
GOBubble(circ,
         title = 'Bubble plot with background colour',
         display = 'multiple',
         bg.col = T, labels = 1)
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
GOBubble(reduced_circ, labels = 1)

# ggsave(GOplot4, filename = "figures/GOplot4.png", width = 15, height = 8)

GOplot5 <- GOCircle(circ, nsub = 6)
print(GOplot5)
# ggsave(GOplot5, filename = "figures/GOplot5.png", width = 10, height = 8)

# tiff(file="figures/GOplot5-GOCircle.tiff",width = 20,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOCircle(circ,table.legend = F,label.size=5,nsub=nrow(david))
# dev.off()

chord <- chord_dat(circ,genelist,david$Term)
GOplot6 <- GOChord(chord, space = 0.02, gene.order = 'logFC',gene.space = 0.25, gene.size = 3)
print(GOplot6)
# ggsave(GOplot6, filename = "pics/GOplot6.png", width = 13, height = 13)

#限定term数目
termNum <- 5
#限定基因数目
geneNum <- nrow(genelist)
chord <- chord_dat(circ, genelist[1:geneNum,],
                   david$Term[1:termNum])
# tiff(file="figures/GOplot6-chord.tiff",width = 32,height = 32,units ="cm",compression="lzw",bg="white",res=300)
GOChord(chord,
        ribbon.col=brewer.pal(termNum,'Set3'),
        gene.size=6)
# dev.off()
# clusterProfiler版
# chord <- chord_dat(circ,
#                    genelist[1:geneNum,],
#                    kegg$Term[1:termNum])
# GOChord(chord,
#         space = 0.001,           #基因之间的间距
#         gene.order = 'logFC',    #按照logFC值对基因排序
#         gene.space = 0.25,       #基因名跟圆圈的相对距离
#         gene.size = 4,           #基因名字体大小
#         border.size = 0.1,       #线条粗细
#         process.label = 7.5)     #term字体大小


GOplot7 <- GOCluster(circ, david$Term, clust.by = 'logFC', term.width = 2)
print(GOplot7)
# ggsave(GOplot7, filename = "figures/GOplot7.png", width = 12, height = 8)

# tiff(file="figures/GOplot7-cluster.tiff",width = 30,height = 30,units ="cm",compression="lzw",bg="white",res=300)
GOCluster(circ, david$Term[1:termNum])
# dev.off()

termCol <- c("#223D6C","#D20A13","#FFD121",
             "#088247","#58CDD9","#7A142C",
             "#5D90BA","#431A3D","#91612D",
             "#6E568C","#E0367A","#D8D155",
             "#64495D","#7CC767")

# GOCluster(circ,
#           go$Term[1:termNum],
#           lfc.space = 0.2,                   #倍数跟树间的空隙大小
#           lfc.width = 1,                     #变化倍数的圆圈宽度
#           term.col = termCol[1:termNum],     #自定义term的颜色
#           term.space = 0.2,                  #倍数跟term间的空隙大小
#           term.width = 1)                    #富集term的圆圈宽度
```

# 配对检验与可视化pairedPlot

``` r
#' Paired plot
#'
#' @param data data for paired plot
#' @param group group vector
#' @return paired plot
#' @author Yangming si
#' @export
PairedPlot <- function(data,group){
  xfactors = as.factor(group[,2])
  xnumsample = as.numeric(xfactors)
  xaxis = levels(xfactors)
  link = group[,3]
  links = unique(group[,3])
  x1data = data[xnumsample==1]
  x2data = data[xnumsample==2]
  pvalue = round(wilcox.test(x1data,x2data)$p.value,3)

  par(las=1)
  plot(1,xlim=c(0.5,2.5),ylim=c(0,max(data)*1.2),type="n",xlab="",ylab="",xaxt="n")
  points(rep(1,length(x1data)),x1data,pch=16,cex=2)
  points(rep(2,length(x2data)),x2data,pch=15,cex=2)
  axis(1,1:2,xaxis)
  for(i in links){
    w = which(link==i)
    x1 = xnumsample[w[1]]
    y1 = data[w[1]]
    x2 = xnumsample[w[2]]
    y2 = data[w[2]]
    segments(x1,y1,x2,y2)
  }
  par(xpd=T)
  arrows(1,max(data)*1.1,2,max(data)*1.1,angle=90,code=3,length=0.1)
  text(1.5,max(data)*1.1,paste("pvalue =",pvalue),pos=3,cex=1)
}
```

``` r
##### 配对检验与可视化pairedPlot #####
rm(list=ls())
#设置图片输出目录的名字
picDir="inst/Figures/pairedPlotpicture"
if(!dir.exists(picDir)){
  dir.create(picDir)
}
load(file = "data/PairedPlot.rda")

m = match(group[,1],rownames(df))
df = df[m,]

#自定义画图函数
#for(i in 1:ncol(df))
for(i in 1:2){
  data = df[,i]
  cell = colnames(df)[i]
  cell = gsub(' ','_',cell)
  outtiff = paste0(picDir,"\\",cell,".tiff")
  tiff(file=outtiff,width = 14,height = 12,
       units ="cm",compression="lzw",bg="white",res=300)
  #pdf(outpdf,width=7,height=6)
  PairedPlot(data,group)
  dev.off()
}
```

# 圈图RCircos_circlize

``` r
##### 圈图circlize #####
#引用包
suppressPackageStartupMessages(library(circlize))
#读取文件1
load("data/CircosData.rda")
#pdf(file="figures/circlize_circos.pdf",width=7,height=7)
circos.par("track.height" = 0.1,
           cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram()
circos.genomicLink(bed1, bed2, col =rainbow(20), border = NA)             #rainbow(20)
circos.clear()
#dev.off()
```

# 突变分析

## 突变分析可视化maftools

``` r
##### 突变分析可视化maftools #####
#http://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html
library(TCGAmutations, quietly = TRUE)
library(maftools)
maf = read.maf(maf = 'inst/maftools.maf')

#pdf(file="results/maftools/summary.pdf",width=7,height=6)
plotmafSummary(maf = maf, rmOutlier = TRUE,
               addStat = 'median', dashboard = TRUE,
               titvRaw = FALSE)
#dev.off()

#pdf(file="results/maftools/waterfall.pdf",width=7,height=6)
oncoplot(maf = maf, top = 30, fontSize = 0.5 ,showTumorSampleBarcodes = F )
#dev.off()

#pdf(file="results/maftools/interaction.pdf",width=7,height=6)
somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.001))
#dev.off()

#pdf(file="results/maftools/Genecloud.pdf",width=7,height=6)
#source("https://gist.githubusercontent.com/PoisonAlien/3f8752a89c1d63f64afe55b441621223/raw/3959548d1a2d70843d562963dd8d9373d4d4f884/geneCloud.r")
geneCloud(input = maf, minMut = 5)
#dev.off()
```

# 瀑布图waterfall_GenVisR

``` r
load(file = "data/waterfallData.rda")
##### 瀑布图waterfall_GenVisR #####
library(GenVisR)

# pdf(file="figures/waterfall.pdf",height=9,width=12)
waterfall(waterfallData1)
# dev.off()

# Melt the clinical data into 'long' format.
library(reshape2)
clinical <- melt(waterfallclinical, id.vars = c("sample"))

# Run waterfall
# pdf(file="waterfall2.pdf",height=10,width=12)
waterfall(waterfallData2, clinDat = clinical,clinLegCol=3)
# dev.off()
```

# 棒棒糖图

## LollipopChart

``` r
# 读取数据
rt <- read.table("test.txt", 
                 sep = "\t",
                 header = TRUE)

dat <- plyr::arrange(rt, Correlation.Coefficient)
dat$order <- as.numeric(rownames(dat))
```

``` r
library(ggplot2)
library(glue)
library(ggtext)

#data(LollipopChart, package = "StatPlotR")
load("data/LollipopChart.rda")
dat <- plyr::arrange(LollipopChart, Correlation.Coefficient)
dat$order <- as.numeric(rownames(dat))

# 处理y轴标签
labels_y_left = dat$Type
highlight = function(x, pat, color="black", family="") {
  ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}
target <- round(dat$pvalue[dat$pvalue<0.05],3)
target_con <- c()
for(i in 1:length(target)){
  ifelse(i==1,
         target_con <- paste(sprintf("%0.3f",target[i]),sprintf("%0.3f",target[i+1]),sep = "|"),
         target_con <- paste(target_con,sprintf("%0.3f",target[i]),sep = "|"))
}
labels_y_right = highlight(sprintf("%.3f",dat$pvalue),
                           target_con,"black")


pillipopplot =
  # 构建画布
  ggplot(data = dat,
         mapping = aes(x = Correlation.Coefficient,
                       y = order)) +
  # 加辅助线
  geom_line(color = NA) +
  geom_segment(aes(x=0,xend=Correlation.Coefficient,
                   y=order,
                   yend=order)) +
  # 气泡图
  geom_point(aes(size=abs(Correlation.Coefficient),
                 color=-log10(pvalue))) +
  guides(color="none", size = guide_legend("Correlation")) +
  scale_color_gradient(low="yellow",high = "green") +
  labs(x="Correlation Coefficient(r)", y="") +
  # 构建并修改第一/二y坐标轴标签
  scale_y_continuous(breaks = 1:length(dat$order),
                     labels = labels_y_left,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(dat$order),
                                         labels = labels_y_right)) +
  # 标题和背景主题调整
  labs(title = "GRB10", subtitle = "pvalue") +
  theme(axis.ticks = element_blank(),
        panel.background = element_rect(fill=NA,size=rel(10)),
        panel.grid.major = element_line(colour = "darkgrey"),
        panel.grid.minor = element_line(colour = "lightgrey"),
        plot.title = element_text(size = 20,hjust = 0.5,margin = margin(b = -2)),
        plot.subtitle = element_text(size = 15,hjust = 1.18,margin = margin(t = -0.5,unit = "cm")),
        axis.text.y.right = element_markdown(margin = margin(l = 0.5,t = 0.2,unit = "cm")),
  )
pillipopplot
#ggsave("figures/pillipop.png")
```

## pillipop1

``` r
rt <- read.table("test.txt", 
                 sep = "\t",
                 header = TRUE)

dat <- plyr::arrange(rt, Correlation.Coefficient)
dat$order <- as.numeric(rownames(dat))
library(ggplot2)

# 生成画布
p = ggplot(data = dat,
           mapping = aes(x = Correlation.Coefficient,
                            y = order))


# 加辅助线
p1 = p + 
  geom_line(color = NA) +
  geom_segment(aes(x=0,xend=Correlation.Coefficient,
                   y=order, 
                   yend=order))

p1
# 绘制气泡
p2 = p1 + 
  geom_point(aes(size=abs(Correlation.Coefficient),
                 color=-log10(pvalue))) + 
  guides(color="none", size = guide_legend("Correlation")) + 
  scale_color_gradient(low="yellow",high = "green") +
  labs(x="Correlation Coefficient(r)", y="") 

p2

# 构建第二坐标轴标签
library(glue)
library(ggtext)
labels_y_left = dat$Type
highlight = function(x, pat, color="black", family="") {
     ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}
target <- round(dat$pvalue[dat$pvalue<0.05],3)
target_con <- c()
for(i in 1:length(target)){
  ifelse(i==1,
         target_con <- paste(sprintf("%0.3f",target[i]),sprintf("%0.3f",target[i+1]),sep = "|"),
         target_con <- paste(target_con,sprintf("%0.3f",target[i]),sep = "|"))
}
labels_y_right = highlight(sprintf("%.3f",dat$pvalue),
                           target_con,"black")


p3 = p2 +  
  scale_y_continuous(breaks = 1:length(dat$order),
                     labels = labels_y_left,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(dat$order),
                                         labels = labels_y_right))
  
  
p3

p3 + 
  labs(title = "GRB10", subtitle = "pvalue") +
  theme(axis.ticks = element_blank(),
        panel.background = element_rect(fill=NA,size=rel(10)),
        panel.grid.major = element_line(colour = "grey"),
        panel.grid.minor = element_line(colour = "lightgrey"),
        plot.title = element_text(size = 20,hjust = 0.5,margin = margin(b = -2)),
        plot.subtitle = element_text(size = 15,hjust = 1.18,margin = margin(t = -0.5,unit = "cm")),
        axis.text.y.right = element_markdown(margin = margin(l = 0.5,t = 0.2,unit = "cm")),
    )
ggsave("pillipop.png")
```