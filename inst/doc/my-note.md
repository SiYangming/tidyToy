---
title: "my-note"
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{my-note}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(tidyToy)
```

# 数据处理

## 批次标准化batchNormalize_sva

``` R
# 加载软件包
library(sva)
library(limma)

rt=read.delim("data-raw/batchNormalize.txt",sep="\t")
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
batchType=c(rep(1,50),rep(2,10))
modType=c(rep("normal",25),rep("tumor",25),rep("normal",5),rep("tumor",5))
mod = model.matrix(~as.factor(modType))
outTab=ComBat(data, batchType, mod, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab),outTab)
# save(outTab, file = "inst/Results/sva_batchNormalize.Rdata")
# write.table(outTab,file="normalize.txt",sep="\t",quote=F,col.names=F)
```

## 数据标准化normalize

``` R
##### 基因表达量标准化方法 #####

library(limma)
data=data[rowMeans(data)>0,]

### 对重复基因名取平均值
GeneExp <- avereps(data)

##### 芯片数据normalizeBetweenArrays标准化 #####

### 在芯片间进行标准化矫正数据
GeneExp <- normalizeBetweenArrays(as.matrix(GeneExp))
normalData <- cbind(id=row.names(GeneExp),GeneExp)


##### voom对转录组进行标准化 #####
#正常样品数目
normalCount=10
#肿瘤样品数目
tumorCount=10
group=c(rep("normal",normalCount),rep("tumor",tumorCount))
design <- model.matrix(~factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(GeneExp)
v <-voom(GeneExp, design = design,
         plot = F, save.plot = F)
out=v$E
# out=rbind(ID=colnames(out),out)
# write.table(out,file="uniq.symbol.txt",
#             sep="\t",quote=F,col.names=F)
```

## 数据补全impute

``` R
##### impute补全缺失值 #####
#最近邻居法（KNN，k-Nearest Neighbor）法：
#此方法是寻找和有缺失值的基因的表达谱相似的其他基因，
#通过这些基因的表达值（依照表达谱相似性加权）来填充
#缺失值

library(impute)
load(file = "data/limmaData.rda")
geneExpMatrix <- as.matrix(limmaData)

# KNN法计算缺失值
geneExpMatrix <- impute.knn(geneExpMatrix,
                            k=10,rowmax=0.5,
                            colmax=0.8,maxp=3000,
                            rng.seed=362436069)
#读出经过缺失值处理的数据
GeneExp <- geneExpMatrix$data
data=GeneExp[rowMeans(GeneExp)>0.05,]

#过滤波动太小的AS
genes=c()
for(i in rownames(data)){
  if(sd(data[i,])>0.01){
    genes=c(genes,i)
  }
}
#输出满足过滤条件的所有AS的表格
# all=rbind(ID=colnames(data[genes,]),data[genes,])
# write.table(all,file="imputeASmatrix.txt",
#             sep="\t",col.names=F,quote=F)
```

## 数据分隔caret

``` R
library(caret)

load(file = "data/caretData.rda")
set.seed(300)
data<-createDataPartition(y=TCGA$id,p=0.50,list=F)
TrainData<-TCGA[data, ]
TestData<-TCGA[-data,]
# write.table(TrainData, "results/trian_data.txt",row.names = F,quote = F,sep = "\t")
```

## 数据合并

``` R
##### intersect合并数据 #####
load(file = "data/MergeData.rda")

# 处理文件名
colnames(cnv)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*",
                   "\\1\\-\\2\\-\\3\\-\\4",colnames(cnv))
colnames(RNA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*",
                   "\\1\\-\\2\\-\\3\\-\\4",colnames(RNA))

group=sapply(strsplit(colnames(cnv),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
Noraml=cnv[,group==1]
methy=cnv[,group==0]

### with和subset取子集
### write table
# diffLab <- allDiff[with(allDiff,
#                         ((logFC>2 |logFC<(-2)) &
#                         adj.P.Val<0.01)),]
# write.table(diffLab,file="diffEXp.xls",
#             sep="\t", quote=F)

# diffUp <- allDiff[with(allDiff,
#                        (logFC>logFoldChange &
#                           adj.P.Val < adjustP )), ]
#
# diffDown <- allDiff[with(allDiff,
#                          (logFC<(-logFoldChange) &
#                             adj.P.Val < adjustP )), ]
```

## 基因名转换和处理

``` R
##### 提取编码蛋白和lncRNA #####
bioType <- read.delim("data-raw/biotype.txt", row.names = 1)
type=sapply(strsplit(rownames(bioType),"\\|"),"[",2)
protein=bioType[type=="protein_coding",]
lncRNA=bioType[type=="lncRNA",]
rownames(bioType)=gsub("(.*?)\\|.*","\\1",rownames(bioType))
rownames(protein)=gsub("(.*?)\\|.*","\\1",rownames(protein))
rownames(lncRNA)=gsub("(.*?)\\|.*","\\1",rownames(lncRNA))
rownames(bioType)=gsub("(.*?)\\|.*","\\1",rownames(bioType))

##### 基因名转换symbol2id #####
library(org.Hs.eg.db)
entrezIDs <- mget(genes,
                  org.Hs.egSYMBOL2EG,
                  ifnotfound=NA)
                  
##### ENSE2Symbol
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

ENSE2Symbol <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)

# save(ENSE2Symbol, file = "data/ENSE2Symbol.rda")
```

## 比对rBlast

``` R
library(rBLAST)
library(dplyr)

##### 学习资料 #####
## blast比对
# 从以下网站下载序列：https://www.gencodegenes.org/human/
# blast比对软件 https://www.ncbi.nlm.nih.gov/books/NBK131777/

## 需要服务器运行
# makeblastdb -in gencode.v35.transcripts.fa -dbtype nucl
# blastn -db gencode.v35.transcripts.fa -query Step05-GSE70880_GPL19748.fasta -out Step05-blastOut.txt -outfmt 6 -num_threads 3 -num_alignments 1 -evalue 1e-10

#### rBLAST比对
# https://github.com/SiYangming/rBLAST
# devtools::install_bioc("Biostrings")
# devtools::install_github("mhahsler/rBLAST")
# Install the BLAST software by following the instructions found in ?blast

##### rBLAST #####

## 构建数据库
makeblastdb("data-raw/gencodev35.fa", dbtype = "nucl", args="-out inst/examples/gencodev35/gencodev35")

## 读取注释文件

seq <- readDNAStringSet("inst/examples/GPL27713-41055.fasta")

#bl <- blast(db = "inst/examples/gencodev37_lncRNA/gencodev37_lncRNA", type = "blastn")

bl <- blast(db = "inst/examples/gencodev35/gencodev35", type = "blastn")

#cl <- predict(bl, seq,
#BLAST_args="-num_threads 3 -num_alignments 1 -evalue 1e-10")

## 进行比对
cll <- list()
cll[[1]] <- predict(bl, seq[1:30000,]
                    # ,
                    # BLAST_args="-num_threads 3 -num_alignments 1 -evalue 0.05"
)
cll[[2]] <- predict(bl, seq[30001:60000,]
                    # ,
                    # BLAST_args="-num_threads 3 -num_alignments 1 -evalue 0.05"
)
cll[[3]] <- predict(bl, seq[60001:90000,]
                    # ,
                    # BLAST_args="-num_threads 3 -num_alignments 1 -evalue 0.05"
)
cll[[4]] <- predict(bl, seq[90001:120000,]
                    # ,
                    # BLAST_args="-num_threads 3 -num_alignments 1 -evalue 0.05"
)
cll[[5]] <- predict(bl, seq[120001:length(seq),]
                    # ,
                    # BLAST_args="-num_threads 3 -num_alignments 1 -evalue 0.05"
)
cl <- do.call(rbind,cll)
write_csv(cl, file = "inst/examples/GPL27713_blast.csv")
namesDF <- as.data.frame(do.call(rbind,str_split(cl$SubjectID,"\\|")))
anno <- namesDF %>%
  transmute(ProbeID = cl$QueryID, EnsemblID = V2, geneSymbol = V5, Alias = V6, bioType = V8)
write_csv(anno, file = "inst/examples/GPL27713_anno.csv")
```

# 差异与检验

## 单多cox筛选基因

``` R
##### 单因素cox回归 #####
pFilter=0.05

library(survival)
rt=read.delim("inst/survival/survivalInput.txt",
              check.names=F,row.names=1)

coxR=data.frame()
# coxf<-function(x){
#   fmla1 <- as.formula(Surv(survival_time,vital_status)~lncRNA[,x])
#   mycox <- coxph(fmla1,data=lncRNA)
# }
# for(a in colnames(lncRNA[,3:ncol(lncRNA)])){
#   mycox=coxf(a)
#   coxResult = summary(mycox)
#   coxR=rbind(coxR,cbind(lncRNAname=a,
#                         HR=coxResult$coefficients[,"exp(coef)"],
#                         P=coxResult$coefficients[,"Pr(>|z|)"]))
# }

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
#输出所有单因素的结果
outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="inst/Results/uniCox.xls",sep="\t",row.names=F,quote=F)

#输出单因素显著的结果
sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
# write.table(sigTab,file="inst/Results/uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)
#输出单因素显著AS的PSI值，用于后续建模
sigGenes=c("futime","fustat")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
# write.table(uniSigExp,file="inst/Results/uniCoxSigExp.txt",
#             sep="\t",row.names=F,quote=F)

##### 单因素森林图 #####
#绘制森林图
library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
#定义森林图颜色
clrs <- fpColors(box="green",
                 line="darkblue",
                 summary="royalblue")
rt=read.delim("inst/survival/uniCox.xls",
              row.names=1,check.names=F)

data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <-
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
# pdf(file="inst/Figures/uniCoxforest.pdf",
#     width = 6,             #图片的宽度
#     height = 4,            #图片的高度
# )
forestplot(tabletext,
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.4,
           xlab="Hazard ratio"
)
# dev.off()


##### 用韦恩表示cox分析结果 #####
library(UpSetR)
rd <-read.delim("inst/survival/uniCoxResult.Sig.txt",
                row.names = 1)
sigGenes <- rownames(rd)
#绘制upset图
gene=sapply(strsplit(sigGenes,"\\|"),"[",1)
asType=sapply(strsplit(sigGenes,"\\|"),"[",3)
upsetList=list(AA=unique(gene[asType=="AA"]),
               AD=unique(gene[asType=="AD"]),
               AP=unique(gene[asType=="AP"]),
               AT=unique(gene[asType=="AT"]),
               ES=unique(gene[asType=="ES"]),
               ME=unique(gene[asType=="ME"]),
               RI=unique(gene[asType=="RI"]) )
upsetData=fromList(upsetList)

# pdf(file="inst/Figures/uniCoxUpset.pdf",onefile = FALSE,width=8,height=5)              #保存图片
upset(upsetData,
      nsets = 7,                                    #展示可变剪切类型个数
      order.by = "freq",                            #按照数目排序
      show.numbers = "yes",                         #柱状图上方是否显示数值
      number.angles = 20,                           #字体角度
      point.size = 1.5,                             #点的大小
      matrix.color="red",                           #交集点颜色
      line.size = 0.8,                              #线条粗线
      mainbar.y.label = "Gene Intersections",
      sets.x.label = "Set Size")
#dev.off()





##### 多因素cox回归 #####
library(survival)

rt=read.delim("inst/survival/lassoSigExp.txt",
              check.names=F,row.names=1)    #读取输入文件

#COX模型构建
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
# Concordance= 0.717  (se = 0.026 )
# C_index
#Concordance= 0.86  (se = 0.08 )
# 0.5      #没有任何预测能力
# 0.51-0.7 #低准确度
# 0.71-0.9 #中等准确度
# >0.9     #高准确度

#lower:0.86-1.96*0.08
#up: 0.86+1.96*0.08
#输出模型参数
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="results/multiCox.xls",
            sep="\t",row.names=F,quote=F)

#输出病人风险值
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
# write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
#             file="inst/Results/multiCoxrisk.txt",
#             sep="\t",
#             quote=F,
#             row.names=F)

##### 多因素森林图 #####
library(survival)
library(survminer)

# pdf(file="inst/Figures/multiCoxforest.pdf",
#     width = 8,             #图片的宽度
#     height = 5,            #图片的高度
# )
ggforest(multiCox)
# ggforest(multiCox,
#          main = "Hazard ratio",
#          cpositions = c(0.02,0.22, 0.4),
#          fontsize = 0.7,
#          refLabel = "reference",
#          noDigits = 2)
# dev.off()
```

# 甲基化

## 甲基化驱动基因MethylMix

``` R
##### 甲基化驱动基因MethylMix #####
adjustpFilter=0.05
logFCfilter=0
corFilter=-0.3

library(MethylMix)
load("data/MethylMixData.rda")
GEcancer=as.matrix(GeneCancer)
METcancer=as.matrix(MethylationCancer)
METnormal=as.matrix(MethylationNormal)

MethylMixResultsFile <- "inst/Results/MethylMixResults.rda"

if(!file.exists(MethylMixResultsFile)){
  library(doParallel)
  cl <- makeCluster(5)
  registerDoParallel(cl)
  MethylMixResults=
    MethylMix(METcancer, GEcancer, METnormal)
  stopCluster(cl)
  save(MethylMixResults,
       file = MethylMixResultsFile)
}else{
  load(file = MethylMixResultsFile)
}

# Found 251 samples with both methylation and expression data.
# Correlating methylation data with gene expression...
#
# Found 9 transcriptionally predictive genes.
#
# Starting Beta mixture modeling.
# Running Beta mixture model on 9 genes and on 251 samples.

MethylMixResults$MethylationStates[, 1:5]

plots <- MethylMix_PlotModel("ERBB2",
                             MethylMixResults,
                             METcancer, GEcancer,
                             METnormal)
plots$MixtureModelPlot
plots$CorrelationPlot

# for (gene in MethylMixResults$MethylationDrivers) {
#   MethylMix_PlotModel(gene, MethylMixResults, METcancer, METnormal = METnormal)
# }

outTab=data.frame()
for (gene in MethylMixResults$MethylationDrivers) {
  wilcoxTest=wilcox.test(METnormal[gene,], METcancer[gene,])
  wilcoxP=wilcoxTest$p.value
  normalGeneMeans=mean(METnormal[gene,])
  tumorGeneMeans=mean(METcancer[gene,])
  logFC=log2(tumorGeneMeans)-log2(normalGeneMeans)
  betaDiff=tumorGeneMeans-normalGeneMeans

  normalMed=median(METnormal[gene,])
  tumorMed=median(METcancer[gene,])
  diffMed=tumorMed-normalMed

  x=as.numeric(METcancer[gene,])
  y=as.numeric(GEcancer[gene,])
  corT=cor.test(x,y)
  z=lm(y~x)
  cor=corT$estimate
  cor=round(cor,3)
  pvalue=corT$p.value

  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
    outTab=rbind(outTab,cbind(gene=gene,normalMean=normalGeneMeans,TumorMean=tumorGeneMeans,logFC=logFC,pValue=wilcoxP,cor=corT$estimate,corPavlue=corT$p.value))
  }
}
pValue=outTab[,"pValue"]
adjustP=p.adjust(as.numeric(as.vector(pValue)),method="bonferroni")
outTab2=cbind(outTab,adjustP=adjustP)
out=outTab2[(outTab2$adjustP<adjustpFilter & abs(as.numeric(as.vector(outTab2$logFC)))>logFCfilter & as.numeric(as.vector(outTab2$cor))<corFilter),]
out=out[c("gene","normalMean","TumorMean","logFC","pValue","adjustP","cor","corPavlue")]
out=out[order(as.numeric(as.vector(out$pValue))),]

# write.table(out,file="results/drivenGene.xls",
#             sep="\t",row.names=F,quote=F)
# write.table(out[,1],file="results/drivenGene.txt",
#             sep="\t",row.names=F,quote=F,col.names=F)


##### 批量绘制甲基化可视化结果 #####
methyDir="methy"
corDir="methycor"
if(!dir.exists(methyDir)){
  dir.create(methyDir)
}

if(!dir.exists(corDir)){
  dir.create(corDir)
}

geneList=out$gene

outTab=data.frame()
for (gene in as.vector(geneList) ) {
  x=as.numeric(METcancer[gene,])
  y=as.numeric(GEcancer[gene,])
  corT=cor.test(x,y)
  z=lm(y~x)
  cor=corT$estimate
  cor=round(cor,3)
  pvalue=corT$p.value
  if(pvalue<0.001){
    pval=signif(pvalue,4)
    pval=format(pval, scientific = TRUE)
  }else{
    pval=round(pvalue,3)}

  tiffFile=paste(gene,".tiff",sep="")
  tiff(file=paste(methyDir,tiffFile,sep="\\"),width = 30, height = 15,
       units ="cm",compression="lzw",bg="white",res=600)
  plots=MethylMix_PlotModel(gene,MethylMixResults,METcancer,GEcancer,METnormal)
  print(plots$MixtureModelPlot)
  dev.off()

  tiffFile=paste(gene,".cor.tiff",sep="")
  outTiff=paste(corDir,tiffFile,sep="\\")
  tiff(file=outTiff,width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
  plot(x,y, type="p",pch=16,main=paste("Cor=",cor," (p-value=",pval,")",sep=""),
       cex=0.9, cex.lab=1.05, cex.main=1.1,cex.axis=1,     #由左到右一次是点、坐标轴、标题、刻度字体大小
       xlab=paste(gene,"methylation"),
       ylab=paste(gene,"expression") )
  lines(x,fitted(z),col=2)
  dev.off()
}

```

# 免疫

## 甲基化免疫浸润EpiDISH

``` r
##### 甲基化免疫浸润EpiDISH #####
library(EpiDISH)
data(centDHSbloodDMC.m)
ref.m <- centDHSbloodDMC.m[,1:6]
load(file = "data/EpiDISHData.rda")
epiOut=epidish(as.matrix(EpiDISHData), ref.m,method="CBS")
out=epiOut[["estF"]]
out=rbind(id=colnames(out),out)
#write.table(out,file="results/EpiDISH.out.txt",sep="\t",quote=F,col.names=F)
```

## 肿瘤纯度打分estimate

``` r
#https://bioinformatics.mdanderson.org/estimate/rpackage.html

library(limma)
library(estimate)
load("data/CIBERSORTData.rda")

#如果一个基因占了多行，取均值
rt=as.matrix(CIBERSORTData)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
### TCGA样本名处理
# group=sapply(strsplit(colnames(data),"\\-"),"[",4)
# group=sapply(strsplit(group,""),"[",1)
# group=gsub("2","1",group)
# data=data[,group==0]
out=data[rowMeans(data)>0,]
out=rbind(ID=colnames(out),out)
#输出整理后的矩阵文件
write.table(out,file="inst/Results/estimate-uniqinput.txt",sep="\t",quote=F,col.names=F)

#运行estimate包
filterCommonGenes(input.f="inst/Results/estimate-uniqinput.txt",
                  output.f="inst/Results/commonGenes.gct",
                  id="GeneSymbol")

estimateScore(input.ds = "inst/Results/commonGenes.gct",
              output.ds="inst/Results/estimateScore.gct",
              platform="illumina")

#输出每个样品的打分
scores=read.delim("inst/Results/estimateScore.gct",skip = 2)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="inst/Results/Estimatescores.txt",
            sep="\t",quote=F,col.names=F)
```

## 独立性检验

``` R
##### 单因子独立性检验 #####

library(survival)
rt=read.delim("inst/survival/unIndepInput.txt",
              check.names=F,row.names=1)

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="results/unindep_uniCox.xls",
            sep="\t",row.names=F,quote=F)


##### 单因子独立性检验 #####
library(survival)
library(forestplot)
#定义森林图颜色
clrs <- fpColors(box="green",line="darkblue", summary="royalblue")
rt=read.delim("inst/Results/unindep_uniCox.xls",
              row.names=1,check.names=F)

data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <-
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
# pdf(file="figures/unindep_forest.pdf",onefile = FALSE,
#     width = 6,             #图片的宽度
#     height = 4,            #图片的高度
# )
forestplot(tabletext,
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Hazard ratio"
)
# dev.off()

##### 多因子森林图 #####
clrs <- fpColors(box="red",line="darkblue", summary="royalblue") #定义森林图颜色
rt=read.delim("inst/survival/unIndepInput.txt",
              check.names=F,row.names=1)
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="results/multiIndep_multiCox.xls",
            sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.delim("inst/Results/multiIndep_multiCox.xls",
              check.names=F,row.names=1)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <-
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
# pdf(file="figures/multiIndep_multiCox_forest.pdf",onefile = FALSE,
#     width = 6,             #图片的宽度
#     height = 4,            #图片的高度
# )
forestplot(tabletext,
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Hazard ratio"
)
# dev.off()
```

## 生存分析survival

``` R
##### 生存分析survival #####

##### 生存文件准备 #####
#gene='BRCA1'   #基因名（需修改）
clinicalFile="inst/survival/time.txt"
expFile="inst/survival/normalizeExp.txt"

library(hash)
rt=read.delim(clinicalFile,check.names=F)
h = hash(keys = rt$id,
         values = paste(rt$futime,rt$fustat,sep="\t"))

exp=read.delim(expFile,check.names=F,row.names=1)

group=sapply(strsplit(colnames(exp),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
exp=exp[,group==0]

#geneExp=t(exp[gene,])
geneExp=t(exp)
write.table("sample\tfutime\tfustat\texpression",
            file="inst/survival/survivalInput.txt",sep="\t"
            ,quote=F,row.names=F,col.names=F)

for(i in rownames(geneExp)){
  j=unlist(strsplit(i,"\\-"))
  if(grepl("^0",j[4])){
    name4=paste(j[1],j[2],j[3],j[4],sep="-")
    name3=paste(j[1],j[2],j[3],sep="-")
    if(has.key(name3,h)){
      write.table(paste(name4,h[[name3]],geneExp[i,],sep="\t"),
                  file="inst/survival/survivalInput.txt",sep="\t",
                  quote=F,append=TRUE,row.names=F,col.names=F)
    }
  }
}

##### 单个基因做差异分析 #####

library(survival)
rt=read.delim("inst/survival/survivaldata.txt")
rt$futime=rt$futime/365       #如果以月为单位，除以30；以年为单位，除以365
a=rt[,"expression"]<median(rt[,"expression"])
diff=survdiff(Surv(futime, fustat) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=round(pValue,5)
fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
summary(fit)    #查看五年生存率
#clusterNum=2#聚类分为几类
pdf(file="figures/survival-one.pdf")
plot(fit, lty = 2:3,
     # col=rainbow(clusterNum),
     col=c("red","blue"),
     xlab="time (day)",ylab="surival rate",
     main=paste("surival curve (p=", pValue ,")",
                sep=""))
legend("topright", c("BRCA1 high expression", "BRCA1 low expression"), lty = 2:3, col=c("red","blue"))
dev.off()

##### 多分组生存分析 #####
#基因名字
gene="rs121913529"

library(survival)
rt=read.delim("inst/survival/survivalmore.txt")
#如果以月为单位，除以30；以年为单位，除以365
rt$futime=rt$futime/365
diff=survdiff(Surv(futime, fustat) ~snp,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=round(pValue,4)
fit=survfit(Surv(futime, fustat) ~ snp, data = rt)
#查看五年生存率
summary(fit)

pdf(file="figures/survivalmore.pdf",width = 7,height =7)
plot(fit,
     lwd=2,
     col=c("red","blue","green","black"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste(gene," (p=", pValue ,")",sep=""))
legend("topright",
       c("CA","CC","CT","TT"),
       lwd=2,
       col=c("red","blue","green","black"))
dev.off()

##### 甲基化和表达量联合生存分析 #####
gene="ADCYAP1"#基因名字

library(survival)
rt=read.delim("inst/survival/survivalExpMethy.txt")
rt$futime=rt$futime/365                  #如果以月为单位，除以30；以年为单位，除以365
a=ifelse(rt[,"methy"]>=median(rt[,"methy"]),0,1)  #高甲基化0，低甲基化是1
b=ifelse(rt[,"exp"]<median(rt[,"exp"]),0,1)       #低表达是0，高表达是1
c=a+b
rt=cbind(rt,score=c)
rt=rt[(rt[,"score"]==0 | rt[,"score"]==2),]

diff=survdiff(Surv(futime, fustat) ~ score,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=round(pValue,3)
fit=survfit(Surv(futime, fustat) ~ score, data = rt)
summary(fit)                            #查看五年生存率

tiff(file="figures/survivalExpMethy.tiff",
     width = 15,
     height =15,
     units ="cm",
     compression="lzw",
     bg="white",
     res=300)
plot(fit,
     lty=2:3,
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     mark.time=T,
     main=paste("Survival curve (p=", pValue ,")",sep=""))
legend("topright",
       c(paste(gene," hyper & low expression",sep=""), paste(gene," hypo & high expression",sep="")),
       lty = 2:3,
       lwd=2,
       col=c("red","blue"))
dev.off()

##### 批量生存分析 #####
rm(list = ls())
library(survival)
pFilter=0.05
rt=read.delim("inst/survival/expTime.txt",
              check.names=F)

rt$futime=rt$futime/365                                               #如果以月为单位，除以30；以年为单位，除以365
outTab=data.frame()

for(gene in colnames(rt[,4:ncol(rt)])){
  if(median(rt[,gene])==0){
    next
  }else{a=factor(ifelse(rt[,gene]<=median(rt[,gene]),
                      "low","high"))}
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,
                            pvalue=pValue))

  fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
  summary(fit)

  if(pValue<pFilter){
    if(pValue<0.001){
      pValue=signif(pValue,4)
      pValue=format(pValue, scientific = TRUE)
    }else{
      pValue=round(pValue,3)
    }
    pdf(file=paste("figures/survival/",gene,".survival.pdf",sep=""),
        width = 5.5,            #图片的宽度
        height =5,              #图片的高度
    )
    plot(fit,
         lwd=2,
         col=c("red","blue"),
         xlab="Time (year)",
         mark.time=T,
         ylab="Survival rate",
         ylim=c(0,1.09),
         main=paste("Survival curve (p=", pValue ,")",
                    sep=""))
    legend("topright",
           c(paste(gene," high expression",sep=""),
             paste(gene," low expression",sep="") ),
           lwd=2,
           col=c("red","blue"))
    dev.off()
  }
}
write.table(outTab,file="results/survival.xls",
            sep="\t",row.names=F,quote=F)

##### survminer生存分析 #####
library(survival)
library(survminer)
rt=read.delim("inst/survival/ASrisk.txt")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

#绘制生存曲线
pdf(file="figures/survminer_survival.pdf",onefile = FALSE,
    width = 5.5,             #图片的宽度
    height =5)             #图片的高度
ggsurvplot(fit,
           data=rt,
           conf.int=TRUE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=TRUE,
           legend.labs=c("High risk", "Low risk"),
           legend.title="Risk",
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("red", "blue"),
           risk.table.height=.25)
dev.off()

summary(fit)    #查看五年生存率

##### 批量生存分析 #####
outTab=data.frame()

picDir="figures/survivalpicture"
dir.create(picDir)

library(survival)
library(qvalue)
rt=read.delim("inst/tumor.time.txt",
              row.names=1,check.names=F)
rt[,"futime"]=rt[,"futime"]/365

for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]

  med=median(rt[,i])
  if(med!=0){
    a=rt[,i]>med
    rt1=rt[a,]
    b=setdiff(rownames(rt),rownames(rt1))
    rt2=rt[b,]
    n1=nrow(rt1)
    n2=nrow(rt2)
    surTab1=summary(survfit(Surv(futime, fustat) ~ 1, data = rt1))
    surTab2=summary(survfit(Surv(futime, fustat) ~ 1, data = rt2))
    medianTab1=surTab1$table
    medianTab2=surTab2$table
    diff=survdiff(Surv(futime, fustat) ~a,data = rt)
    fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    outTab=rbind(outTab,cbind(gene=i,coxSummary$coefficients,coxSummary$conf.int,KM=pValue,
                              H_med=medianTab1["median"],H_0.95LCL=medianTab1["0.95LCL"],H_0.95UCL=medianTab1["0.95UCL"],
                              L_med=medianTab2["median"],L_0.95LCL=medianTab2["0.95LCL"],L_0.95UCL=medianTab2["0.95UCL"]))
    pval=0
    if(pValue<0.05){
      pval=signif(pValue,4)
      pval=format(pval, scientific = TRUE)
    }else{
      pval=round(pValue,3)
    }
    if(pValue<0.05){
      geneName=unlist(strsplit(i,"\\|",))[1]
      tiffFile=paste(geneName,".survival.tiff",sep="")
      outTiff=paste(picDir,tiffFile,sep="\\")
      tiff(file=outTiff,width = 15,height = 15,units ="cm",compression="lzw",bg="white",res=600)
      plot(fit, col=c("blue","red"), xlab="Time (years)", ylab="Overall survival",
           main=paste(geneName,"(p=",pval, ")", sep=""),mark.time=T,ylim=c(0,1.1),
           lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
      legend("topright", c(paste("Low expression"),
                           paste("High expression")),
             col=c("blue","red"), bty="n", lwd = 2, cex=0.8)
      dev.off()
    }
  }
}
```

## 相关性Cor

``` R
##### 基因表达量与甲基化的相关性 #####
#输入文件
inputFile="data-raw/methyGeneCor.txt"
#甲基化的基因
methyGene="ADCYAP1|methy"
#甲基化的位点
#methyGene="cg07211875|methy"
#表达的基因
expGene="ADCYAP1|exp"

rt <- read.delim(inputFile,
                 row.names = 1)
i=rt[methyGene,]
j=rt[expGene,]
x=as.numeric(i)
y=log2(as.numeric(j)+1)
corT=cor.test(x,y)
methyGeneName=unlist(strsplit(methyGene,"\\|",))[1]
expGeneName=unlist(strsplit(expGene,"\\|",))[1]
z=lm(y~x)
cor=corT$estimate
cor=round(cor,3)
pvalue=corT$p.value
pval=signif(pvalue,4)
pval=format(pval, scientific = TRUE)

#outTiff="cor.tiff"
#tiff(file=outTiff,width =15,height = 15,units ="cm",compression="lzw",bg="white",res=300)
plot(x,y,
     type="p",
     pch=16,
     main=paste("Cor=",cor," (p-value=",pval,")",sep=""),
     xlab=paste(methyGeneName,"methylation"),
     ylab=paste(expGeneName,"expression") )
lines(x,fitted(z),col=2)
#dev.off()
library(ggplot2)
p <- qplot(x, y,
           #data = NULL, color = clarity,
           xlab = paste(methyGeneName,"methylation"),
           ylab = paste(expGeneName,"expression"),
           main = paste("Cor=",cor," (p-value=",pval,")",sep=""))
p

##### 一个基因与其他基因的相关性检验 #####
library(limma)

# 输入文件
inputFile="data-raw/CorInput.txt"
#lncRNA名字
gene="LINC00526|lincRNA"
#相关系数过滤值
corFilter=0.5
#统计学p值过滤值
pFilter=0.001

picDir="CorPicture"
if(!dir.exists(picDir)){
  dir.create(picDir)
}

rt=read.delim(inputFile)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
exp=matrix(as.numeric(as.matrix(exp)),
            nrow=nrow(exp),dimnames=dimnames)
exp=avereps(exp)
group=sapply(strsplit(colnames(exp),"\\."),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
dat=exp[,group==0]
rt=dat[rowMeans(dat)>0.5,]

x=log2(as.numeric(rt[gene,])+1)
gene1=unlist(strsplit(gene,"\\|",))[1]
outputFile=paste("results/",gene1,
                 ".cor.xls",sep="")
outTab=data.frame()
for(j in rownames(rt)){
  y=log2(as.numeric(rt[j,])+1)
  gene2=unlist(strsplit(j,"\\|",))[1]
  gene2Type=unlist(strsplit(j,"\\|",))[2]
  if(gene2Type=="protein_coding"){
    corT=cor.test(x,y)
    gene1Name=unlist(strsplit(gene1,"\\|",))[1]
    gene2Name=unlist(strsplit(gene2,"\\|",))[1]

    z=lm(y~x)
    cor=corT$estimate
    cor=round(cor,3)
    pvalue=corT$p.value
    if(pvalue<0.001){
      pval=signif(pvalue,4)
      pval=format(pval, scientific = TRUE)
    }else{
      pval=round(pvalue,3)}

    #输出相关性图片
    if((abs(cor)>corFilter) & (pvalue<pFilter)){
      tiffFile=paste(gene1Name,"_",gene2Name,".cor.tiff",sep="")
      outTiff=paste(picDir,tiffFile,sep="\\")
      tiff(file=outTiff,width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
      plot(x,y, type="p",pch=16,col="blue",main=paste("Cor=",cor," (p-value=",pval,")",sep=""),
           cex=1, cex.lab=1, cex.main=1,cex.axis=1,
           xlab=paste(gene1Name,"expression"),
           ylab=paste(gene2Name,"expression") )
      lines(x,fitted(z),col=2)
      dev.off()
      outTab=rbind(outTab,
                   cbind(gene1,
                         gene2,
                         gene2Type,
                         cor,
                         pvalue))
      }

  }
}
save(outTab,
     file = "inst/Results/corTestResult.rda")
# write.table(file=outputFile,
#             outTab,sep="\t",
#             quote=F,row.names=F)

##### 可变剪接因子与事件相关性 #####
#相关系数过滤标准
corFilter=0.6
#p值过滤标准
pvalueFilter=0.001

#读取文件，并对样品取交集
#读取SF表达文件
SF = read.delim("data-raw/SFexp.txt", row.names=1)
#读取生存相关的AS文件
AS = read.delim("data-raw/uniSigExp.txt", row.names=1)
AS=t(AS[,3:ncol(AS)])
rownames(AS)=gsub("\\.","\\-",rownames(AS))
group=sapply(strsplit(colnames(SF),"\\."),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
SF=SF[,group==0]
colnames(SF)=gsub("(.*?)\\.(.*?)\\.(.*?)\\.(.*?)\\..*","\\1\\-\\2\\-\\3",colnames(SF))
sameSample=intersect(colnames(SF),
                     colnames(AS))
SF1=SF[,sameSample]
AS1=AS[,sameSample]

#相关性检验
outTab=data.frame()
for(i in row.names(SF1)){
  if(sd(SF1[i,])>1){
    for(j in row.names(AS1)){
      x=as.numeric(SF1[i,])
      y=as.numeric(AS1[j,])
      corT=cor.test(x,y)
      cor=corT$estimate
      pvalue=corT$p.value
      if((cor>corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(SF=i,AS=j,cor,pvalue,Regulation="postive"))
      }
      if((cor< -corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(SF=i,AS=j,cor,pvalue,Regulation="negative"))
      }
    }
  }
}
save(outTab,
     file = "results/AS-SF_Cor.Rdata")
#输出相关性结果
# write.table(file="corResult.txt",
#             outTab,sep="\t",
#             quote=F, row.names=F)

##### 免疫细胞与临床信息相关性 #####
inputFile="data-raw/immuneClinical.txt"
#输入文件
rt=read.delim(inputFile,
              row.names = 1)
#定义临床类型
clinical="grade"

picDir="CorPicture"
if(!dir.exists(picDir)){
  dir.create(picDir)
}

xlabel=vector()
tab1=table(rt[,clinical])
labelNum=length(tab1)
dotCol=c(2,3)
if(labelNum==3){
  dotCol=c(2,3,4)
}
if(labelNum==4){
  dotCol=c(2,3,4,5)
}
if(labelNum>4){
  dotCol=rainbow(labelNum)
}
for(i in 1:labelNum){
  xlabel=c(xlabel,names(tab1[i]) )
}

outTab=data.frame()

for(i in colnames(rt[,3:ncol(rt)])){
  rt1=data.frame(expression=rt[,i],
                 clinical=rt[,clinical])
  if(labelNum==2){
    wilcoxTest<-wilcox.test(expression ~ clinical, data=rt1)
  }else{
    wilcoxTest<-kruskal.test(expression ~ clinical, data = rt1)}
  pValue=wilcoxTest$p.value
  outTab=rbind(outTab,cbind(gene=i,pVal=pValue))
  pval=0
  if(pValue<0.001){
    pval=signif(pValue,4)
    pval=format(pval, scientific = TRUE)
  }else{
    pval=round(pValue,3)
  }

  b = boxplot(expression ~ clinical, data = rt1,outline = FALSE, plot=F)
  yMin=min(b$stats)
  yMax = max(b$stats/5+b$stats)
  ySeg = max(b$stats/10+b$stats)
  ySeg2 = max(b$stats/12+b$stats)
  n = ncol(b$stats)

  tiffFile=paste0(i,".",clinical,".tiff")
  outTiff=paste(picDir,tiffFile,sep="\\")
  tiff(file=outTiff,width = 25,height = 15,
       units ="cm",compression="lzw",bg="white",res=600)
  par(mar = c(4,7,3,3))
  boxplot(expression ~ clinical, data = rt1,names=xlabel,
          ylab = paste0(i," fraction"),col=dotCol,
          cex.main=1.6, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax),outline = FALSE)
  segments(1,ySeg, n,ySeg);
  segments(1,ySeg, 1,ySeg2)
  segments(n,ySeg, n,ySeg2)
  text((1+n)/2,ySeg,labels=paste0("p=",pval),cex=1.5,pos=3)
  dev.off()

}
save(outTab,
     file = "results/immuneCellCor.Rdata")

# write.table(outTab,
#             file=paste0("data", clinical,".xls"),
#             sep="\t",row.names=F,quote=F)

##### corrplot绘制相关性热图 #####
library(corrplot)
rt <- read.delim("data-raw/CorrplotInput.txt",
                 row.names = 1)
rt <- t(rt)
res1 <- cor.mtest(rt, conf.level = 0.95)
corrplot(corr=cor(rt),
         # method = "color",
         method = "circle",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         p.mat = res1$p,
         sig.level = 0.001,
         insig = "pch",
         number.cex = 1,
         type = "upper",
         col=colorRampPalette(c("blue", "white", "red"))(50),
)
M <- WGCNA::bicor(rt)
#type: "full", "upper" or "lower"
corrplot(M,type = "upper")
corrplot.mixed(M)
```

## 主成分分析PCA

``` R
# ==========================================================
#
#      PCA打分图
#      •   主成分分析
#      •   降维分析
#      •   用于区分离群样本
#
# ==========================================================

##### 主成分分析PCA #####

##### prcomp做主成分分析 #####
#读取表格
data=read.delim("data-raw/PCAdata.txt",
                row.names=1)
data=t(as.matrix(data))   #矩阵转置
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)   #PCA分析
# write.table(data.pca$rotation,file="results/PCA_PC.xls",
#             quote=F,sep="\t")   #输出特征向量
# write.table(predict(data.pca),file="results/PCA_newTab.xls",
#             quote=F,sep="\t")   #输出新表
pca.sum=summary(data.pca)
# write.table(pca.sum$importance,file="results/PCA_importance.xls",
#             quote=F,sep="\t")   #输出PC比重
# pdf(file="figures/PCABarplot.pdf",width=15)   #柱状图
barplot(pca.sum$importance[2,]*100,xlab="PC",ylab="percent",col="skyblue")
dev.off()
pdf(file="figures/PCAPlot.pdf",width=15)   #碎石图
plot(pca.sum$importance[2,]*100,type="o",col="red",xlab="PC",ylab="percent")
# dev.off()

### PCA 2d plot
### 免疫细胞信号值PCA聚类
library(ggplot2)
group=c(rep("con",5),rep("A",5),rep("B",3))  #需要大家输入
pcaPredict=predict(data.pca)
PCA = data.frame(PCA1 = pcaPredict[,1], PCA2 = pcaPredict[,2],group=group)
PCA.mean=aggregate(PCA[,1:2],list(group=PCA$group),mean)
#定义椭圆
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(as.factor(PCA$group))){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$group==g,],
                  veganCovEllipse(cov.wt(cbind(PCA1,PCA2),
                  wt=rep(1/length(PCA1),length(PCA1)))$cov,
                  center=c(mean(PCA1),mean(PCA2))))),group=g))
}

### 不画圈
ggplot(data = PCA, aes(PCA1, PCA2)) +
  geom_point(aes(color = group)) +
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# pdf(file="figures/PCA2d.pdf")
ggplot(data = PCA, aes(PCA1, PCA2)) + geom_point(aes(color = group)) +
  geom_path(data=df_ell, aes(x=PCA1, y=PCA2,colour=group), size=1, linetype=2)+
  annotate("text",x=PCA.mean$PCA1,y=PCA.mean$PCA2,label=PCA.mean$group)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# dev.off()

#pca 3d plot
library(pca3d)
library(rgl)
pca3d(data.pca, components = 1:3,
      group = group, show.centroids=TRUE,
      show.group.labels =TRUE)  #画3d图
#保存3d图形
rgl.snapshot("figures/PCA3d.png",fmt="png")

##### gmodel进行主成分分析 #####
library(gmodels)

## 输入
# input file name **
# inname = "data-raw/rpkm.txt"

# out PCA figure name **
# outname = "figures/all_DEGs.PCA.pdf"

# define the color for points  **
mycolors <- c(rep("red",2),rep("green",2),rep("yellow",2),rep("blue",2))

## step 1: 数据的读取和处理
# read the expr data
expr1 <- read.table(inname, header=T, row.names=1)
expr = t(scale(t(expr1)))


# transpose the data
data <- t(expr)


## step2：PCA分析
# do PCA
data.pca <- fast.prcomp(data)
#data.pca <- fast.prcomp(data,retx=T,scale=F,center=T)

## step3： PCA结果解析
# fetch the proportion of PC1 and PC2
# 一般情况下PC1 + PC2 > 70% 二维PCA散点图才有效
a <- summary(data.pca)
tmp <- a[4]$importance
pro1 <- as.numeric(sprintf("%.3f",tmp[2,1]))*100
pro2 <- as.numeric(sprintf("%.3f",tmp[2,2]))*100

# fetch the x axis min and max value (PC1)
xmax <- max(data.pca$x[,1])
xmin <- min(data.pca$x[,1])

# fetch the y axis min and max value (PC2)
ymax <- max(data.pca$x[,2])
ymin <- min(data.pca$x[,2])

# fetch sample names
samples =rownames(data.pca$x)

## step 4: 绘图
# draw PCA plot figure
# pdf(outname)

plot(
  data.pca$x[,1],
  data.pca$x[,2],
  xlab=paste("PC1","(",pro1,"%)",sep=""),
  ylab=paste("PC2","(",pro2,"%)",sep=""),
  main="PCA",
  xlim=c(1.1*xmin,1.1*xmax),
  ylim=c(1.1*ymin,1.1*ymax),
  pch=16,col=mycolors)

# 添加辅助线
abline(h=0,col="gray")
abline(v=0,col="gray")

# 添加样品名
text(data.pca$x[,1],data.pca$x[,2],labels=samples)

#dev.off()


##### FactoMineR和factoextra做主成分分析 #####

##### 各组织基因表达量的PCA分析
a <- read.table("data-raw/mcp.M113.035600-2.txt",
                head = TRUE, sep = "\t",
                check.names = FALSE,
                row.names = 1)
##第一列为基因ID，去掉第一列，
##把剩余的列保存到x里
x <- a[,-ncol(a)]
##转置
x <- t(x)
x <- log2(x+1)

# 加载PCA扩展包
library(FactoMineR)
library(factoextra)

# 进行PCA分析
dat.pca <- PCA(x, graph = FALSE)

# plot绘制PCA图
plot(dat.pca,choix="ind")

library(RColorBrewer)

fviz_pca_ind(dat.pca)

#ggsave('figs/PCA1.png')

##### airway数据示例
load("data/airway.rda")

data <- t(airway)
data <- as.data.frame(data)
group <- rep(c("control","treat"), c(4,4))
data <- cbind(data, group)

# pca分析
pca <- PCA(data[,-ncol(data)], graph = FALSE)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = data$group,
             addEllipses = TRUE,
             legend.title = "Groups"
)

pca <- PCA(data[,-ncol(data)], graph = FALSE)
fviz_pca_ind(pca,
             #axes=c(2,3),
             geom.ind = "point",
             col.ind = data$group,
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE,
             legend.title = "Groups"
)

##### iris数据集示例
#选iris来做pca分析示例，是因为他的分组间有差异，有研究意义
#去掉最后一列（分组信息）
dat <- iris[,-ncol(iris)]
table(iris$Species)
pdata=data.frame(Species=iris$Species)
rownames(pdata)=rownames(dat)

#观察这个图
pheatmap::pheatmap(dat)

# 选择要分析的主成分（一般是前两个）
#一个函数帮你完成复杂的pca统计分析
pca <- PCA(dat, graph = FALSE)
eig.val <- get_eigenvalue(pca)

###所有组成分的方差等信息
eig.val
###eigenvalue-方差；percentage of variance-方差占比；
###cumulative percentage of variance方差累积占比；
pca$eig
pca$eig[,2]
#特征值/方差/累计

#图1：碎石图
barplot(pca$eig[,2])
lines(x = 1:nrow(pca$eig),
      pca$eig[, 2],
      type="b", pch=19, col = "red")

fviz_eig(pca, addlabels = TRUE,
         ylim = c(0, 100))

#图2：样本聚类
fviz_pca_ind(pca, label="none",
             habillage=iris$Species,
             addEllipses=TRUE,
             ellipse.level=0.95,
             palette = "Dark2")
# Read more: http://www.sthda.com/english/wiki/ggplot2-colors

##图3：变量聚类
#####对样本进行作图，对应有对变量（当前场景，即基因水平）作图
#####通过设置axes可以进行其他主成分的展示
fviz_pca_var(pca, col.var = "contrib",
             gradient.cols = c("white", "blue", "red"),
             ggtheme = theme_minimal())
# 这个颜色是根据变量的贡献值赋值的

## 图4：结合变量和观测值，变量太多时不适用
fviz_pca_biplot(pca, label = "var", habillage=iris$Species,
                addEllipses=TRUE, ellipse.level=0.95,
                ggtheme = theme_minimal())
#只看对PC1
fviz_contrib(pca, choice = "var", axes = 1)
#只看对PC2
fviz_contrib(pca, choice = "var", axes = 2)
#综合看PC1+PC2
fviz_contrib(pca, choice = "var", axes = 1:2)
#变量数多时加参数：top = n，表示前n


```

## Diff_ChiTest

``` R
load(file = "data/ChiTestData.rda")
##### 卡方检验差异分析 #####
rt = ChiTestData[1:10,]

outTab=data.frame()
for(i in 1:nrow(rt)){
  x=matrix(c(rt[i,1],rt[i,2],rt[i,3],rt[i,4]), ncol = 2)
  chiTest=chisq.test(x)
  normalRatio=rt[i,1]/(rt[i,1]+rt[i,2])
  tumorRatio=rt[i,3]/(rt[i,3]+rt[i,4])
  Gene=row.names(rt[i,])
  Stat=chiTest$statistic
  Pvalue=chiTest$p.value
  outTab=rbind(outTab,
               cbind(Gene,
                     normalRatio,
                     tumorRatio,
                     Stat,
                     Pvalue))
}
pvalue=as.numeric(as.vector(outTab[,"Pvalue"]))
adjP=p.adjust(pvalue,method ="bonferroni")
outTab=cbind(outTab,adjPvalue=adjP)
#write.table(outTab,file="chiResult.txt",sep="\t",quote=F,row.names=F)
diffTab=outTab[outTab[,"adjPvalue"]<0.05,]
#write.table(diffTab,file="diff.txt",sep="\t",quote=F,row.names=F)

##### 卡方检验筛选临床信息1 #####
field="Cluster"
flag1="cluster1"
flag2="cluster2"

rt=clusterCliGroup
trainFlag=rt[rt[,field]==flag1,]
trainFlag=cbind(trainFlag,flag="Group1")
testFlag=rt[rt[,field]==flag2,]
testFlag=cbind(testFlag,flag="Group2")
newTable=rbind(trainFlag,testFlag)

newLabels=c("id")
for(i in 2:(ncol(rt)-1) ){
  nameStat=colnames(newTable)[i]
  tableStat=table(newTable[,c(nameStat,"flag")])
  pStat=chisq.test(tableStat)
  pvalue=pStat$p.value
  if(pvalue<0.001){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"***"))
  }else if(pvalue<0.01){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"**"))
  }else if(pvalue<0.05){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"*"))
  }else{
    newLabels=c(newLabels,colnames(newTable)[i])
  }
  print(paste(colnames(newTable)[i],pvalue,sep=" "))
}
newLabels=c(newLabels,colnames(newTable)[ncol(rt)])
colnames(rt)=newLabels
save()
# write.table(rt,file="results/clusterCliGroup.Sig.txt",
#             sep="\t",row.names=F,quote=F)


##### 卡方检验筛选临床信息2 #####
field="Risk"
flag1="low"
flag2="high"

rt=riskCliGroup
trainFlag=rt[rt[,field]==flag1,]
trainFlag=cbind(trainFlag,flag="Group1")
testFlag=rt[rt[,field]==flag2,]
testFlag=cbind(testFlag,flag="Group2")
newTable=rbind(trainFlag,testFlag)

newLabels=c("id")
for(i in 2:(ncol(rt)-1) ){
  nameStat=colnames(newTable)[i]
  tableStat=table(newTable[,c(nameStat,"flag")])
  pStat=chisq.test(tableStat)
  pvalue=pStat$p.value
  if(pvalue<0.001){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"***"))
  }else if(pvalue<0.01){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"**"))
  }else if(pvalue<0.05){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"*"))
  }else{
    newLabels=c(newLabels,colnames(newTable)[i])
  }
  print(paste(colnames(newTable)[i],pvalue,sep=" "))
}
newLabels=c(newLabels,colnames(newTable)[ncol(rt)])
colnames(rt)=newLabels
# write.table(rt,file="results/riskCliGroup.sig.txt",
#             sep="\t",row.names=F,quote=F)
```

## Diff_DESeq2

``` R
##### DESeq2数据差异分析 #####
library(DESeq2,quietly = TRUE)
library(limma,quietly = TRUE)
load(file = "data/DiffInputData.rda")

#读取mRNA表达矩阵数据
### 当行名有重复时
rt <- DiffInputData
# 有重复行名时，使用矩阵类型
#将数据框转化矩阵
rt=as.matrix(rt)
#将第一列基因做行名
rownames(rt)=rt[,1]
#第二列到最后一列作为基因表达数据
exp=rt[,2:ncol(rt)]
# 提取行名和列名
dimnames=list(rownames(exp),colnames(exp))
#将表达量转成数值型
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
# 对基因名相同的行取平均值
#对重复的基因取平均值
data=avereps(data)
# 提取基因表达量高于1的样本
#过滤表达量低的基因
data=data[rowMeans(data)>1,]

# 对数据取整
data=round(data,0)

#设置每个样本的实验条件,
#按照癌症和正常样品数目修改
group=c(rep("normal",4),
        rep("tumor",178))

# 将分组信息因子化
group_list = factor(group)
# 组合表达量与分组信息
colData <- data.frame(row.names=colnames(data),
                      group=group_list)
# 构建DESeq数据集用于后续分析
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = colData,
                              design = ~group)

###标准化 ###
dds2 <- DESeq(dds)
## 得到经过DESeq2软件normlization的表达矩阵！
# 50个样本以下的基因标准化
# rld <- rlogTransformation(dds2)
# data_new=assay(rld)
# 50个样本以上的基因标准化
vst <- vst(dds2)
data_new=assay(vst)
# 将标准化后的结果写入文件
# write.csv(data_new,file="gene_Norlization.csv",
#           quote=FALSE)

# 绘制直方图展示标准化前后结果
# par(cex = 0.7)
# n.sample=ncol(data)
# if(n.sample>40) par(cex = 0.5)
# cols <- rainbow(n.sample*1.2)
# par(mfrow=c(2,2))
# boxplot(data, col = cols,main="expression value",las=2)
# boxplot(data_new, col = cols,main="expression value",las=2)
# hist(data)
# hist(data_new)

###差异分析 ####
resultsNames(dds2)
# 设置要比较的分组
res <-  results(dds2,
                contrast=c("group","tumor","normal"),
                lfcThreshold = 1, alpha = 0.05,
                pAdjustMethod = "fdr")

## 提取想要的差异分析结果
resOrdered <- res[order(res$padj),]
# 转换为dataframe格式
resOrdered=as.data.frame(resOrdered)
# 去除na值
resOrdered=na.omit(resOrdered)
# 将所有差异分析结果写入文件
#write.csv(resOrdered,file="tumor_vs_normal.csv",quote = FALSE)
save(resOrdered, file = "data/DESeq2_dds2.rda")
#查看结果前六行
head(resOrdered)
```

## Diff_edgeR

``` R
##### edgeR数据差异分析 #####
library(edgeR, quietly = TRUE)
load(file = "data/DiffInputData.rda")

#读取mRNA表达矩阵数据
### 当行名有重复时
# 有重复行名时，使用矩阵类型
#将数据框转化矩阵
rt=as.matrix(DiffInputData)
#将第一列基因做行名
rownames(rt)=rt[,1]
#第二列到最后一列作为基因表达数据
exp=rt[,2:ncol(rt)]
#将表达量转成数值型
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
# 对基因名相同的行取平均值
#对重复的基因取平均值
data=avereps(data)
# 提取基因表达量高于1的样本
#过滤表达量低的基因
data=data[rowMeans(data)>1,]

#设置每个样本的实验条件,
#按照癌症和正常样品数目修改
group=c(rep("normal",4),rep("tumor",178))

# 差异分析
#构建设计矩阵
design <- model.matrix(~0+group)
#构建对比矩阵
y <- DGEList(counts=data,group=group)
#计算标准化因子
y <- calcNormFactors(y)
#估计散度
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
#计算差异表达
# fit <- glmFit(y, design)
# et <- glmLRT(fit, coef=2)
#两组数据间进行双尾检验
et <- exactTest(y, pair = c("normal","tumor"))
#对P值进行矫正并且提取出所有的基因
ordered_tags <- topTags(et, n=Inf,
                        adjust.method = "fdr",
                        sort.by = "PValue",
                        p.value = 0.05)
save(y, ordered_tags, file = "data/edgeR_y.rda")

# 提取差异基因
allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]

# edgeR标准化后的表达量
newData=y$pseudo.counts

Diffgene = allDiff[(allDiff$FDR < 0.05 &
                      (allDiff$logFC>2 |
                         allDiff$logFC<(-2))),]
#保存上调的差异基因
Upgene = allDiff[(allDiff$FDR < 0.05 &
                    (allDiff$logFC>2)),]

#保存下调的差异基因
Downgene = allDiff[(allDiff$FDR < 0.05 &
                      (allDiff$logFC<(-2))),]
```

## Diff_impute_limma

``` R
##### limma芯片数据差异分析 #####
library(limma)
#http://bioinf.wehi.edu.au/limma/data/Yoruba.RData
#http://eqtl.uchicago.edu/RNA_Seq_data/list_lanes_pickrell_2010_plosgenetics
#load("data/limma_data/Yoruba.RData")
#y

### 根据实际情况取log
# GeneExp <- log2(GeneExp+1)
GeneExp <- imputeData
### differential
### 样本数目需要修改
class <- c(rep("normal",10),rep("tumor",10))
design <- model.matrix(~0+factor(class))
colnames(design) <- c("normal","tumor")
rownames(class) <- colnames(GeneExp)
fit <- lmFit(GeneExp, design)
contrast.matrix<-makeContrasts(paste0(unique(class),
                                      collapse = "-"),
                               levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
# coef = "normal"
# adjust.method = "holm"
allDiff <- topTable(fit2, coef = 2, number = Inf,
                    adjust.method = 'fdr',
                    p.value = 0.05, lfc = 1)
save(fit2, allDiff, file = "data/limma_fit2.Rdata")
# write.table(allDiff, file="limmaTab.xls",
#             sep="\t", quote=F)



```

## Diff_WilcoxTest

``` R
##### Wilcox 检验做差异分析 #####

##### 批量做wilcox检验 #####
### 加载预处理包
library(limma)
load(file = "data/WilcoxTestData.rda")
#读取输入文件
rt <- as.matrix(WilcoxTestData)
rownames(rt) <- rt[,1]
exp <- rt[,2:ncol(rt)]
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data <- avereps(data)
data <- data[rowMeans(data)>0.5,]


# 根据中位数设置分组
score <- EstimateScores
Stromalmed <- median(score$StromalScore)
ImmuneScoremed <- median(score$ImmuneScore)
conTab <- score[score$StromalScore <
                  Stromalmed, ]
treatTab <- score[score$StromalScore >
                    Stromalmed, ]
con=as.vector(rownames(conTab))
treat=as.vector(rownames(treatTab))
conNum=length(con)
treatNum=length(treat)
#修改正常和癌症样品数目
# Type=c(rep("low",conNum),
#        rep("high",treatNum))
# Type=as.data.frame(Type)
# rownames(Type)=rownames(score)


### 设置分组
conNum <- 202 #normal组样品数目
treatNum <- 240 #tumor组样品数目
grade <- c(rep(1,conNum),rep(2,treatNum))

### 差异分析
fdrFilter <- 0.05 #fdr临界值
logFCfilter <- 2 #logFC临界值
outTab <- data.frame()


for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) |
      ((logFC<0) & (diffMed<0)) ){
    outTab=rbind(outTab,
                 cbind(gene=i,
                       conMean=conGeneMeans,
                       treatMean=treatGeneMeans,
                       logFC=logFC,
                       pValue=pvalue))
  }
}

# 根据p值标识显著性
pValue <- outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),
             method="fdr")
outTab <- cbind(outTab, fdr=fdr)

newGeneLists=c()

for (i in 1:nrow(outTab)) {
  gene <- outTab$gene[i]
  pvalue <- outTab$fdr[i]
  if(pvalue<0.001){
    newGeneLists=c(newGeneLists,paste0(gene,"***"))
  }else if(pvalue<0.01){
    newGeneLists=c(newGeneLists,paste0(gene,"**"))
  }else if(pvalue<0.05){
    newGeneLists=c(newGeneLists,paste0(gene,"*"))
  }else{
    newGeneLists=c(newGeneLists,gene)
  }
}

outTab <- cbind(newGeneLists, outTab)


#输出所有基因的差异情况
# write.table(outTab,file="data/all.xls",sep="\t",row.names=F,quote=F)

#输出差异表格
# outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
# write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)

#绘制热图需要的文件
# heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
# write.table(heatmap,file="diffGeneExp.txt",sep="\t",col.names=F,quote=F)


##### 单独做wilcox检验 #####
# 某基因与其snp
SNP="TP53|snp"
EXP="TP53|exp"
geneName="TP53"
load(file = "data/WilcoxTestEXPSNP.rda")
rt1=rt[c(EXP,SNP),]
rownames(rt1)=c("expression","snp")
rt=as.data.frame(t(rt1))

xlabel=vector()
tab1=table(rt[,"snp"])
labelNum=length(tab1)

for(i in 1:labelNum){
  xlabel=c(xlabel,paste(names(tab1[i]),"(n=",tab1[i],")",sep=""))
}

rt[,1]=as.numeric(as.vector(rt[,1]))

wilcoxTest<-wilcox.test(expression~snp,data=rt)
wilcoxP=wilcoxTest$p.value
pvalue=signif(wilcoxP,4)
pval=0
if(pvalue<0.001){
  pval=signif(pvalue,4)
  pval=format(pval, scientific = TRUE)
}else{
  pval=round(pvalue,3)
}

b = boxplot(expression ~ snp, data = rt,outline = FALSE, plot=F)
yMin=min(b$stats)
yMax = max(b$stats/5+b$stats)
ySeg = max(b$stats/10+b$stats)
ySeg2 = max(b$stats/12+b$stats)
n = ncol(b$stats)

#pdfFile=paste(geneName,".pdf",sep="")
#pdf(file=pdfFile,width=8,height=7)
par(mar = c(4,7,3,3))
boxplot(expression ~ snp, data = rt,names=xlabel,
        ylab = paste(geneName," expression",sep=""),col=c("green","red"),
        cex.main=1.6, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax),outline = FALSE)
segments(1,ySeg, n,ySeg);
segments(1,ySeg, 1,ySeg2)
segments(n,ySeg, n,ySeg2)
text((1+n)/2,ySeg,labels=paste("p=",pval,sep=""),cex=1.5,pos=3)
#dev.off()
```

## kruskaltest分析

``` R
##### kruskaltest分析 #####
load(file = "data/kruskaltestData.rda")

rt=kruskaltestData

#删除正常样品
group=sapply(strsplit(colnames(rt),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
rt=rt[,group==0]

#统计分析，输出文件
outTab=data.frame()
for(i in 1:(nrow(rt)/2)){
  cnv=rownames(rt[i,])
  exp=gsub("cnv","exp",cnv)
  geneName=gsub("\\|cnv","",cnv)
  data=rbind(cnv=rt[cnv,],exp=log2(rt[exp,]+1))
  data=t(data)
  ksTest<-kruskal.test(exp ~ cnv, data=data)
  ksPval=ksTest$p.value
  outTab=rbind(outTab,cbind(geneName,pvalue=ksPval))
}
# write.table(outTab,file="results/kruskaltest_ks.xls",
#             sep="\t",row.names=F,quote=F)

##### 多分组kruskaltest分析 #####
SNP="rs121913529|KRAS"
EXP="KRAS|exp"
geneName="KRAS"
snpName="rs121913529"

rt1=kruskalMergeData[c(EXP,SNP),]
rownames(rt1)=c("expression","snp")
rt=as.data.frame(t(rt1))

xlabel=vector()
tab1=table(rt[,"snp"])
labelNum=length(tab1)
for(i in 1:labelNum){
  xlabel=c(xlabel,paste(names(tab1[i]),"(n=",tab1[i],")",sep=""))
}

rt[,1]=as.numeric(as.vector(rt[,1]))
wilcoxTest=kruskal.test(expression~snp,rt)
wilcoxP=wilcoxTest$p.value
pvalue=signif(wilcoxP,4)
pval=0
if(pvalue<0.001){
  pval=signif(pvalue,4)
  pval=format(pval, scientific = TRUE)
}else{
  pval=round(pvalue,3)
}

b = boxplot(expression ~ snp, data = rt,outline = FALSE, plot=F)
yMin=min(b$stats)
yMax = max(b$stats/5+b$stats)
ySeg = max(b$stats/10+b$stats)
ySeg2 = max(b$stats/12+b$stats)
n = ncol(b$stats)

pdfFile=paste("figures/kruskal_",geneName,".pdf",sep="")
pdf(file=pdfFile,width=10,height=7)
par(mar = c(4.5,7,3,3))
boxplot(expression ~ snp, data = rt,names=xlabel,
        ylab = paste(geneName," expression",sep=""),col=c("green","red","blue","black"),xlab=snpName,
        cex.main=1.6, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax),outline = FALSE)
segments(1,ySeg, n,ySeg);
segments(1,ySeg, 1,ySeg2)
segments(n,ySeg, n,ySeg2)
text((1+n)/2,ySeg,labels=paste("p=",pval,sep=""),cex=1.5,pos=3)
dev.off()
```

## lasso_cox

``` R
##### 按照剪接类型区分 #####
asType="ME"
#(AA AD AP AT ES RI ME)中的一种,
#如果7种一起做，asType设置为空即可
rt=read.delim("inst/survival/uniSigExp.txt",
              row.names=1,check.names=F)#读取文件
rt$futime[rt$futime<=0]=1
rt$futime=rt$futime/365
genes=colnames(rt)
gene=grep(paste0("\\|",asType),genes,value=T)
geneLength=ifelse(length(gene)>20,20,length(gene))
rt=rt[,c("futime","fustat",gene[1:geneLength])]

##### lasso回归筛选模型 #####
library(glmnet)
library(survival)

rt=read.delim("data-raw/lassoInput.txt",row.names=1)#读取文件

### cox回归筛选基因
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
# pdf("inst/Figures/lasso_lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
# dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
# pdf("figures/lasso_cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
# dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
# write.table(lassoGene,file="inst/Results/lassoGene.txt",
#             sep="\t",quote=F,row.names=F,col.names=F)

riskScore=predict(cvfit, newx = x, s = "lambda.min",type="response")
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(riskScore),risk)

# write.table(cbind(id=rownames(outTab),outTab),
#             file="inst/Results/lassoRisk.txt",
#             sep="\t",
#             quote=F,
#             row.names=F)

##### riskplot图绘制1 #####
library(pheatmap)
rt=read.delim("inst/survival/ASrisk.txt",
              row.names=1,check.names=F)       #读取输入文件
rt=rt[order(rt$riskScore),]                                     #按照riskScore对样品排序

#绘制风险曲线
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
# pdf(file="inst/Figures/riskScore.pdf",width = 12,height = 5)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
abline(h=median(rt$riskScore),
       v=lowLength,lty=2)
#dev.off()

#绘制生存状态图
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="green"
#pdf(file="inst/Figures/risksurvStat.pdf",width = 12,height = 5)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
#dev.off()

#绘制风险热图
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
# pdf(file="inst/Figures/riskheatmap.pdf",width = 12,height = 5)
pheatmap(rt1,
         annotation=annotation,
         cluster_cols = FALSE,
         fontsize_row=11,
         fontsize_col=3,
         color = colorRampPalette(c("green", "black", "red"))(50) )
#dev.off()

##### riskplot图绘制2 #####
rt=read.delim("inst/survival/ASrisk.txt",
              row.names=1,check.names=F)       #读取输入文件
rt=rt[order(rt$riskScore),]
ggplot(rt,aes(x=risk,y=riskScore,
              fill=risk))+
        geom_boxplot()+
        xlab("Risk")+ylab("Risk Score")+
        labs(fill="Risk")


ggplot(rt,aes(x=risk, y=riskScore,
               fill=risk))+
        geom_dotplot(binaxis = "y",
                     binwidth = 0.2,
                     stackdir = "center")+
        xlab("Risk")+ylab("Risk Score")+
        labs(fill="Risk")


attach(rt)
plot(riskScore,futime,
     col=ifelse(fustat==1,"red","green"),
     xlab = "",ylab="Survival time",
     pch=16)
legend("topright", c("Death", "Alive"),
       pch=16, col=c("red","green"))



ggplot(rt,aes(x=riskScore,
                y=futime))+
        geom_area(aes(fill=risk))+
        geom_line()+
        geom_hline(yintercept = 0)+
        xlab("Risk Score")+
        ylab("Survival time")+
        labs(fill="Risk")
```

## ROC曲线

``` R
##### ROC曲线 #####
library(survivalROC)

load(file = "data/SurvivalData.rda")
rt <- ASrisk
#pdf(file="figures/ROC.pdf")
par(oma=c(0.5,1,0,1),
    font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore,
                predict.time =3, method="KM")
plot(roc$FP, roc$TP,
     type="l", xlim=c(0,1),
     ylim=c(0,1),col='red',
     xlab="False positive rate",
     ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",
                round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3,
     cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
#dev.off()
#AUC
# 0.5      #没有任何预测能力
# 0.51-0.7 #低准确度
# 0.71-0.9 #中等准确度
# >0.9     #高准确度


##### 多条ROC曲线 #####
library(survival)
library(timeROC)

predict_3_year<- 3
predict_5_year<- 5

ROC<-timeROC(T=rt$futime,delta=rt$fustat,
             marker=rt$riskScore,cause=1,
             weighting="marginal",
             times=c(predict_3_year,
                     predict_5_year),
             ROC=TRUE)

#pdf("figures/ROC2.pdf")
plot(ROC,time=predict_3_year,title=FALSE,lwd=3)
plot(ROC,time=predict_5_year,col="blue",add=TRUE,title=FALSE,lwd=3)
legend("bottomright",
       c(paste("AUC of 3 year survival: ",round(ROC$AUC[1],3)),
         paste("AUC of 5 year survival: ",round(ROC$AUC[2],3))),col=c("red","blue"),lwd=3)
#dev.off()
```

## 
