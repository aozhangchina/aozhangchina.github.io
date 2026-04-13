if (!exists('inputing_method')) {inputing_method <- "possibility"} #? ML/REML/possibility未添加，预留
if (!exists('method')) {method <- "percentage"}   #? "percentage"/"leave_one"/"k-fold"
if (!exists('trainingSize')) {trainingSize <- 0.5}   # training population size
if (!exists('cycles')) {cycles <- 100}   # 设置循环次数
if (!exists('frequency_NA')) {frequency_NA <- 0.15}   # 允许最大缺失值
if (!exists('polymorphism')) {polymorphism <- TRUE}   # 多态性检测
if (!exists('frequency_MAF')) {frequency_MAF <- 0.05}   # 最小等位基因需要默认大于0.05
if (!exists('output_filter_HMP')) {output_filter_HMP <- FALSE}   # 输出筛选后的基因型
if (!exists('Geno_Probability_imputation')) {Geno_Probability_imputation <- FALSE}   # 杂交种合成前用概率方法补缺失
if (!exists('Synthetic_ignore_gene')) {Synthetic_ignore_gene <- TRUE}   # 合成忽略NA
if (!exists('output_numeric_geno')) {output_numeric_geno <- TRUE}   # 输出数值化基因型
if (!exists('Geno_before_imputation')) {Geno_before_imputation <- "ML"}   # 分析前的补缺失"ML"/"REML"/"Probability"
if (!exists('package_choose')) {package_choose <- "rrBLUP"}   #  选择方法：目前只支持rrBLUP
if (!exists('predict_type')) {predict_type <- "hybrid/SCA"}   # "hybrid/SCA" or "Line/GCA"
if (!exists('predict_inpop')) {predict_inpop <- TRUE}   # 是否为群体内预测TRUE or FALSE
if (!exists('k_fold')) {k_fold <- 5}   # # 所有个体平均分成k份，每次预测其中的一份
if (!exists('k_fold_Phenogroup')) {k_fold_Phenogroup <- TRUE}   ### 将材料按照表型分成（材料数/k组），向上取整

  
# rrBLUP杂交种及combining ability预测
# 编写：张敖
# 更新：2020-08-14 16:39:13
# https://datahold.cn


#↓ 函数部分开始 ↓#
#! 读取文件
readFiles <- function(header=TRUE,choose=FALSE,fname){
    if(!"readr" %in% installed.packages()) {
        install.packages("readr")
        library(readr)
    }else{
        library(readr)
    }
    if (choose==TRUE){
        myFileName <- file.choose()
    }else{
        myFileName <- fname
    }
    if (grepl("\\.csv$",myFileName)){
        if (header==TRUE){
            myFile <- read_csv(myFileName, col_names = TRUE)
        }else{
            myFile <- read_csv(myFileName, col_names = FALSE)
        }
    }else{
        if (header==TRUE){
            myFile <- read_tsv(myFileName, col_names = TRUE)
        }else{
            myFile <- read_tsv(myFileName, col_names = FALSE)
        }
    }
    return(myFile)
}

#! 功能函数：变为双碱基形式
douobleBase <- function(x){
  tema <- x
  if (!length(which(is.na(tema)==T))==length(tema)){
    tema <- sub("A","AA",tema)
    tema <- sub("T","TT",tema)
    tema <- sub("C","CC",tema)
    tema <- sub("G","GG",tema)
    tema <- sub("R","AG",tema)
    tema <- sub("Y","CT",tema)
    tema <- sub("S","CG",tema)
    tema <- sub("W","AT",tema)
    tema <- sub("K","GT",tema)
    tema <- sub("M","AC",tema)
    b <- tema
    b[which(is.na(tema))] <- ""   # NA替换成""
    b <- paste(b,collapse = "")   # 合并字母
    b <- unlist(strsplit(b,split=""))   # 拆分字母
    cb <- as.data.frame(table(b))   # 获取字母个数
    cb <- cb[order(cb[,2],decreasing = TRUE),]   # 倒序排列
    A <- as.character(cb[1,1])
    a <- as.character(cb[2,1])
    tema <- gsub(A,1,tema)
    tema <- gsub(a,0,tema)
    if (nrow(cb)>2){
      tema <- gsub("[^1|^0]",NA,tema)
    }
    tema <- gsub(1,"A",tema)
    tema <- gsub(0,"a",tema)
    for(xv in 1:length(tema)){
      if (!is.na(tema[xv])){
        temb <- unlist(strsplit(tema[xv],split=""))
        temb <- sort(temb,decreasing=T)
        tema[xv] <- paste(temb,collapse = "")
      }else{
        tema[xv] <- NA
      }
    }
  }
  return(tema)
}

#! 功能函数： 最大等位基因2，最小等位基因0，杂合1，缺失值NA
genoImputeM <- function(x) {
  len <- c( length(x[x=="A"]), length(x[x=="T"]), length(x[x=="C"]),length(x[x=="G"]))   # ATCG个数
  len1 <- order(len,na.last=F, decreasing = T)   # 降序排列
  refer <- c("A","T","C","G")
  MajorA <- NULL
  MinorA <- NULL
  MajorA <- refer[len1[1]]
  MinorA <- refer[len1[2]]
  heto <- c("R","Y","S","W","K","M")
  index1 <- which(!x %in% c(MajorA,MinorA,heto))
  x[index1] <- NA
  x <- sub(MajorA,2,x)
  x <- sub(MinorA,0,x)
  index2 <- which(x %in% heto)
  x[index2]= 1
  return(as.numeric(x))
}  

#! 安装和加载包
libraryAndInstall <- function(X_name){
  tt <- paste0("library(",X_name,")")
  if(!X_name %in% installed.packages()) {
    install.packages(X_name)
    eval(parse(text=tt))
  }else{
    eval(parse(text=tt))
  }
}

#! 删除无多态性的标记
delet_monomorphic <- function(myVector){
  monomorphic <- numeric()
  if(length(as.numeric(table(myVector)))<=1)
  {
    monomorphic <- 1   # 1表示无多态性
  }else{
    monomorphic <- 0  # 0表示有多态性
  }
  return(monomorphic)
}
#! 在向量中利用可能性补缺失
probability_impute <- function(myVector){
  if (any(is.na(myVector)==TRUE)){
    table_vector <- table(myVector)
    tempA <- names(table_vector)
    impute_value <- sample(x=tempA,size=sum(is.na(myVector)),replace=T,prob=table_vector/length(myVector))
    myVector[which(is.na(myVector))] <- impute_value
  }
  return(myVector)
}
#! 结合无多态性标记和可能性补缺失两个步骤的自动化函数
impute_prob <- function(myMatrix){
  cat("Status:正在去除无多态性标记!\n")
  system.time(index_monomorphic <- apply(geno,2,delet_monomorphic))  # 获得无多态性标记索引
  X_without_monomorphic <- geno[,index_monomorphic==0]
  print_info <- data.frame("Individual"=dim(X_without_monomorphic)[1],"Markers"=dim(X_without_monomorphic)[2])
  print(print_info)
  numbersOfNA <- data.frame("NA"=length(which(is.na(X_without_monomorphic))))
  print(numbersOfNA)
  cat("Status:正在执行Imputation!\n")
  impu_prob_value <- apply(X_without_monomorphic,2,probability_impute)
  cat("Status:Imputation complete!\n")
  return(impu_prob_value)
}
#! 结合无多态性标记和似然补缺失两个步骤的自动化函数
impute_rrBLUP <- function(myMatrix,inputing_method){
  cat("Status:正在去除无多态性标记!\n")
  system.time(index_monomorphic <- apply(geno,2,delet_monomorphic))  # 获得无多态性标记索引
  X_without_monomorphic <- geno[,index_monomorphic==0]
  print_info <- data.frame("Individual"=dim(X_without_monomorphic)[1],"Markers"=dim(X_without_monomorphic)[2])
  print(print_info)
  numbersOfNA <- data.frame("NA"=length(which(is.na(X_without_monomorphic))))
  print(numbersOfNA)

  return(X_without_monomorphic)
}
#! 功能函数，字母排序联用函数1
Aasort1 <- function(Aav){
  Aav <- as.matrix(Aav)
  sAa <- apply(Aav,1,Aasort2)
  return(sAa)
}
#! 功能函数，字母排序联用函数2
Aasort2 <- function (Aa) {
  if (!is.na(Aa)){
    Aax <- unlist(strsplit(Aa,split=""))
    Aax <- sort(Aax,decreasing=T)
    Aax <- paste(Aax,collapse = "")
  }else{
    Aax <- NA
  }
  return(Aax)
}
#! 功能函数：杂交种基因型合成
Synthetic_hybrids <- function(x){
  if (length(grep("@@",x))>0){
    xx <- strsplit(x,split="@@")
    if (Synthetic_ignore_gene==TRUE){
      x1 <- GD_D[,xx[[1]][1]]
      x2 <- GD_D[,xx[[1]][2]]
      if (any(is.na(x1)==TRUE)){
        x1[which(is.na(x1)==TRUE)] <- ""
      }
      if (any(is.na(x2)==TRUE)){
        x2[which(is.na(x2)==TRUE)] <- ""
      }
      xresu <- paste0(x1,x2)
    }else{
      xresu <- paste0(GD_D[,xx[[1]][1]],GD_D[,xx[[1]][2]])
      xresu <- gsub("NA",NA,xresu)
    }
  }else{
    xx <- strsplit(x,split="@")
    if (Synthetic_ignore_gene==TRUE){
      x1 <- GD_D[,xx[[1]][1]]
      x2 <- GD_D[,xx[[1]][2]]
      if (any(is.na(x1)==TRUE)){
        x1[which(is.na(x1)==TRUE)] <- ""
      }
      if (any(is.na(x2)==TRUE)){
        x2[which(is.na(x2)==TRUE)] <- ""
      }
      xresu <- paste0(x1,x2)
    }else{
      xresu <- paste0(GD_D[,xx[[1]][1]],GD_D[,xx[[1]][2]])
      xresu <- gsub("NA",NA,xresu)
    }
  }
  return(xresu)
}
#! 功能函数：杂交种基因型数值化联用函数1
numericH <- function(x){
  x <- as.matrix(x)
  nHv <- apply(x,1,numericH2)
  return(nHv)
}

#! 功能函数：杂交种基因型数值化联用函数1(单个值，每个A代表的值)
numericH2 <- function(Hx){
  if (!is.na(Hx)){
    xmax <- nchar(max(Hx,na.rm=T))
    eachA <- 1/xmax
    eachCell <- unlist(strsplit(Hx,split="")) 
    countA <- length(which(eachCell=="A"))
    Hvalue <- countA*eachA
  }else{
    Hvalue <- Hx
  }
  Hvalue <- as.numeric(Hvalue)
  return(Hvalue)
}
#↑ 函数部分结束 ↑#

#↓ 分析部分开始 ↓#
## 加载rrBLUP
libraryAndInstall("rrBLUP")
## 加载dplyr
libraryAndInstall("dplyr")

if (!exists("myGenoName")){
  cat("请选择基因型文件...\n Please choose a genotipic file in HMP format...\n")
  myGenoName <- choose.files()
}else{
  if (length(myGenoName)==0){
    cat("请选择基因型文件...\n Please choose a genotipic file in HMP format...\n")
    myGenoName <- choose.files()
  }
}
if (!exists("myPhenoName")){
    cat("请选择表型文件...\n Please choose a phenotipic file in TXT format...\n")
    myPhenoName <- choose.files()
}else{
  if (length(myPhenoName)==0){
    cat("请选择表型文件...\n Please choose a phenotipic file in TXT format...\n")
    myPhenoName <- choose.files()
  }
}


myGeno <- readFiles(header=FALSE,fname=myGenoName)   # Import the Genotipic file.

## 检查基因型数据是否正确
bit=nchar(as.character(myGeno[2,12])) 
if(!bit==1) {stop("Data format is hapmap with haploid allele! Please Check Genotype data!")}

myPheno <- readFiles(header=TRUE,fname=myPhenoName)   # Import the Phenotipic file.
if (length(which(is.na(myPheno$TESTER)==TRUE))==nrow(myPheno)) {
  noTester <- TRUE
}else{
  noTester <- FALSE
  # 判断TESTER是否为杂交种,并分开杂交种
  if (length(grep("@",myPheno$TESTER))>0){
    Tester0 <- strsplit(myPheno$TESTER[grep("@",myPheno$TESTER)],split="@")
    Tester1 <- NULL
    Tester2 <- NULL
    for (li in 1:length(Tester0)){
      Tester1 <- c(Tester1,Tester0[[li]][1])
      Tester2 <- c(Tester2,Tester0[[li]][2])
    }
    Tester1 <- unique(Tester1)
    Tester2 <- unique(Tester2)
    Testerx <- myPheno$TESTER[-c(grep("@",myPheno$TESTER))]
    TESTERS <- unique(c(Tester1,Tester2,Testerx))
  }else{
    TESTERS <- unique(c(myPheno$TESTER))
  }
}
# 判断Line是否为杂交种,并分开杂交种
if (length(grep("@",myPheno$LINE))>0){
  Line0 <- strsplit(myPheno$LINE[grep("@",myPheno$LINE)],split="@")
  Line1 <- NULL
  Line2 <- NULL
  for (li in 1:length(Line0)){
    Line1 <- c(Line1,Line0[[li]][1])
    Line2 <- c(Line2,Line0[[li]][2])
  }
  Line1 <- unique(Line1)
  Line2 <- unique(Line2)
  Linex <- myPheno$LINE[-c(grep("@",myPheno$LINE))]
  LINES <- unique(c(Line1,Line2,Linex))
}else{
  LINES <- unique(c(myPheno$LINE))
}


## 获取性状名称
index <- names(myPheno) %in% c("LOC","ENV","LINE","TESTER","YEAR","P1P2")
fstName <- names(myPheno)[!index]   # The trait name.

## 自动设置输出目录
    newDirectory <- paste0(dirname(myPhenoName),"/","GenoInfo")
    if(!dir.exists(newDirectory)){
        dir.create(newDirectory)
    }
    setwd(newDirectory)
cat("Output Directory:",getwd(),"\n")

## 基因型和表型对应##################################################
myGenoInformation <- myGeno[,1:11]   # 基因型数据基本信息
myGenoMarkers <- myGeno[,12:ncol(myGeno)]   # 基因型数据标记信息
myGenoIndividuals <- as.character(myGenoMarkers[1,])   # 提取所有的材料名
myGenoMarkers <- myGenoMarkers[,order(myGenoIndividuals)]   # 按照材料名排序标记信息部分

## 基因型与表型材料名称对应
if (noTester==T){
  myPhenoN <- LINES
}else{
  myPhenoN <- unique(c(LINES,TESTERS))
}
myGenoN <- unique(myGenoIndividuals)
bothGP <- intersect(myPhenoN,myGenoN)
### 有表型的基因型数据
index <- myGenoIndividuals %in% bothGP
cat("Unmatched genotype: \n")
print(unique(myGenoIndividuals[!index]))
myGenoMarkers <- as.matrix(myGenoMarkers[,index])
myGenoIndividuals <- as.character(myGenoMarkers[1,])   # 提取所有的材料名
myGenoMarkers <- myGenoMarkers[,order(myGenoIndividuals)]   # 按照材料名排序标记信息部分
myGeno <- cbind(myGenoInformation,myGenoMarkers)   # 重新合成基因型数据（排序后）
#myGenoIndividualsT <- t(myGeno[1,-(1:11)])   # 获得排序后的材料名（包括Lines和Testers（亲本1或亲本2））
### 有基因型的表型数据
#### LINE
if (length(grep("@",myPheno$LINE))){
  Line00 <- strsplit(myPheno$LINE,split="@")
  index <- myPheno$LINE %in% bothGP   # 查看Line中基因型和表型的对应索引(非杂交种的匹配情况)
  for (li in grep("@",myPheno$LINE)){
    index[li] <- all((Line00[[li]] %in% bothGP)==TRUE)   # 修正杂交种的匹配情况
  }
  cat("Unmatched Line: \n")
  print(unique(myPheno$LINE[!index]))   # 查看是否有未匹配的项并列出
}else{
  index <- myPheno$LINE %in% bothGP   # 查看Line中基因型和表型的对应索引
  cat("Unmatched Line: \n")
  print(unique(myPheno$LINE[!index]))   # 查看是否有未匹配的项并列出
}
myPheno <- myPheno[index,]   # 去掉未对应表型数据
#### TESTER
if (noTester==F){
  if (length(grep("@",myPheno$TESTER))){
    Tester00 <- strsplit(myPheno$TESTER,split="@")
    index <- myPheno$TESTER %in% bothGP   # 表型Tester中基因型和表型的对应索引
    for (li in grep("@",myPheno$TESTER)){
      index[li] <- all((Tester00[[li]] %in% bothGP)==TRUE)   # 修正杂交种的匹配情况
    }
    cat("Unmatched Line: \n")
    print(unique(myPheno$LINE[!index]))   # 查看是否有未匹配的项并列出
  }else{
    index <- myPheno$TESTER %in% bothGP   # 表型Tester中基因型和表型的对应索引
    cat("Unmatched Tester: \n")
    print(unique(myPheno$TESTER[!index]))   # 查看是否有未匹配的项并列出
  }
  myPheno <- myPheno[index,]   # 去掉未对应的表型数据
}
####################################################

if (myGeno[1,1]=="rs#"){
  names(myGeno) <- myGeno[1,]   # 为数据框添加标题
  myGeno <- myGeno[-1,]   # 去掉第一行
}

GD <- as.matrix(myGeno[,-(1:11)])   # get genotype matrix   获得基因型矩阵
## "N"换成NA
if (any(GD=="N")){
  GD <- sub("N",NA,GD)
}
myGeno[,12:ncol(myGeno)] <- GD   # 将myGeno替换成有NA形式
Marker <- myGeno[,c(1,3,4)]   # 获得标记信息数据框
colnames(Marker) <- c("SNP","Chromosome","Position")   # 修改标记数据框的列名，对应相应数据
rownames(GD) <- Marker$SNP   # 增加标记名称
No.markers <- character()
No.markers["original"] <- nrow(GD)

## 多态性检测
if (polymorphism == TRUE){
  cat("Status: Deleting no polymorphism loci!\n")
  system.time(index_monomorphic <- apply(GD,1,delet_monomorphic))  # 获得无多态性标记索引
  myGeno <- myGeno[index_monomorphic==0,]
  GD <- GD[index_monomorphic==0,]
  Marker <- Marker[index_monomorphic==0,]   # 获得标记信息数据框
  rm(index_monomorphic)
}

## 单碱基变双碱基
No.markers["polymorphism"] <- nrow(GD)   # 获得筛选后的标记数量
cat("Getting genotype of doube-base type: \n")
system.time(GD_D <- apply(GD,1,douobleBase))   # 单碱基变双碱基
GD_D <- t(as.matrix(GD_D))   # 转置，所有分析基于该矩阵


## NA频率检测
freqNA <- rowSums(is.na(GD_D))/ncol(GD_D)    # NA频率
pdf("freqNA.pdf")
hist(freqNA,xlim=c(0,1.0),xlab="Missing value frequency")   # NA直方图
dev.off()
freqNAInfo <- as.character(summary(freqNA))
names(freqNAInfo) <- names(summary(freqNA))
freqNAInfo["No.markers"] <- No.markers[length(No.markers)]
write.table(freqNAInfo,"freqNA.txt",sep="\t",col.names=F,quote=F)  # NA基本信息

## 缺失值筛选
index <- freqNA < frequency_NA   #? 设置缺失值筛选条件
GD <- GD[index,]   # 筛选缺失值
myGeno <- myGeno[index,]   # 筛选缺失值
Marker <- Marker[index,]   # 筛选缺失值
GD_D <- GD_D[index,]   # 筛选缺失值
freqNA <- rowSums(is.na(GD_D))/ncol(GD_D)    # NA频率
pdf("freqNA2.pdf")
hist(freqNA,xlim=c(0,1.0),xlab="Missing value frequency")   # NA直方图
dev.off()
freqNAInfo <- as.character(summary(freqNA))
names(freqNAInfo) <- names(summary(freqNA))
No.markers["filterNA"] <- nrow(GD)
freqNAInfo["No.markers"] <- No.markers["filterNA"]
write.table(freqNAInfo,"freqNA2.txt",sep="\t",col.names=F,quote=F)  # NA基本信息

## MAF筛选
GD[which(is.na(GD_D)==T)] <- NA   # 将超过2种碱基变异转换的NA赋值给GD
## 数字转换0,1,2
GD2 <- NULL   # 初始化GD2变量，用于存放数值型基因型数据
system.time({
GD2 <- apply(GD, 1, genoImputeM)   # 数字转换
})
GD2 <- t(matrix(as.numeric(as.matrix(GD2)),nrow=nrow(GD2)))   # 转为数值型并转置
rownames(GD2) <- Marker[,1]   # 增加标记名称
colnames(GD2) <- names(myGeno[-1:-11])   # 增加材料名称(包括所有TESTER和LINE)
## 计算MAF,等位基因频率需要0,1,2编码
phat <- rowMeans(GD2,na.rm=T)/2   # 计算每列的均均值
maf <- ifelse(phat<0.5,phat,1-phat)   # 计算最小等位基因频
pdf("MAF.pdf")
hist(maf,xlab="Minor allele frequency",xlim=c(0,0.5))   # 绘制最小等位基因频率图
dev.off()
MAFInfo <- as.character(summary(maf))   # MAF基本信息
names(MAFInfo) <- names(summary(maf))
MAFInfo["No.markers"] <- No.markers[length(No.markers)]
write.table(MAFInfo,"maf.txt",sep="\t",col.names=F,quote=F)   # NA基本信息
## 根据MAF筛选
index <- maf > frequency_MAF   #? 获得MAF筛选索引
GD2 <- GD2[index,]   # 按MAF筛选基因型数据
GD <- GD[index,]   # 按MAF筛选基因型数据
myGeno <- myGeno[index,]   # 按MAF筛选基因型数据
Marker <- Marker[index,]   # 按MAF筛选基因型数据
GD_D <- GD_D[index,]   # 按MAF筛选基因型数据
maf <- maf[index]   # 按MAF筛选maf
phat <- rowMeans(GD2,na.rm=T)/2   # 计算每列的均均值
maf <- ifelse(phat<0.5,phat,1-phat)   # 计算最小等位基因频
pdf("MAF2.pdf")
hist(maf,xlab="Minor allele frequency",xlim=c(0,0.5))   # 绘制最小等位基因频率图
dev.off()
MAFInfo <- as.character(summary(maf))   # MAF基本信息
names(MAFInfo) <- names(summary(maf))
No.markers["filterMAF"] <- nrow(GD2)
MAFInfo["No.markers"] <- No.markers["filterMAF"]
write.table(MAFInfo,"maf2.txt",sep="\t",col.names=F,quote=F)   # NA基本信息

### 输出hmp格式筛选后的基因型
if (output_filter_HMP==TRUE){
  myGenoNewName <- paste0(gsub(".txt|.hmp","",basename(myGenoName)),"_Available.hmp.txt")
  myGeno2 <- as.matrix(myGeno)
  GD3 <- myGeno2[,12:ncol(myGeno2)]
  GD3[which(is.na(GD3)==T)] <- "N"
  myGeno2[,12:ncol(myGeno2)] <- GD3
  write.table(myGeno2,myGenoNewName,sep="\t",quote=F,row.names=F)
  rm(myGeno2,GD3)
}

### 概率方案补缺失
if (Geno_Probability_imputation==TRUE){
  cat("Status: Imputating!\n")
  system.time(impu_prob_value <- apply(GD_D,1,probability_impute))
  GD_D <- t(impu_prob_value)
}


### 杂交种合成
## Line的杂交种合成
if (length(grep("@",myPheno$LINE))>0){
  HLine <- as.matrix(myPheno$LINE[grep("@",myPheno$LINE)])
  system.time(HLGeno <- apply(HLine,1,Synthetic_hybrids))  # 获得无多态性标记索引
  rownames(HLGeno) <- rownames(GD_D)
  colnames(HLGeno) <- HLine[,1]
  system.time(test <- apply(HLGeno,1,Aasort1))
  HLGeno <- t(test)
  if (nrow(HLGeno)==1){
    HLGeno <- t(HLGeno)
    colnames(HLGeno) <- HLine[,1]
  }
  if (!HLine[1,1] %in% colnames(GD_D)){
    GD_D <- cbind(GD_D,HLGeno)
  }
}
## Tester的杂交种合成
if (length(grep("@",myPheno$TESTER))>0){
  HTester <- as.matrix(myPheno$TESTER[grep("@",myPheno$TESTER)])
  system.time(HTGeno <- apply(HTester,1,Synthetic_hybrids))  # 获得无多态性标记索引
  rownames(HTGeno) <- rownames(GD_D)
  colnames(HTGeno) <- HTester[,1]
  system.time(test <- apply(HTGeno,1,Aasort1))
  HTGeno <- t(test)
  if (nrow(HTGeno)==1){
    HTGeno <- t(HTGeno)   # 只有一个时转换
    colnames(HTGeno) <- HTester[,1]
  }
  if (!HTester[1,1] %in% colnames(GD_D)){
    GD_D <- cbind(GD_D,HTGeno)
  }
}
### Line×Tester合成
if (noTester==F){
  HH <- as.matrix(paste0(myPheno$LINE,"@@",myPheno$TESTER))
  myPheno$LINExTESTER <- HH[,1]
  system.time(HHGeno <- apply(HH,1,Synthetic_hybrids))
  rownames(HHGeno) <- rownames(GD_D)
  colnames(HHGeno) <- HH[,1]
  system.time(test <- apply(HHGeno,1,Aasort1))
  HHGeno <- t(test)
  if (nrow(HHGeno)==1){
    HHGeno <- t(HHGeno)
    colnames(HHGeno) <- HLine[,1]
  }
  if (!HHGeno[1,1] %in% colnames(GD_D)){
    GD_D <- cbind(GD_D,HHGeno)
  }
}


### 杂交种基因型数值化
GD_DN <- apply(GD_D,2,numericH)

### 输出数值化基因型
if (output_numeric_geno==TRUE){
  myGenoNewName <- paste0(gsub(".txt|.hmp","",basename(myGenoName)),"_numeric_Geno.txt")
  GD_DN2 <- as.data.frame(GD_DN)
  
  if (!names(GD_DN2)[1]=="Markers"){
    GD_DN2 <- cbind(row.names(GD_DN2),GD_DN2)
    names(GD_DN2)[1] <- "Markers"
  }
  write.table(GD_DN2,myGenoNewName,sep="\t",quote=F,row.names=F)
  rm(GD_DN2)
}

### 输出标记筛选变化
write.table(No.markers,"No.markers.txt",quote=F,col.names=F)

# 获取基因型名称
finalGenoMaterial <- colnames(GD_DN)

####? 将基因型矩阵转换为rrBLUP/BGLR适用的格式
if (package_choose=="rrBLUP"){
  if (!min(GD_DN,na.rm=T)==-1){
    GD_DN <- GD_DN*2-1
  }
  if (colnames(GD_DN)[1] %in% finalGenoMaterial){
    GD_DN <- t(GD_DN)
  }
}

##### 最后一次补缺失
if (Geno_before_imputation =="ML"){
  impute <- A.mat(GD_DN,impute.method="ML",return.imputed=T)   #### not remove missing
  str(impute)
  Markers_impute <- impute$imputed
}else if (Geno_before_imputation =="REML"){
  impute <- A.mat(GD_DN,impute.method="REML",return.imputed=T)   #### not remove missing
  str(impute)
  Markers_impute <- impute$imputed
}else{
  system.time(impu_prob_value <- apply(GD_DN,2,probability_impute))
  Markers_impute <- apply(impu_prob_value,2,as.numeric)
  rownames(Markers_impute) <- finalGenoMaterial
}
## 检查是否还有NA
if(any(is.na(Markers_impute))) {stop("")}else{cat("No NA! Pass!\n")}



## 检查基因型和表型的对应情况
### 根据预测类型选择表型是Line还是Hybrid
if (predict_type=="hybrid/SCA"){
  Phenos <- unique(HH[,1])   # 获取表型名称
}else{
  Phenos <- unique(myPheno$LINE)   # 获取表型名称
}

resultNames <- intersect(finalGenoMaterial,Phenos)   # 基因型表型取交集
length(finalGenoMaterial)   # 基因型
length(Phenos)   # 表型
length(resultNames)   # 交集

### 
if (predict_inpop==TRUE){
  pop <- 1
}else{
  pop <- NA
}
assistCol <- cbind(resultNames,pop)   # 获得辅助列
colnames(assistCol) <- c("names","pop")   # 增加列名
assistCol <- assistCol[order(as.numeric(assistCol[,2])),]   # 按照数字顺序排序
if (!predict_inpop==TRUE){
  write.csv(assistCol,"finalTaxa.csv",row.names = F)   # 生成文件，手动编号
  cat("Please assign training population and test population.(1/2)\n1 = Training population, \n2 = Test population.\nModify file: ",getwd(),"/finalTaxa.csv\nThen choose this file again in the pop-up window.\n",sep="")
  assistCol <- read.csv(choose.files(),header=T)   # 读取编号后的文件
}

#### 基因型对应
index　<-  rownames(Markers_impute) %in% assistCol[,"names"]
Markers_impute <- Markers_impute[index,]   # 按材料名筛选
#### 表型对应
if (predict_type=="hybrid/SCA"){
  index <- myPheno$LINExTESTER %in% assistCol[,"names"]
}else{
  index <- myPheno$LINE %in% assistCol[,"names"]
}
myPheno <- myPheno[index,]   # 按材料名筛选
##### 表型基因型排序
if (predict_type=="hybrid/SCA"){
  myPheno <- myPheno[order(myPheno$LINExTESTER),]
}else{
  myPheno <- myPheno[order(myPheno$LINE),]
}
Markers_impute <- Markers_impute[order(rownames(Markers_impute)),]
# if(!all(Phenos==finalGenoMaterial)) {stop("There are mismatched material names!")}


#? 功能函数，获得rrBLUP的预测参数
calculateRRBLUP <- function(x){
  ans <- mixed.solve(x, Z=m_train)   # 计算
  e <- as.matrix(ans$u)   # 获得BLUP的随机效应
  pred_temp_valid <-  m_valid %*% e   # 与验证群体进行矩阵乘法运算
  vv <- list("ans"=ans,"e"=e,"pred_temp_valid"=pred_temp_valid)   # 估计育种值转表型值
  return(vv)
}

# 所有性状循环
for (i in 1:length(fstName)){
  cat ("processing",fstName[i],"......\n")
  if(!dir.exists(fstName[i])){
    dir.create(fstName[i])
  }
  ### 群体内验证 ###
  ####↓ 留一法 method="leave_one"
  if (method=="leave_one" & predict_inpop==TRUE){
    ## 初始化建模和预测等变量
    Pheno_train <- NULL   # 建模群体的表型
    m_train <- NULL   # 建模群体的基因型
    Pheno_valid <- NULL   # 验证群体的表型
    m_valid <- NULL   # 验证群体的基因型
    r_MG <- NULL   # 预测精度
    e <- NULL   # 初始化e
    pred_valid_r <- NULL   # 初始化最终预测值
    pred_test_value_r <- NULL   # 初始化最终预测表型值
    
    ### 材料个数计算
    nS <- length(resultNames)    #  Line或杂交种的材料总个数
    
    
    
    for(r in 1:nS){
      test <- r   # 唯一用作验证的材料
      train <- setdiff(1:nS,test)   # 剩余部分建模
      test_index <- rep(FALSE,nS)
      test_index[r] <- TRUE   # 用于筛选验证群体的子集
      train_index <- test_index==FALSE   # 用于筛选建模群体的子集
      Pheno_valid <- as.matrix(subset(myPheno[,fstName[i]],test_index))   # 验证群体的表型
      m_valid <- as.matrix(subset(Markers_impute, test_index))   # 验证群体的基因型
      
      Pheno_train <- as.matrix(subset(myPheno[,fstName[i]],train_index))  # 建模群体的表型
      m_train <- as.matrix(subset(Markers_impute, train_index))   # 建模群体的基因型
      
      PredictRRBLUP <- calculateRRBLUP(Pheno_train)
      pred_valid_r[r] <- PredictRRBLUP$pred_temp_valid   # 添加每一个估计育种值
    }
    pred_test_value_r <- pred_valid_r + as.numeric(PredictRRBLUP$ans$beta)   # 育种值转换为表型值
    r_MG <- as.matrix(cor(pred_valid_r, as.numeric(unlist(myPheno[,fstName[i]])), use="complete" ))   # 预测精度
    colnames(r_MG)<- paste(fstName[i],"Accuracy",sep=".")   # 为预测精度添加列名
    rownames(r_MG)<- "Accuracy"   # 为预测精度添加行名
    if (predict_type=="hybrid/SCA"){
      names(pred_valid_r) <- myPheno$LINExTESTER
      names(pred_test_value_r) <- myPheno$LINExTESTER   # 为预测表型值添加材料名
    }else{
      names(pred_valid_r) <- myPheno$LINE
      names(pred_test_value_r) <- myPheno$LINE   # 为预测表型值添加材料名
    }
    ## 输出预测精度
    write.csv(r_MG,paste(fstName[i],"/Accuracy-",method,"-",fstName[i],".csv",sep=''),row.names=T)
    ## 输出预测的表型值
    write.table(pred_test_value_r,paste0(fstName[i],"/",fstName[i],"-",method,"-Predictive pheno value",".csv"),sep=",",col.names=F)  
    ## 输出预测育种值
    write.table(pred_valid_r,paste0(fstName[i],"/",fstName[i],"-",method,"-Predictive breeding value",".csv"),sep=",",col.names=F) 
  }
  ####↑ 留一法 method="leave_one"
  ####↓ 百分比预测 method="percentage"
  if (method=="percentage" & predict_inpop==TRUE){
    ## 初始化建模和预测等变量
    Pheno_train <- NULL   # 建模群体的表型
    m_train <- NULL   # 建模群体的基因型
    Pheno_valid <- NULL   # 验证群体的表型
    m_valid <- NULL   # 验证群体的基因型
    r_MG <- matrix()   # 预测精度
    e <- NULL   # 初始化e
    pred_valid_r <- NULL   # 初始化最终预测值
    pred_test_value_r <- NULL   # 初始化最终预测表型值
    
    ### 材料个数计算
    nS <- length(resultNames)    #  Line或杂交种的材料总个数
    if (predict_type=="hybrid/SCA"){
      all_material <-  as.data.frame(myPheno$LINExTESTER)   # 设置材料名用于获取每一次的预测结果
      all_beta <- as.data.frame(myPheno$LINExTESTER)   # 设置材料名用于获取每一次的预测表型结果
    }else{
      all_material <-  as.data.frame(myPheno$LINE)   # 设置材料名用于获取每一次的预测结果
      all_beta <- as.data.frame(myPheno$LINE)   # 设置材料名用于获取每一次的预测表型结果
    }
    colnames(all_material) <- "Material"
    colnames(all_beta) <- "Material"
    for (r in 1:cycles) {
      train <- sample(1:nS, ceiling(nS*trainingSize))   # 将小数形式的百分比换算成具体个数,向上取整并抽取用于建模
      test <- setdiff(1:nS,train)   # 区域部分用来验证
      train_index <- rep(FALSE,nS)
      train_index[train] <- TRUE   # 用于筛选建模群体的子集
      test_index <- train_index==FALSE   # 用于筛选验证群体的子集
      
      Pheno_valid <- as.matrix(subset(myPheno[,fstName[i]],test_index))   # 验证群体的表型
      m_valid <- as.matrix(subset(Markers_impute, test_index))   # 验证群体的基因型
      
      Pheno_train <- as.matrix(subset(myPheno[,fstName[i]],train_index))  # 建模群体的表型
      m_train <- as.matrix(subset(Markers_impute, train_index))   # 建模群体的基因型
      
      PredictRRBLUP <- calculateRRBLUP(Pheno_train)
      cycName <- paste0("Cyc",r)   # 循环列名，避免程序警告
      txt <- paste0("pred_temp_valid2 <- data.frame('Material'=row.names(PredictRRBLUP$pred_temp_valid),'",cycName,"'=PredictRRBLUP$pred_temp_valid)")
      eval(parse(text = txt))   # 生成数据框用于合并
      all_material <- merge(all_material,pred_temp_valid2,by='Material',all.x=TRUE)
      txt <- paste0("all_beta$",cycName,"<-as.numeric(unlist(all_material$",cycName,"))+as.numeric(PredictRRBLUP$ans$beta)") 
      eval(parse(text = txt))   # 生成数据框用于合并
      r_MG[r] <- cor(PredictRRBLUP$pred_temp_valid, Pheno_valid, use="complete" )   # 预测精度
      names(r_MG)[r] <- cycName
    }
    pred_valid_r <- rowMeans(all_material[,-1],na.rm=T)   # 估计育种值均值
    all_material$mean <- pred_valid_r
    pred_test_value_r <- rowMeans(all_beta[,-1],na.rm=T)   # 估计育种值转表型值
    all_beta$mean <- pred_test_value_r
    
    r_MG <- as.matrix(r_MG)
    colnames(r_MG)<- paste(fstName[1],"Accuracy",sep=".")   # 为预测精度添加列名
    if (predict_type=="hybrid/SCA"){
      names(pred_valid_r) <- myPheno$LINExTESTER
      names(pred_test_value_r) <- myPheno$LINExTESTER   # 为预测表型值添加材料名
    }else{
      names(pred_valid_r) <- myPheno$LINE
      names(pred_test_value_r) <- myPheno$LINE   # 为预测表型值添加材料名
    }
    
    ## 输出预测精度
    write.csv(r_MG,paste(fstName[i],"/Accuracy-",method,"-",paste0(trainingSize*100,"%"),"-",fstName[i],".csv",sep=''),row.names=T)
    ## 输出平均预测的表型值
    write.table(all_beta,paste0(fstName[i],"/",fstName[i],"-",method,"-",paste0(trainingSize*100,"%"),"-Predictive pheno value",".csv"),sep=",",col.names=T,row.names = F)  
    ## 输出平均预测育种值
    write.table(all_material,paste0(fstName[i],"/",fstName[i],"-",method,"-",paste0(trainingSize*100,"%"),"-Predictive breeding value",".csv"),sep=",",col.names=T,row.names = F)
  }
  ####↑ 百分比预测 method="percentage"
  ####↓ k倍交叉验证 method="k-fold"
  if (method=="k-fold" & predict_inpop==TRUE){
    ## 初始化建模和预测等变量
    Pheno_train <- NULL   # 建模群体的表型
    m_train <- NULL   # 建模群体的基因型
    Pheno_valid <- NULL   # 验证群体的表型
    m_valid <- NULL   # 验证群体的基因型
    r_MG <- matrix()   # 预测精度
    e <- NULL   # 初始化e
    pred_valid_r <- NULL   # 初始化最终预测值
    pred_test_value_r <- NULL   # 初始化最终预测表型值
    all_material <- data.frame()
    ### 材料个数计算
    nS <- length(resultNames)    #  Line或杂交种的材料总个数
    if (predict_type=="hybrid/SCA"){
      all_material <-  as.data.frame(myPheno$LINExTESTER)   # 设置材料名用于获取每一次的预测结果(k倍)
      all_beta <- as.data.frame(myPheno$LINExTESTER)   # 设置材料名用于获取每一次的预测表型结果(k倍)
      all_material2 <-  as.data.frame(myPheno$LINExTESTER)   # 设置材料名用于获取每一次的预测结果（循环）
      all_beta2 <- as.data.frame(myPheno$LINExTESTER)   # 设置材料名用于获取每一次的预测表型结果（循环）
    }else{
      all_material <-  as.data.frame(myPheno$LINE)   # 设置材料名用于获取每一次的预测结果(k倍)
      all_beta <- as.data.frame(myPheno$LINE)   # 设置材料名用于获取每一次的预测表型结果(k倍)
      all_material2 <-  as.data.frame(myPheno$LINE)   # 设置材料名用于获取每一次的预测结果（循环）
      all_beta2 <- as.data.frame(myPheno$LINE)   # 设置材料名用于获取每一次的预测表型结果（循环）
    }
    colnames(all_material) <- "Material"
    colnames(all_beta) <- "Material"
    for (r in 1:cycles) {
      if (k_fold_Phenogroup==TRUE){
        PhenoOrder <- order(myPheno[,fstName[i]])   # 表型排序
        phenoGroup <- ceiling(nS/k_fold)   # 分成多少个组
        tempDF <- data.frame(name=1:nS,phenoGroup=NA,ran=NA)   # 初始化分组数据框
        tempDF <- tempDF[PhenoOrder,]
        if (nS%%(phenoGroup)!=0){
          tempDF$phenoGroup <- sort(c(rep(1:phenoGroup,nS/phenoGroup),1:(nS%%(phenoGroup))))   # 分组编号
        }else{
          tempDF$phenoGroup <- sort(rep(1:phenoGroup,nS/phenoGroup))
        }

        for (phi in 1:phenoGroup){
          tempDF$ran[which(tempDF$phenoGroup==phi)] <- sample(1:length(which(tempDF$phenoGroup==phi)))   # 分组抽取并填充数据框
        }
        tempDF <- tempDF[order(tempDF$name),]
        k_time <- tempDF$ran   # 按照
      }else{
        if (nS%%(k_fold)!=0){
          k_time <- c(rep(1:k_fold,nS/k_fold),1:(nS%%(k_fold)))   # 获得分组编号
        }else{
          k_time <- rep(1:k_fold,nS/k_fold)
        }
        k_time_random <- sample(1:nS,nS)   # 获得随机数
        k_time <- k_time[k_time_random]
      }

      for (k in 1:k_fold){
        test <- which(k_time==k)   # 根据分组，选择验证群体
        train <- setdiff(1:nS,test)   # 区域部分用来验证
        train_index <- rep(FALSE,nS)
        train_index[train] <- TRUE   # 用于筛选建模群体的子集
        test_index <- train_index==FALSE   # 用于筛选验证群体的子集
      
        Pheno_valid <- as.matrix(subset(myPheno[,fstName[i]],test_index))   # 验证群体的表型
        m_valid <- as.matrix(subset(Markers_impute, test_index))   # 验证群体的基因型
      
        Pheno_train <- as.matrix(subset(myPheno[,fstName[i]],train_index))  # 建模群体的表型
        m_train <- as.matrix(subset(Markers_impute, train_index))   # 建模群体的基因型
      
        PredictRRBLUP <- calculateRRBLUP(Pheno_train)
        cycName <- paste0("Cyc",r)   # 循环列名，避免程序警告
        
        if (!cycName %in% names(all_material)){
          txt <- paste0("all_material <- cbind(all_material,",cycName,"=NA)")
          eval(parse(text = txt))   # 填写数据
        }
        
        for (kk in 1:length(test)){
          txt <- paste0("all_material$",cycName,"[test[",kk,"]]<- as.numeric(PredictRRBLUP$pred_temp_valid[",kk,"])")
          eval(parse(text = txt))   # 填写数据
        }
        txt <- paste0("all_beta$",cycName,"<-as.numeric(unlist(all_material$",cycName,"))+as.numeric(PredictRRBLUP$ans$beta)") 
        eval(parse(text = txt))   # 生成数据框用于合并
      }
      r_MG[r] <- cor(all_material[,(1+r)],as.numeric(unlist(myPheno[,fstName[i]])),use="complete")
      names(r_MG)[r] <- cycName
    }
    pred_valid_r <- rowMeans(all_material[,-1],na.rm=T)   # 估计育种值均值
    all_material$mean <- pred_valid_r
    pred_test_value_r <- rowMeans(all_beta[,-1],na.rm=T)   # 估计育种值转表型值
    all_beta$mean <- pred_test_value_r
    r_MG <- as.matrix(r_MG)
    colnames(r_MG)<- paste(fstName[1],"Accuracy",sep=".")   # 为预测精度添加列名
    if ("Mean" %in% rownames(r_MG)==FALSE){
      r_MG_r <- mean(r_MG[,1])
      r_MG <- rbind(r_MG,"Mean"=r_MG_r)
    }
    if (predict_type=="hybrid/SCA"){
      names(pred_valid_r) <- myPheno$LINExTESTER
      names(pred_test_value_r) <- myPheno$LINExTESTER   # 为预测表型值添加材料名
    }else{
      names(pred_valid_r) <- myPheno$LINE
      names(pred_test_value_r) <- myPheno$LINE   # 为预测表型值添加材料名
    }
    
    ## 输出预测精度
    write.csv(r_MG,paste(fstName[i],"/Accuracy-",method,"-",k_fold,"-",fstName[i],".csv",sep=''),row.names=T)
    ## 输出平均预测的表型值
    write.table(all_beta,paste0(fstName[i],"/",fstName[i],"-",method,"-",k_fold,"-Predictive pheno value",".csv"),sep=",",col.names=T,row.names = F)  
    ## 输出平均预测育种值
    write.table(all_material,paste0(fstName[i],"/",fstName[i],"-",method,"-",k_fold,"-Predictive breeding value",".csv"),sep=",",col.names=T,row.names = F)
  }
  ####↓ k倍交叉验证 method="k-fold"
  ### 群体间验证 ###
  if (predict_inpop==FALSE){
    ## 初始化建模和预测等变量
    Pheno_train <- NULL   # 建模群体的表型
    m_train <- NULL   # 建模群体的基因型
    Pheno_valid <- NULL   # 验证群体的表型
    m_valid <- NULL   # 验证群体的基因型
    r_MG <- matrix()   # 预测精度
    e <- NULL   # 初始化e
    pred_valid_r <- NULL   # 初始化最终预测值
    pred_test_value_r <- NULL   # 初始化最终预测表型值
    pred_temp_valid2 <- NULL   # 初始化最终结果数据框
  
    
    ### 材料个数计算
    nS <- length(resultNames)    #  Line或杂交种的材料总个数


      train <- which(assistCol$pop==1)   # 获取建模群体
      test <- which(assistCol$pop==1)   # 获取测试群体
      train_index <- rep(FALSE,nS)
      train_index[train] <- TRUE   # 用于筛选建模群体的子集
      test_index <- rep(FALSE,nS)
      test_index[test] <- TRUE   # 用于筛选验证群体的子集
      
      Pheno_valid <- as.matrix(subset(myPheno[,fstName[i]],test_index))   # 验证群体的表型
      m_valid <- as.matrix(subset(Markers_impute, test_index))   # 验证群体的基因型
      
      Pheno_train <- as.matrix(subset(myPheno[,fstName[i]],train_index))  # 建模群体的表型
      m_train <- as.matrix(subset(Markers_impute, train_index))   # 建模群体的基因型
      
      PredictRRBLUP <- calculateRRBLUP(Pheno_train)
      
      # 
      pred_temp_valid2 <- data.frame('Material'=row.names(PredictRRBLUP$pred_temp_valid),'GEBV'=PredictRRBLUP$pred_temp_valid)   # GEBV
      pred_temp_valid2$PPV <- as.numeric(unlist(pred_temp_valid2$GEBV))+as.numeric(PredictRRBLUP$ans$beta)   # predicted phenotypic value
      r_MG[i] <- cor(PredictRRBLUP$pred_temp_valid, Pheno_valid, use="complete" )   # 预测精度

    r_MG <- as.matrix(r_MG)
    colnames(r_MG)<- paste(fstName[1],"Accuracy",sep=".")   # 为预测精度添加列名
    
    ## 输出预测精度
    write.csv(r_MG,paste(fstName[i],"/Accuracy-",method,"-",paste0("between_pop"),"-",fstName[i],".csv",sep=''),row.names=T)
    ## 输出平均预测的表型值
    write.table(pred_temp_valid2[c("Material","PPV")],paste0(fstName[i],"/",fstName[i],"-",method,"-",paste0("between_pop"),"-Predictive pheno value",".csv"),sep=",",col.names=T,row.names = F)  
    ## 输出平均预测育种值
    write.table(pred_temp_valid2[c("Material","GEBV")],paste0(fstName[i],"/",fstName[i],"-",method,"-",paste0("between_pop"),"-Predictive breeding value",".csv"),sep=",",col.names=T,row.names = F)
  }
}

#↑ 分析部分结束 ↑#
