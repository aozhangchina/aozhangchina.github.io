rm(list=ls())
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

## 加载rrBLUP
libraryAndInstall("rrBLUP")
## 加载dplyr
libraryAndInstall("dplyr")

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

polymorphism <- TRUE

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
## 多态性检测
if (polymorphism == TRUE){
  cat("Status: Deleting no polymorphism loci!\n")
  system.time(index_monomorphic <- apply(GD,1,delet_monomorphic))  # 获得无多态性标记索引
  myGeno <- myGeno[index_monomorphic==0,]
  GD <- GD[index_monomorphic==0,]
  Marker <- Marker[index_monomorphic==0,]   # 获得标记信息数据框
  rm(index_monomorphic)
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
## 单碱基变双碱基
No.markers["polymorphism"] <- nrow(GD)   # 获得筛选后的标记数量
cat("Getting genotype of doube-base type: \n")
system.time(GD_D <- apply(GD,1,douobleBase))   # 单碱基变双碱基
GD_D <- t(as.matrix(GD_D))   # 转置，所有分析基于该矩阵

Geno_Probability_imputation <- TRUE
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
### 概率方案补缺失
if (Geno_Probability_imputation==TRUE){
  cat("Status: Imputating!\n")
  system.time(impu_prob_value <- apply(GD_D,1,probability_impute))
  GD_D <- t(impu_prob_value)
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
### 杂交种合成
Synthetic_ignore_gene <- TRUE
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
### 杂交种基因型数值化
GD_DN <- apply(GD_D,2,numericH)

### 输出数值化基因型
output_numeric_geno <- TRUE
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
