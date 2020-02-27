# Read files
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

# get file information
getFileInfo <- function(){
myFileName <- choose.files()
f.name <- basename(myFileName)
f.dir <- dirname(myFileName)
f <- list(f.name=f.name,f.dir=f.dir)
return(f)
}

#! 功能函数： 缺失值转换为NA
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
  x <- sub(MajorA,1,x)
  x <- sub(MinorA,0,x)
  index2 <- which(x %in% heto)
  x[index2]= 0.5
  return(as.numeric(x))
}
#  删除无多态性的标记
delet_monomorphic <- function(myVector){
  monomorphic <- numeric()
  if(length(as.numeric(table(myVector)))==1)
  {
    monomorphic <- 1   # 1表示无多态性
  }else{
    monomorphic <- 0  # 0表示有多态性
  }
  return(monomorphic)
}
# 在向量中利用可能性补缺失
probability_impute <- function(myVector){
  table_vector <- table(myVector)
  tempA <- as.numeric(names(table_vector))
  table_vector/sum(myVector)
  impute_value <- sample(x=tempA,size=sum(is.na(myVector)),replace=T,prob=table_vector/sum(table_vector))
  myVector[which(is.na(myVector))] <- impute_value
  return(myVector)
}
# 结合上面两个步骤的自动化函数
impute_prob <- function(myMatrix){
  cat("Status: Removing polymorphism markers.\n")
  system.time(index_monomorphic <- apply(geno,2,delet_monomorphic))  # 获得无多态性标记索引
  X_without_monomorphic <- geno[,index_monomorphic==0]
  print_info <- data.frame("Individual"=dim(X_without_monomorphic)[1],"Markers"=dim(X_without_monomorphic)[2])
  print(print_info)
  numbersOfNA <- data.frame("NA"=length(which(is.na(X_without_monomorphic))))
  print(numbersOfNA)
  cat("Status: Imputating!\n")
  impu_prob_value <- apply(X_without_monomorphic,2,probability_impute)
  cat("Status: Imputation complete!\n")
  return(impu_prob_value)
}

# code start
cat("Please choose a file of hmp.\n")
fileInfo <- getFileInfo()   # get file info

cat("Please choose a group file.\n")
fgroup <- choose.files()
setwd(fileInfo$f.dir[1])   # set workspace
cat(paste0("Work Space:\n",getwd(),"\n"))   # echo workspace

# read file
geno <- readFiles(header=FALSE,fname=fileInfo$f.name[1])

bit=nchar(as.character(geno[2,12])) 
if(!bit==1) {stop("Data format is hapmap with haploid allele! Please Check Genotype data!")}

Gx <- geno[,1:11]
Gy <- geno[,12:ncol(geno)]
Gystr <- as.character(Gy[1,])
Gy <- Gy[,order(Gystr)]
geno <- cbind(Gx,Gy)

GT= t(geno[1,-(1:11)])   # 获得性状数据框
colnames(GT)="Taxa"   # 性状数据框列名修改
GT <- sub(":.*","",GT)   ### 此处可以修改，简化群体名称
Marker= geno[-1,c(1,3,4)]   # 获得标记数据框 
colnames(Marker)=c("SNP","Chromosome","Position")   # 修改标记数据框列名
GD= as.matrix(geno[-1,-(1:11)])   # get genotype matrix   获得基因型矩阵

## 数字转换部分
GD2 = NULL   # 初始化GD2变量
system.time({
  GD2<-apply(GD, 1, genoImputeM)   # 数字转换
})
GD2<-matrix(as.numeric(as.matrix(GD2)),nrow=nrow(GD2))   # 转为数值型
GD3<- GD2
GD3 <- as.data.frame(GD3) 

Marker2 <- t(Marker[,1])
names(GD3) <- Marker2
GD3 <- cbind(GT,GD3)
#system.time(write.table(GD3, paste(fileInfo[2],".NumericGenoForGS.csv",sep=""), 
#            quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE, na=""))
write_csv(GD3,paste0(fileInfo[1],".NumericGenoForGS.csv"),na="")

geno <- GD3
genoID <- as.matrix(geno[,1])   # 获取材料名称
geno <- geno[,-1]   # 去掉材料名
geno <- as.matrix(geno)   # 变为矩阵
row.names(geno) <- genoID   # 获得行名

geno=2*geno   #? 如果矩阵是0,0.5,1，则变为0,1,2
freqNA <- colSums(is.na(geno))/nrow(geno)   # NA频率
hist(freqNA)   # NA直方图
summary(freqNA)   # NA基本信息

# 补缺失
set.seed(123)   # 种子点
index=freqNA<0.15   #? 设置缺失值筛选条件
geno=geno[,index]   # 筛选缺失值
system.time(X <- impute_prob(geno))   # 补缺

# 检查是否还有NA
if(any(is.na(X))) {stop("")}else{cat("No NA! Pass!\n")}

# 计算最小等位基因频率
phat <- colMeans(X)/2   # 计算每列的均均值
maf <- ifelse(phat<0.5,phat,1-phat)   # 计算最小等位基因频
hist(maf,xlab="Minor allele frequency")   # 绘制最小等位基因频率图
summary(maf)   # MAF基本信息

# 根据MAF筛选
index <- maf>0.05   #? 获得MAF筛选索引
X <- X[,index]   # 按MAF筛选基因型数据
maf <- maf[index]   # 按MAF筛选maf

# 检查NA情况
if(any(is.na(X))) stop("You still have missing values\n")

# 遗传相关矩阵G，用标记的方差-协方差矩阵表示，=ZZ'
Z <- scale(X,center=TRUE,scale=TRUE)   # 筛选后的数据进行中心化和标准化
G <- tcrossprod(Z)/ncol(Z)   # 标记矩阵自乘，变成方形矩阵,降维


# PCA
pcageno <- prcomp(G,scale=TRUE)
str(pcageno$rotation)


if (length(fgroup) == 0){
  group <- data.frame(Gystr,rep(1,nrow(G)))
}else{
  group <- read.table(fgroup,header=F)
}

pcagroup <- as.character(group[,2])

mycolor <- c("black","#E94B35","#2C97De","#00BD9C","9C56B8")
mypch <- c(21,22,23,24,25)

colour_group <- mycolor[1:length(unique(pcagroup))]
pch_group <- mypch[1:length(unique(pcagroup))]

colour<-colour_group[as.numeric(factor(pcagroup))]
pch <- pch_group[as.numeric(factor(pcagroup))]

pve <- pcageno$sdev^2/sum(pcageno$sdev^2)
plot(pcageno$x[,1:2],col=colour,pch=pch,xlab=paste0("PC1 (",round(pve[1],2)*100,"%)"),ylab=paste0("PC2 (",round(pve[2],2)*100,"%)"))

if (length(fgroup)!=0){
  legend("topleft", pch=pch_group, horiz=TRUE, bty="n",col=colour_group, legend=unique(pcagroup)) 
}
