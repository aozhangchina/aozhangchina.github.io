if (!exists('inputing_method')) {inputing_method <- "possibility"} #? ML/REML/possibility
if (!exists('method')) {method <- "percentage"}   #? "percentage"/"leave_one"/"k-fold"
if (!exists('trainingSize')) {trainingSize <- 0.5}   # training population size
if (!exists('cycles')) {cycles <- 10}   # 设置循环次数
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



#↑ 函数部分结束 ↑#

#↓ 分析部分开始 ↓#
## 加载rrBLUP
libraryAndInstall("rrBLUP")
## 加载dplyr
libraryAndInstall("dplyr")




























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
