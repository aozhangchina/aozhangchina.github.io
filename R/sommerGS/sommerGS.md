# SommerGS

 

[TOC]

**制作：张敖、陈珊、阮燕晔 [https://datahold.cn](https://datahold.cn/)** 

# 鸣谢

张学才，CIMMYT，指导

候智馨，沈阳农业大学，测试

# 公开测试

该版本公开测试，可能有很多BUG，如有问题请联系作者。

## 最新更新

【2023-11-03 11:00:29】

1. 修正了有时avgGRM无法计算预测精度的问题。

【2023-03-03 15:53:29】

1. 公开测试1.0版本。
2. 整合新版本的GBIT。

【2021-09-07 22:24:40】

1. 修正了不适用CDmean出错的BUG。

【2021-07-13 16:23:55】

1. 修正了CDmean的群体大小修改无效的BUG。
2. 修正了有时rrBLUP方法输出文件错误的BUG。

【2021-07-11 22:21:58】

1. 修正了avgGRM的路径问题。
2. 增加了预测群体名到avgGRM和CDmean方法的预测精度文件名。

【2021-07-11 12:48:06】

1. 增加了利用已经分析好的基因型和表型中间文件跳过前面若干步骤，大幅提高计算效率。
2. 增加CDmean和avgGRM的预测精度。

【2021-07-08 18:03:20】

1. 增加CDmean方法的群体优化。
2. 增加avgGRM方法的群体优化。

【2021-07-04 11:24:03】

1. 修正BUG。

【2021-06-29 22:33:45】

1. 完成功能。
2. 当sommer包的rrBLUP无法使用时，可以选择使用rrBLUP包的rrBLUP建模算法。

## 功能

1. HMP数值化。
4. 利用sommer包的GBLUP和rrBLUP方法计算全基因组预测值。
3. 将环境效应作为固定效应，做GS。
4. 自定义k倍交叉验证的k值。
5. 整合了rrBLUP包，可以使用rrBLUP包做预测。
6. 使用CDmean和avgGRM方法优化建模群体。

## 使用方法

### 普通预测

1. 准备基因型文件，HMP格式，推荐先使用TASSEL进行表型和基因型的筛选。

2. 准备表型文件，可以直接使用CIMMYT开发的META-R软件的**跨环境**结果作为表型文件。

3. 多个群体组合在一起可以增加group列，用于区分每个群体分组，做群体间预测。

4. ```R
   Env <- NULL   # 这项必须设置为NULL
   NTrn.Optmize <- NULL   # 这项必须设置为NULL
   ```
   
5. 结果会生成在【GBIT/NoGbyE】目录下。

### 群体A预测群体B

1. 设置如下内容。

   ```R
   A_pred_B <- TRUE   # A群预测B
   ```

2. 其他设置同普通预测。

### 基因型与环境互作预测

#### 环境效应作为固定效应

1. 该预测假设，每个环境须有相同的方差。

2. 表型文件使用META-R软件对个体计算BLUP、BLUE的结果，或类似格式。

3. 若要使用不同群体进行预测，首先要保证每个群体都有相同的地点。

4. ```R
   Env <- "Environment"   # 表型文件中对应列的名称
   GbyEMain <- TRUE   # 使用sommer包，环境作为固定效应
   group <- c("DH1","DH2","DH3","CML","DTMA","SYNDH")   # 根据要求选择群体
   	train_pop <- "SYNDH"   # 指定建模群体
   	test_pop <- "DTMA"   # 指定验证群体
   ```

5. 结果会生成在【GBIT/GbyE】目录下。

6. 每个循环的预测结果在【GbyE_GEBV.csv】中，可以手动计算相关性（预测精度）。

7. 其他设置同普通预测。

#### 环境效应方差不同（完善中）

1. 该预测假设，每个环境的方差不同，在每个环境拟合方差组分。缺点是环境间没有协方差。

2. 表型文件使用META-R软件对个体计算BLUP、BLUE的结果，或类似格式。

3. ```R
   Env <- "Environment"   # 表型文件中对应列的名称
   GbyEDG <- TRUE   # 使用sommer包，对角线模型
   group <- c("DH1","DH2","DH3","CML","DTMA","SYNDH")   # 根据要求选择群体
   	train_pop <- "SYNDH"   # 指定建模群体
   	test_pop <- "DTMA"   # 指定验证群体
   ```

4. 结果会生成在【GBIT/GbyE】目录下。

5. 每个循环的预测结果在【GbyE_GEBV.csv】中，可以手动，见

6. 其他设置同普通预测。

### CDmean建模群体优化

1. 设置如下内容。

   ```R
   A_pred_B <- TRUE   # A群预测B群
   optimization <- "CDmean"   # 建模群体优化方法
   NTrn.Optmize <- 100   # 建模群体大小，CDmean方法有效
   group <- c("DH1","DH2","DH3","CML","DTMA","SYNDH")   # 多个群体合并成大群体时，可以利用该列进行分组，如填写c("SYNDH","DTMA")
   test_pop <- "CML"   # 非抽取情况下指定验证群体
   ```

2. 如果有错误信息，请忽略，直接查看【GBLUP\_方法\_群体大小.cor.csv】和【rrBLUPr\_方法\_群体大小.cor.csv】文件中的相关系数，出错时，通常GBLUP的结果为NA，所有预测结果均为0，原因不明，推测是由于sommer包使用的直接反演算法，有局限性。此时可参考rrBLUP的结果，两者实际预测结果非常接近。

### avgGRM建模群体优化

1. 设置如下内容。train_pop里不能包含test_pop。

   ```R
   A_pred_B <- TRUE   # A群预测B群
   optimization <- "avgGRM"   # 建模群体优化方法
   NTrn.Optmize <- 100   # 建模群体大小，CDmean方法有效
   group <- c("DH1","DH2","DH3","CML","DTMA","SYNDH")   # 多个群体合并成大群体时，可以利用该列进行分组，如填写c("SYNDH","DTMA")
   test_pop <- "CML"   # 非抽取情况下指定验证群体
   ```

2. 运行后会获得【optimization\avgGRM_DTMA.csv】文件，`_DTMA`表示测试群体为DTMA。

3. 再次运行，不会重复生成该文件，直接读取。

4. 如果有错误信息，请忽略，直接查看【GBLUP\_方法\_群体大小.cor.csv】和【rrBLUPr\_方法\_群体大小.cor.csv】文件中的相关系数，出错时，通常GBLUP的结果为NA，所有预测结果均为0，原因不明，推测是由于sommer包使用的直接反演算法，有局限性。此时可参考rrBLUP的结果，两者实际预测结果非常接近。

### 计算相关性

1. 打开【GbyE_GEBV.csv】文件，全选并复制。

2. 运行下面代码。

   ```R
   source("https://dataholdcn.cn/R/GBIT/GBIT.R")
   setwd(choose.dir())   # 选择GbyE_GEBV.csv所在目录
   GBIT.library("reshape2")
   a <- read.table("clipboard",header=F)   # 需要先在csv文件中复制所有内容
   b <- dcast(a,V1~V3,value.var = "V2")
   write.table(b,"GbyE_BLUE.csv",sep=",",row.names = F,col.names = F
   ```

3. 运行完成会生成【GbyE_BLUE.csv】，另存为【GbyE_BLUP.xlsx】。

4. 将BLUP或BLUE数据拷贝到该数据旁，最好空几行。

5. 按照材料名称完全对应，可使用排序功能，多余材料删除。

6. 利用函数`=CORREL(B1:B282,DC1:DC282)`计算相关性。

## 代码

1. 运行代码前，请按照说明调整参数，否则会出错。
2. 将下面代码复制到RStudio中，修改参数。
3. WD设置后，会生成【GBIT】和【optimization】。
4. 全选运行。

国内：

```R
rm(list=ls())   # 清空内存
traitName <- "BLUE"   # 表型文件中，数据的列名
taxa_name_in_Pheno <- "Genotype"   # 表型文件中，基因型名称的列名
Env <- NULL   # 计算GbyE时，为表型文件的的环境名称，如"Environment"
GBLUP <- TRUE   # 使用sommer包的GBLUP模型
rrBLUP <- FALSE   # 使用sommer包的rrBLUP模型，数据量大时会出错
rrBLUPr <- TRUE   # 使用rrBLUP包的rrBLUP模型
GbyEMain <- FALSE   # 使用sommer包，环境作为固定效应
GbyEDG <- FALSE   # 使用sommer包，环境对角线模型
cycles <- 20   # 随机抽取建模群体预测验证群体的重复次数
get.train_test <- FALSE   # 生成每次抽取的建模群体和验证群体材料名
A_pred_B <- FALSE   # A群体预测B群体，即确定了建模群体和预测群体时用TRUE
group <- NULL   # 多个群体合并成大群体时，可以利用该列进行分组，如填写c("SYNDH","DTMA")
	train_pop <- "SYNDH"   # 非抽取情况下，指定建模群体
	test_pop <- "DTMA"   # 非抽取情况下，指定验证群体
K_fold <- 5   # K倍交叉验证
optimization <- FALSE   # FALSE，可选择"CDmean"或"avgGRM"，一次只能运行一个
	NTrn.Optmize <- 20   # 建模群体大小，CDmean、avgGRM方法有效
	nIter <- 10000   # 迭代次数，CDmean方法有效
WD <- NULL   # 工作目录，若为NULL，则在程序运行时选择，否则为目录。"E:\\data\\FER data\\results_za\\Raw_phenotyic_data_and_analysis\\data Geno and Pheno\\GS\\SYNDH\\oneLocation"
sommerGenof <- NULL   # HMP格式的基因型数据，为NULL则在运行时选择。"E:\\data\\FER data\\results_za\\Raw_phenotyic_data_and_analysis\\data Geno and Pheno\\Genotipic data\\6Pop_1190_121135.hmp.txt"
sommerPhenof <- NULL   # 表型文件，格式参照META-R的结果，NULL表示程序运行时选择。"E:\\data\\FER data\\results_za\\Raw_phenotyic_data_and_analysis\\data Geno and Pheno\\Meta-R analysis\\Output_META-R\\SYNDH\\TL16A.csv"
source("https://dataholdcn.cn/R/sommerGS/sommerGS.R")   # 加载程序文件，需要联网
```

