## getICIMappingFiles

[TOC]

**制作：张敖** **[https://datahold.cn](https://datahold.cn/)** 

## 鸣谢

张学才，CIMMYT，指导

孙恺悦，沈阳农业大学，测试

## 最新更新

【20220401 18:12:34】

增加.bip文件的生成。

【20220331 09:41:31】

基本功能实现。

## 功能

QTL IciMapping是中国农业科学院开发的优秀基因定位软件，本程序通过处理HMP文件，快速获得能够导入QTL IciMapping的.snp文件；通过处理HMP和表型文件，快速获得.bip文件。HMP文件需要现在TASSEL软件中进行筛选，去掉+/-和0号染色体。

1. 本工具需要使用TASSEL筛选后的基因型。
2. 本工具只能识别双等位基因，第三等位基因会被转变为NA。
3. 生成的.snp文件可以导入QTL IciMapping软件，用于绘制连锁图谱，具体过程请参看IciMapping说明书。
4. 本程序只适合于RIL和DH群体，不适用于其他类型的群体！
5. 用到的表型文件需要使用META-R软件生成的格式。

## 使用方法

准备好输入文件，然后在R Studio中运行下列代码。第一个弹出的对话框选择工作目录，第二个选择HMP格式的基因型文件。

海外：

```r
rm(list=ls())   # Clear the Environment Variables 清空环境变量
first_change <- FALSE   # The first cycle to change letters. 两步变换的第一步，默认不开启，数据量大时非常耗时。无特殊需求无需打开
Eliminate_heterozygous_bases <- TRUE   # the default is TRUE, the RIL and DH population type have to use TRUE. RIL群体和DH群体必须选择TRUE
f_name <- "my_icimapping_data"   # Customaize the name of output file. 自定义输出文件名称。
pheno_colname <- "BLUE_FV_"   # The name of trait in phenotypic file. 表型文件中性状的名称。
source("https://aozhangchina.github.io/R/getICIMappingFiles/getIciMappingFiles.R")   # 加载程序文件，需要联网
```

国内：

```R
rm(list=ls())   # Clear the Environment Variables 清空环境变量 
first_change <- FALSE   # The first cycle to change letters. 两步变换的第一步，默认不开启，数据量大时非常耗时。无特殊需求无需打开
Eliminate_heterozygous_bases <- TRUE   # the default is TRUE, the RIL and DH population type have to use TRUE. RIL群体和DH群体必须选择TRUE
f_name <- "my_icimapping_data"   # Customaize the name of output file. 自定义输出文件名称。
pheno_colname <- "BLUE_FV_"   # The name of trait in phenotypic file.
source("https://dataholdcn.cn/R/getICIMappingFiles/getIciMappingFiles.R")   # 加载程序文件，需要联网
```

