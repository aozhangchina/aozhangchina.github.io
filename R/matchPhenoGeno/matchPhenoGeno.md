# matchPhenoGeno

[TOC]

**制作：张敖 [https://datahold.cn](https://datahold.cn/)** 

# 鸣谢

张学才，CIMMYT，指导

## 最新更新

【2021-06-20 12:55:41】

修正文件名中材料数不正确的Bug。

【2021-06-19 23:44:23】

功能完成。

## 功能

1. 获得同时拥有基因型和表型的基因型数据和表型数据。
2. 可以根据需要缩短材料名。

## 使用方法

1. 计算全过程，请保证网络通畅。
2. 将下列代码复制到Rstudio，然后全选运行。期间会弹出对话框3次，第一次选择工作目录，第二次选择基因型文件（*.hmp.txt），第三次选择表型文件。

海外：

```R
fName.geno <- "myGeno"   # 修改为要生成的基因型文件名
fName.pheno <- "myPheno META-R import"   # 修改为要生成的表型文件名
ID <- "GenoID"   # 设置基因型和表型共有的ID名称，即表型文件中的列名
subtext <- ":.*|-.*"   # 缩短材料名，去除无用信息用，多个处理用|隔开
source("https://aozhangchina.github.io/R/matchPhenoGeno/matchPhenoGeno.r")   # 加载程序文件，需要联网
```

国内：

```R
fName.geno <- "myGeno"   # 修改为要生成的基因型文件名
fName.pheno <- "myPheno META-R import"   # 修改为要生成的表型文件名
ID <- "GenoID"   # 设置基因型和表型共有的ID名称，即表型文件中的列名
subtext <- ":.*|-.*"   # 缩短材料名，去除无用信息用，多个处理用|隔开
source("https://dataholdcn.cn/R/matchPhenoGeno/matchPhenoGeno.r")   # 加载程序文件，需要联网
```

