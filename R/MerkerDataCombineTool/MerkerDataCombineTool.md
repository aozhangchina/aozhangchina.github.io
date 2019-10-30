# MerkerDataCombineTool

 

[TOC]

**制作：张敖 [https://datahold.cn](https://datahold.cn/)** 

## 最新更新

【2019-10-30】

编写完成。 可将两个基因型文件整合成一个。

## 功能

1. 将两个HMP文件或其他标准基因型合成为一个文件，并对齐标记名。
2. 两个文件必须有相同的标记名和不同的材料个体名。
3. 可以选择是否生成TASSEL能打开的HMP文件，只针对HMP文件。

## 使用方法

将下列代码复制到Rstudio，请修改下面全选运行。

```
HMP <- TRUE   #? TRUE of FLASE, TRUE is generating a usable format for TASSEL
# Please choose two marker files once
source("https://aozhangchina.github.io/R/MerkerDataCombineTool/MerkerDataCombineTool.R")   # 加载程序文件，需要联网
```

