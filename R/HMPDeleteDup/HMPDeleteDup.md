# HMPDeleteRep

[TOC]

**制作：张敖 [https://datahold.cn](https://datahold.cn/)** 

# 鸣谢

张学才，CIMMYT，指导

## 最新更新

【2021-06-19 20:21:44】

功能完成。

## 功能

1. HMP基因型文件，长文件名的缩短，如：GBPSVD22:2-B-B-3缩短为GBPSVD22<可选>。
2. 找到重复的材料名，剔除其中一个。
3. 自定义新的HMP文件名，自动补全.hmp.txt，只需提供前面的部分。
4. 生成重复材料列表。
5. 生成TASSEL可打开的HMP文件。

## 使用方法

1. 计算全过程，请保证网络通畅。
2. 将下列代码复制到Rstudio，然后全选运行。期间会弹出对话框2次，第一次选择工作目录，第二次选择基因型文件（*.hmp.txt）。

海外：

```R
fName <- "newMyGeno"   # 修改为自己的文件名称
source("https://aozhangchina.github.io/R/HMPDeleteDup/HMPDeleteDup.r")   # 加载程序文件，需要联网
```

国内：

```R
source("https://dataholdcn.cn/R/HMPDeleteDup/HMPDeleteDup.r")   # 加载程序文件，需要联网
```

