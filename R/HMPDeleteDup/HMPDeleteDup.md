# HMPDeleteRep

[TOC]

**制作：张敖 [https://datahold.cn](https://datahold.cn/)** 

# 鸣谢

张学才，CIMMYT，指导

## 最新更新

【20220330 15:12:30】

1. 整合了计算重复性功能。
2. 增加了两个数据合并的功能，但很耗时。

【2021-06-19 20:21:44】

1. 功能完成。

## 功能

测序时，可能会设置相同的材料测定两次，以检测测序的稳定性。这些材料会对分析产生影响，需要先进行处理。

1. HMP基因型文件，长文件名的缩短，如：GBPSVD22:2-B-B-3缩短为GBPSVD22<可选>。
2. 找到重复的材料名，剔除其中一个或将两个材料的结果合并成一个（很耗时）。
3. 自定义新的HMP文件名，自动补全.hmp.txt，只需提供前面的部分。
4. 生成重复材料列表，并计算重复性，计算公式为$r=v/v_{a}$，这里$r$表示重复性，$v$表示重复碱基数量，$v_a$表示总碱基数量，有缺失的数据不计算在内。
5. 生成TASSEL可打开的HMP文件。

## 使用方法

1. 计算全过程，请保证网络通畅。
2. 将下列代码复制到Rstudio，然后全选运行。期间会弹出对话框2次，第一次选择工作目录，第二次选择基因型文件（*.hmp.txt）。

海外：

```R
fName <- "newMyGeno"   # 修改为自己的文件名称
subtext <- ":.*|-.*"   # 缩短材料名，去除无用信息用，多个处理用|隔开
missingSign <- "N"   # 缺失值符号，有些HMP格式的文件中用字母N表示缺失
chooseOne <- TRUE   # 默认TRUE，当两个材料测定两次时，保留一个材料；FALSE表示合并两次测定，但会消耗很长的时间。
source("https://aozhangchina.github.io/R/HMPDeleteDup/HMPDeleteDupV2.r")   # 加载程序文件，需要联网
```

国内：

```R
fName <- "newMyGeno"   # 修改为自己的文件名称
subtext <- ":.*|-.*"   # 缩短材料名，去除无用信息用，多个处理用|隔开
missingSign <- "N"   # 缺失值符号，有些HMP格式的文件中用字母N表示缺失
chooseOne <- TRUE   # 默认TRUE，当两个材料测定两次时，保留一个材料；FALSE表示合并两次测定，但会消耗很长的时间。
source("https://dataholdcn.cn/R/HMPDeleteDup/HMPDeleteDupV2.r")   # 加载程序文件，需要联网
```
