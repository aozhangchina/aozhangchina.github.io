## HMP2Flapjack

[TOC]

**制作：张敖** **[https://datahold.cn](https://datahold.cn/)** 

## 鸣谢

张学才，CIMMYT，指导

## 最新更新

【2020-06-24】

选择一个HMP的文件，自动生成Flapjack可以导入的MAP文件和GENOTYPE文件。

支持1.16.10.x及以后版本。

## 功能

1. 选择一个HMP的文件，自动生成Flapjack可以导入的MAP文件和GENOTYPE文件
2. 自动删除0号染色体的信息
3. 生成文件和源文件保持在相同目录

## 使用方法

在R Studio中运行下列代码。

```r
rm(list=ls())   # Clear the Environment Variables 
source("https://aozhangchina.github.io/R/HMP2Flapjack/HMP2Flapjack.r")
```

> 选择需要转换的HMP文件。

