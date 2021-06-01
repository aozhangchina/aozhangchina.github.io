# GCA_SCATool 

[TOC]

**制作：张敖 [https://datahold.cn](https://datahold.cn/)** 

# 鸣谢

张学才，CIMMYT，指导

## 最新更新

【2021-06-01 12:55:48】

功能完成，新鲜发布。

## 功能

1. 利用跨环境的BLUE或BLUP值计算一般配合力（GCA）和特殊配合力（SCA）。

## 使用方法

1. 计算全过程，请保证网络通畅。

2. 将表型数据整理成下表形式，多环境数据可以使用META-R计算BLUP或BLUE。【[示例文件.txt](https://aozhangchina.github.io/R/GCA_SCATool/demoGCA.txt)】

   | ID    | Tester1 | Tester2 | Tester3 | Tester4 | Tester5 |
   | ----- | ------- | ------- | ------- | ------- | ------- |
   | Line1 | 10.3    | 10.8    | 9.5     | 10.9    | 9.3     |
   | Line2 | 9.9     | 10.7    | 9       | 8.8     | 9.7     |
   | Line3 | 10.4    | 10.7    | 8.8     | 10.3    | 8.7     |

3. 将下列代码复制到Rstudio，然后全选运行。期间会弹出对话框2次，第一次选择工作目录（表型文件目录），第二次选择要计算一般配合力（GCA）和特殊配合力（SCA）的文件。

海外：

```
source("https://aozhangchina.github.io/R/GCA_SCATool/GCA_SCATool.r")   # 加载程序文件，需要联网
```

国内：

```R
source("https://dataholdcn.cn/R/GCA_SCATool/GCA_SCATool.r")   # 加载程序文件，需要联网
```

