## ViolinPlotTool

[TOC]

**制作：张敖 https://datahold.cn**

## 鸣谢

张学才，CIMMYT，指导

宋俊乔，安阳市农科院，测试

任姣姣，新疆农业大学，测试

## 最新更新

【2020-08-17】

增加国内线路。

【2020-07-16】

增加可以为每一个子图设置颜色。

增加可以自定义Y轴的范围。

修正了Xlab、Ylab、plotcolor设置可能无效的BUG。

【2019-10-09】

编写完成，可多个数据放在同一张图片中。

## 功能

1. 多个数据放在一个文件中，并画在一张图片里。
2. 可以自定义X轴标题和Y轴标题。
3. X轴刻度由数据中的第一行决定。
4. 可以自定义Y轴的范围。
5. 可以自定义每个子图的颜色，通过c()函数添加，多个颜色用,分开。如：c("lightgreen", "lightblue")
6. 操作简单，运算快。

## 生成图展示

![1570652098079](image\1570652098079.png)

## 准备输入文件

用Excel将所有数据组织成下表形式，第一行是数据标题，也就是图片横轴刻度名。数据长度可以不一致，画图时，多余的部分会自动删除掉。数据准备好后，另存为TXT格式。

| AT       | TL       | AL       |
| -------- | -------- | -------- |
| 52.5446  | 64.1877  | 52.55286 |
| 57.32    | 43.20504 | 49.22077 |
| 48.08599 | 41.81292 | 45.50928 |
| 51.66822 | 52.28344 | 49.33801 |
| 52.05925 | 38.60952 | 46.14632 |
| 47.37684 | 40.49051 | 44.89529 |
| 53.51524 | 44.1197  | 48.00747 |
| 51.31609 | 41.56952 | 46.59532 |
| 52.40458 | 44.64032 | 47.70334 |

## 使用方法

将下列代码复制到RStudio，根据需要，自行修改下面代码的第2到第4行，全选运行。

海外：

```r
rm(list=ls())   # remove objects
Xlab <- "Locations"   #? X轴标题
Ylab <- "SRF sorse"   #? Y轴标题
Ylim <- c(35,70)   #? Y轴范围，程序自动计算请写NULL
plotcolor <- "gray"   #? 颜色 "tomato" c("lightgreen", "lightblue")
source("https://aozhangchina.github.io/R/Violin/Vioplot.R", encoding = "utf-8")   # 加载程序文件，需要联网
```

```r
source("https://aozhangchina.github.io/R/Violin/Vioplot.R")   # 非RStudio运行
```

国内：

```R
rm(list=ls())   # remove objects
Xlab <- "Locations"   #? X轴标题
Ylab <- "SRF sorse"   #? Y轴标题
Ylim <- c(35,70)   #? Y轴范围，程序自动计算请写NULL
plotcolor <- "gray"   #? 颜色 "tomato" c("lightgreen", "lightblue")
source("https://dataholdcn.cn/R/Violin/Vioplot.R", encoding = "utf-8")   # 加载程序文件，需要联网
```

```R
source("https://dataholdcn.cn/R/Violin/Vioplot.R")   # 非RStudio运行
```



> 程序运行时，会提示选择数据文件，请选择之前准备好的文件。
>
> 程序自动运行直至全部完毕。

