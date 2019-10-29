## LDdecayPlotTool

[TOC]

**制作：张敖** **[https://datahold.cn](https://datahold.cn/)** 

## 鸣谢

张学才，CIMMYT，指导

曹士亮 ，黑龙江农科院，测试玉米基因组

袁一冰，四川省农科院，测试玉米基因组

宋俊乔，安阳市农科院，测试玉米基因组

Rahul， Indian Agricultural Research Institute，测试小麦基因组

## 最新更新

【2019-10-29】

增加了对小麦基因组的支持。理论上，适合任何基因组的绘制。可能需要多次调节参数。

增加了错误返回值，用于判断如何设置LD阈值。

【2019-10-08】

修正了有时刻度显示不全的BUG。

修正了有时单位显示错误的BUG。

【2016-11-23】

完善x轴刻度。

【2016-11-21】

可选择单位kb或bp。

在图像中心输出r^2^值。

去掉ggplot2包，完全使用基础绘图功能，提高兼容性。

完善细节。

【2016-11-09】

在图中标出LD平均LD衰减值。

为每一个chr计算LD衰减值，并输出在表LDDecayChr.csv中。

去掉分染色体LD衰减图的背景及网络线。

【2016-11-02】

增加简单线性方程计算，获得比较精确的LD衰减值，仅供参考。

【2016-10-31】

增加分染色体画图。

增加图例。

【2016-10-26】

根据衰减阈值，画出局部衰减图。

增加参数自定义。

【2016-10-21】

雏形完成，直接读取TESSEL的LD分析文件进行画图。

## 功能

1. 直接利用TASSEL的LD分析生成的文件
2. 生成整体的LD衰减图（自行保存为PDF、TIFF、JPG、PNG等格式）
3. 生成每条染色体及总体均值的LD衰减图（可自行保存为PDF、TIFF、JPG、PNG等格式）
4. 生成可供Excel画图的表格（文件名：群体名称+LD Decay Plot .csv）
5. 快速读取大文件，1.4G文件1min内读取完毕。
6. 根据简单线性模型，估算比较精确的LD衰减值。

## 生成图展示

### 均值图像

  ![DTMA LD Decay Plot](image\DTMA LD Decay Plot.png)
  ###分染色体图像（显示r^2^）
  ![DTMA LD Decay Plot chr](image\DTMA LD Decay Plot chr.png)

### 分染色体图像（不显示r^2^）

![DTMA LD Decay Plot chr without r2](image\DTMA LD Decay Plot chr without r2.png)

## 生成表展示

| Scale Distance  (kb) | Chr1 | Chr2 | Chr3 | Chr4 | Chr5 | Chr6 | Chr7 | Chr8 | Chr9 | Chr10 | Mean |
| -------------------- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ----- | ---- |
| 0                    | 1    | 1    | 1    | 1    | 1    | 1    | 1    | 1    | 1    | 1     | 1    |
| 0.5                  | 0.36 | 0.35 | 0.35 | 0.38 | 0.36 | 0.35 | 0.34 | 0.36 | 0.36 | 0.36  | 0.36 |
| 1                    | 0.16 | 0.14 | 0.15 | 0.15 | 0.14 | 0.14 | 0.14 | 0.16 | 0.15 | 0.16  | 0.15 |
| 2                    | 0.13 | 0.11 | 0.12 | 0.13 | 0.12 | 0.11 | 0.11 | 0.13 | 0.12 | 0.13  | 0.12 |
| 5                    | 0.11 | 0.1  | 0.1  | 0.1  | 0.09 | 0.09 | 0.09 | 0.11 | 0.1  | 0.1   | 0.1  |
| 10                   | 0.08 | 0.07 | 0.07 | 0.08 | 0.07 | 0.07 | 0.07 | 0.08 | 0.08 | 0.08  | 0.08 |
| 20                   | 0.07 | 0.06 | 0.07 | 0.06 | 0.06 | 0.06 | 0.06 | 0.08 | 0.06 | 0.07  | 0.06 |
| 30                   | 0.07 | 0.05 | 0.06 | 0.06 | 0.05 | 0.06 | 0.06 | 0.06 | 0.05 | 0.06  | 0.06 |
| 50                   | 0.06 | 0.05 | 0.06 | 0.06 | 0.05 | 0.05 | 0.05 | 0.06 | 0.06 | 0.06  | 0.06 |

| Chromosome | LD Decay Distance (kb) |
| ---------- | ---------------------- |
| Chr1       | 6.28                   |
| Chr2       | 4.55                   |
| Chr3       | 5.73                   |
| Chr4       | 5.12                   |
| Chr5       | 4.17                   |
| Chr6       | 3.74                   |
| Chr7       | 3.54                   |
| Chr8       | 6.62                   |
| Chr9       | 4.55                   |
| Chr10      | 4.91                   |
| Mean       | 4.91                   |

## 使用方法

推荐将所有文件（包括r文件和LD衰减文件）放在同一个目录下，然后在R Studio中运行下列代码。

```r
rm(list=ls())   # Clear the Environment Variables 清空环境变量 
coordinateScale <- c(0, 500, 1000, 2000, 5000, 10000, 20000, 30000, 50000)   # x coordinate axis, unit is bp 坐标x轴,按照pb填写
unit <- "Kb"   # Units on the image: "Kb" or "bp" 图像单位
pop.name <- "DTMA"    # The population's name 群体名称
plotTitle <- "LD Decay Plot"   # The title of plot  图表标题
lineWidth <- 2   # the width of line  线宽
lineColor <- "blue"   # The color of line 线颜色
#lineColor <- "blue","brown","green","gray","gold","darkblue","lightblue","red","orange","tan","yellowgreen")   # The color of line 线颜色
dot <- "o"   # "o" is true, "l" is false 线型："l"是折线，"o"是带点的折线
dottedLine <- TRUE # TRUE or FALSE 
dividingLine <- c(0.1)   # 0.1 or 0.2 or c(0.1,0.2)
plotChrR2 <- TRUE   # R^2 value display 全染色体图例中是否显示r^2值，TURE：显示，FLASE：不显示
source("https://aozhangchina.github.io/R/LDdecay/LDdecayV4.r",encoding = "utf-8")   # 加载程序文件，需要联网
```

> 然后选择LD数据文件，如：allLD.csv。
>
> 程序自动运行直至全部完毕。
>
> 最后一行会显示估算的比较精确的LD衰减值，并显示两个图像。

