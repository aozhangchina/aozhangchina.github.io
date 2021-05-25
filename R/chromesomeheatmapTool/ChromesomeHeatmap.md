## ChromesomeHeatmapTool

[toc]

**制作：张敖 [https://datahold.cn](https://datahold.cn/)**

## 鸣谢

张学才，CIMMYT，指导

## 最新更新

【2021-05-25】

修正了未安装RIdeogram无法绘图的BUG。

修正了染色体出现NA无法绘图的BUG。

【2020-08-17】

增加国内线路。

【2020-05-29】

编写完成，选择一个HMP格式的文件即可生成图片。

## 功能

1. 使用RIdeogram包画图，默认在输入文件所在目录生成【chromosome.svg】【chromosome.tiff】【chromosome.png】文件。
2. 内置玉米基因组基本信息，信息来自http://ensembl.gramene.org/Zea_mays/Location/Genome?r=10:1-1000。
3. 暂不支持其他基因组的绘制。

## 生成图片预览

![](img\chromosome.jpg)

## 使用方法

将下列代码复制到Rstudio，根据需要，自行修改下面代码的第二行，全选运行。

海外：

```R
rm(list=ls())
plotcolor <- c("#4575b4", "#ffffbf", "#d73027")   #? 颜色 c("#e5f5f9", "#99d8c9", "#2ca25f")
source("https://aozhangchina.github.io/R/chromesomeheatmapTool/chromesomeheatmap.r")   # 加载程序文件，需要联网r
```

国内：

```R
rm(list=ls())
plotcolor <- c("#4575b4", "#ffffbf", "#d73027")   #? 颜色 c("#e5f5f9", "#99d8c9", "#2ca25f")
source("https://dataholdcn.cn/R/chromesomeheatmapTool/chromesomeheatmap.r")   # 加载程序文件，需要联网r
```

