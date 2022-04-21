## ChromesomeHeatmapTool

[toc]

**制作：张敖 [https://datahold.cn](https://datahold.cn/)**

## 鸣谢

张学才，CIMMYT，指导

## 最新更新

【2022-04-21】

增加了自定义染色体大小和着丝粒位置功能。

增加了根据基因型数据自动适应染色体大小功能。

增加了白边修正功能，出现白边的可能性降低。

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

3. 可以自定义染色体长度和着丝粒位置（需准备info_chr.txt），若无需自定义，则在选择文件时点取消。无需修改的部分写NA即可。

   | Chr  | Chr_strat | Ch_end    | CE_ start | CE_end |
   | ---- | --------- | --------- | --------- | ------ |
   | 1    | NA        | 307041717 | NA        | NA     |
   | 2    | NA        | 244442276 | NA        | NA     |
   | 3    | NA        | 235667834 | NA        | NA     |
   | 4    | NA        | 246994605 | NA        | NA     |
   | 5    | NA        | 223902240 | NA        | NA     |
   | 6    | NA        | 174033170 | NA        | NA     |
   | 7    | NA        | 182381542 | NA        | NA     |
   | 8    | NA        | 181122637 | NA        | NA     |
   | 9    | NA        | 159769782 | NA        | NA     |
   | 10   | NA        | 150982314 | NA        | NA     |

4. 暂不支持其他基因组的绘制。

## 生成图片预览

![](img\chromosome.jpg)

## 使用方法

将下列代码复制到Rstudio，根据需要，自行修改下面代码的第二行，全选运行。

海外：

```R
rm(list=ls())
plotcolor <- c("#4575b4", "#ffffbf", "#d73027")   #? 颜色 c("#e5f5f9", "#99d8c9", "#2ca25f")
chr_size <- "value"   # "value"表示给定的值，即参考基因组给定大小；"auto"表示根据基因型自己设置值。
source("https://aozhangchina.github.io/R/chromesomeheatmapTool/chromesomeheatmap.r")   # 加载程序文件，需要联网r
```

国内：

```R
rm(list=ls())
plotcolor <- c("#4575b4", "#ffffbf", "#d73027")   #? 颜色 c("#e5f5f9", "#99d8c9", "#2ca25f")
chr_size <- "value"   # "value"表示给定的值，即参考基因组给定大小；"auto"表示根据基因型自己设置值。
source("https://dataholdcn.cn/R/chromesomeheatmapTool/chromesomeheatmap.r")   # 加载程序文件，需要联网r
```

