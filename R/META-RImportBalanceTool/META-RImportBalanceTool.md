## META-RImportBalanceTool

[TOC]

**制作：张敖 https://datahold.cn**

## 鸣谢

张学才，CIMMYT，指导

陈珊，沈阳农业大学，测试

## 最新更新

【2020-07-21】

第一版完成。

## 功能

利用META-R计算跨环境的BLUE和BLUP时，经常会遇到闪退的情况。根据排查，是因为各个环境中的材料名称不一致所致。多环境田间调查经常会受到这样或那样的因素影响，导致一些本应调查的材料无法获得结果，简言之，该情况具有普遍性。本程序会自动对各个环境下的材料取交集，并最终生成META-R可以正确计算的输入文件。

## 准备输入文件

用Excel打开之前准备好的CSV文件。插入1行或2行，随便命名（本例中插入了一列，命名为k），并在下面的`使用方法`中修改对应的变量名称（列名）。

1. 增加`名称变量名`：合并ENTRY和REP，如下表第二行：NEC002+1=NEC0021，k列D2单元格的Excel公式为`=B2&C2`，下拉填充。
2. 增加`分类变量名`（可选）：修改如果环境包含YEAR、LOCATION和COUNTRY等信息，需要新建一列，并按照上面方式合并。本例中只有一个环境信息En，因此无需做处理。
3. 保存保存。

![snap_screen_20200721210743](images\snap_screen_20200721210743.png)

## 使用方法

将下列代码复制到RStudio，根据需要，自行修改下面代码的第2到第3行，要与输入文件的列名对应。然后全选运行。

```r
rm(list=ls())   # remove objects
var_env <- 'En'   # 分类变量名
var_material <- 'k'   # 名称变量名
source("https://aozhangchina.github.io/R/META-RImportBalanceTool/META-RImportBalanceTool.R")   # 加载程序文件，需要联网
```

> 程序运行时，会提示选择数据文件，请选择之前准备好的文件。
>
> 程序自动运行直至全部完毕。

