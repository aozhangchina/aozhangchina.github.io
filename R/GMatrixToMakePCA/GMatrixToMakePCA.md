# GMatrixToMakePCA

 

[TOC]

**制作：张敖 [https://datahold.cn](https://datahold.cn/)** 

# 鸣谢

张学才，CIMMYT，指导

## 最新更新

【2020-02-27】

完成基本功能。

## 功能

1. HMP数值化。根据等位基因频率修改，最大等位基因频率1，最小等位基因频率0，杂合0.5。
2. 去除没有多态性的标记。
3. 根据可能性补缺失。
4. 可以选择是否分组，分组后，不同的组会用不同颜色的点和图形显示在图上。（最大支持5组）。
5. 生成PC1和PC2的图，括号里显示解释的变异百分比。

## 使用方法

准备hmp文件和分组文件（可选），分组文件有两列，第一列为材料名，第二列为分组标记。

将下列代码复制到Rstudio，然后全选运行。期间会弹出对话框，先选择hmp格式的基因型文件，再选择分组文件（**注意**：没有分组文件点取消）。

```
source("https://aozhangchina.github.io/R/GMatrixToMakePCA/GMatrixToMakePCA.r")   # 加载程序文件，需要联网
```
