# getTwoTaxaAccuracy 

[TOC]

**制作：张敖 [https://datahold.cn](https://datahold.cn/)** 

# 鸣谢

张学才，CIMMYT，指导

## 最新更新

【2021-06-01 18:59:29】

功能完成，新鲜发布。

## 功能

1. 测序时，有时会设置2个相同的材料用于评估测序的准确性。本程序计算相同碱基/总碱基座的比率（只计算非缺失位点）。

## 使用方法

1. 计算全过程，请保证网络通畅。

2. 将相同的材料列表保存成CSV文件。【[示例文件.csv](https://dataholdcn.cn/R/getTwoTaxaAccuracy/vlist.csv)】

   | Markers                      | GID    |
   | ---------------------------- | ------ |
   | GMPGBS24-P12:MRG:3:250142093 | 275823 |
   | GMPGBS24-P13:MRG:3:250142210 | 275823 |
   | GMPGBS25-P12:MRG:3:250142094 | 34867  |
   | GMPGBS25-P13:MRG:3:250142209 | 34867  |
   | GMPGBS541:MRG:3:250142096    | 279158 |
   | GMPGBS541:MRG:3:250142163    | 279158 |
   | GMPGBS542:MRG:3:250142097    | 279159 |
   | GMPGBS542:MRG:3:250142164    | 279159 |
   | GMPGBS543:MRG:3:250142098    | 279160 |
   | GMPGBS543:MRG:3:250142165    | 279160 |
   | GMPGBS544:MRG:3:250142212    | 279161 |
   | GMPGBS544:MRG:3:250142278    | 279161 |
   | GMPGBS545:MRG:3:250142213    | 279162 |
   | GMPGBS545:MRG:3:250142279    | 279162 |
   | GMPGBS546:MRG:3:250142214    | 279163 |
   | GMPGBS546:MRG:3:250142280    | 279163 |

3. 将下列代码复制到Rstudio，然后全选运行。期间会弹出对话框2次，第一次选择hmp.txt格式的基因型文件，第二次选择包含相同材料名的列表文件。

4. 程序运行完，结果在剪贴板，直接在需要的地方粘贴即可。第一列是ID，第二列是准确性。

海外：

```R
source("https://aozhangchina.github.io/R/getTwoTaxaAccuracy/getTwoTaxaAccuracy.r")   # 加载程序文件，需要联网
```

国内：

```R
source("https://dataholdcn.cn/R/getTwoTaxaAccuracy/getTwoTaxaAccuracy.r")   # 加载程序文件，需要联网
```

