## GOMultiRowintoOneRow

[TOC]

**制作：张敖** **[https://datahold.cn](https://datahold.cn/)** 

## 鸣谢

赵成昊，沈阳农业大学，测试

## 最新更新

【2020-08-17】

增加国内线路。

【2020-01-02】

功能实现。

## 功能

1. 将GO数据多行数据改成一行形式。

2. 需要使用CSV导入，结果文件可在Excel中以逗号分列。

   转换前：

   | RGZM2G010871 | GO:0003700 |
   | ------------ | ---------- |
   |              | GO:0003701 |
   |              | GO:0003702 |
   | RGZM2G010872 | GO:0003707 |
   |              | GO:0003708 |

   转换后：
   
   | RGZM2G010871 | GO:0003700,GO:0003701,GO:0003702 |
   | ------------ | -------------------------------- |
   | RGZM2G010872 | GO:0003707,GO:0003708            |

## 使用方法

在R Studio中运行下列代码。

海外：

```r
source("https://aozhangchina.github.io/R/GOMultiRowintoOneRow/GOMultiRowintoOneRow.r")   # 加载程序文件，需要联网
```

国内：

```R
source("https://dataholdcn.cn/R/GOMultiRowintoOneRow/GOMultiRowintoOneRow.r")   # 加载程序文件，需要联网
```



> 运行后需要选择你要转换的CSV文件，实例文件下载：[【点此下载】](https://aozhangchina.github.io/R/GOMultiRowintoOneRow/test.csv)。
>
> 程序自动运行直至全部完毕，结果文件及生成目录会在console提示。
>

