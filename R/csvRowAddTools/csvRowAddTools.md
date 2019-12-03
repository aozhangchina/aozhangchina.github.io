# csvRowAddTools

[toc]

**制作：张敖 [https://datahold.cn](https://datahold.cn/)**

## 最新更新

【2019-12-3】

编写完成，一次性选择多个文件合并。

## 功能

1. 将列相同的多个csv文件合并成一个文件。
2. 注意：只能按行合并，请确保每一列都是数据都能对应上。

例如，csv文件1：

| Loc  | GID     | BLUP_Fe  | BLUE_Fe  | pedigree                       |
| ---- | ------- | -------- | -------- | ------------------------------ |
| 1    | 1005209 | 17.28632 | 16.45146 | CLWQZ1716//CLWQHZN19/CLWQHZN77 |
| 1    | 1005212 | 18.2155  | 18.25768 | CLWQZ1710//CML503/CML557       |

csv文件2：

| Loc  | GID     | BLUP_Fe  | BLUE_Fe  | pedigree                       |
| ---- | ------- | -------- | -------- | ------------------------------ |
| 2    | 1005209 | 17.15545 | 16.91246 | CLWQZ1716//CLWQHZN19/CLWQHZN77 |
| 2    | 1005212 | 18.48736 | 18.59424 | CLWQZ1710//CML503/CML557       |

合并后：

| Loc  | GID     | BLUP_Fe  | BLUE_Fe  | pedigree                       |
| ---- | ------- | -------- | -------- | ------------------------------ |
| 1    | 1005209 | 17.28632 | 16.45146 | CLWQZ1716//CLWQHZN19/CLWQHZN77 |
| 1    | 1005212 | 18.2155  | 18.25768 | CLWQZ1710//CML503/CML557       |
| 2    | 1005209 | 17.15545 | 16.91246 | CLWQZ1716//CLWQHZN19/CLWQHZN77 |
| 2    | 1005212 | 18.48736 | 18.59424 | CLWQZ1710//CML503/CML557       |

## 使用方法

将下列代码复制到Rstudio，请按需求修改下面第一行，然后全选运行。期间会弹出对话框，一次性选择所有需要合并的文件文件即可。

```r
source("https://aozhangchina.github.io/R/csvRowAddTools/csvRowAddTools.r")   # 加载程序文件，需要联网
xx1 <- combineFile()
write.csv(xx1,paste0(file=choose.dir(),"/","GY_11&12_2017.csv"),row.names = F)   #? 请修改【GY_11&12_2017.csv】为你自己的文件名。
```

