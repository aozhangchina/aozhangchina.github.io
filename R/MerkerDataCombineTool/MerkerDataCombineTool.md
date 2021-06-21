# MerkerDataCombineTool

 

[TOC]

**制作：张敖 [https://datahold.cn](https://datahold.cn/)** 

# 鸣谢

张学才，CIMMYT，指导

扈光辉，黑龙江农科院，测试

## 最新更新

【2021-06-20】

修正了生成的HMP有科学计数法而无法被TASSEL正常读取的问题。

【2020-11-25】

修正了有时合并只合并前十几行的BUG。

【2020-08-17】

增加国内线路。

【2020-02-25】

修正只产生一条数据的BUG。

增加显示每个材料的标记数量。

增加显示重复的材料名，用于手动删除。

增加可选项，delRep <- TRUE，删除重复材料（保留第一个）。

【2019-10-30】

编写完成。 可将两个基因型文件整合成一个。

## 功能

1. 将多个**HMP**文件或其他标准基因型合成为一个文件，并对齐标记名。
2. 文件的标记必须有相同的部分。例如A文件标记名：A1,A2,A3,A5；B文件标记名：A2,A3,A4,A5。
3. 可以选择是否将重复的材料删掉（保留第一个）。
4. 可以选择是否生成TASSEL能打开的HMP文件，只针对HMP文件。

## 使用方法

将下列代码复制到Rstudio，请按需求修改下面第一行，然后全选运行。期间会弹出对话框，选择2个要合并的基因型文件即可。注意，基因型文件需要是hmp格式文件。

合并后的文件命名为【combineGenoData.txt】或【combineGenoData.hmp.txt】。

海外：

```R
HMP <- TRUE   #? TRUE or FALSE, TRUE is generating a usable format for TASSEL
delRep <- TRUE  #? TRUE or FALSE, TRUE is Auto-delete duplicates
# Please choose two marker files once
source("https://aozhangchina.github.io/R/MerkerDataCombineTool/MerkerDataCombineTool.r")   # 加载程序文件，需要联网
```

国内：

```R
HMP <- TRUE   #? TRUE or FALSE, TRUE is generating a usable format for TASSEL
delRep <- TRUE  #? TRUE or FALSE, TRUE is Auto-delete duplicates
# Please choose two marker files once
source("https://dataholdcn.cn/R/MerkerDataCombineTool/MerkerDataCombineTool.r")   # 加载程序文件，需要联网
```

