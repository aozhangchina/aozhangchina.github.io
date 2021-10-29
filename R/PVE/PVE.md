# PVE 

[TOC]

**制作：张敖 [https://datahold.cn](https://datahold.cn/)** 

# 鸣谢

张学才，CIMMYT，指导

刘新春，四川农业大学，测试

## 最新更新

【2021-10-29 10:54:37】

修正方法描述不正确的问题。

【2021-08-05 23:43:53】

修正因函数名问题导致的错误。

【2021-06-02 13:58:37】

功能完成。

## 功能

1. 根据GWAS获得的超过阈值的显著标记计算PVE（Phenotypic variance explained）=R<sup>2</sup>。

## 使用方法

1. 计算全过程，请保证网络通畅。

2. 准备表型文件，两列，包含标题，第一列是标记名，第二列是表型值。【[示例文件.txt](https://dataholdcn.cn/R/PVE/pheno_PVE.txt)】

   Trait	CombinedENV
   CML161	0.405449507
   CML165	0.190992175
   CML206	0.438263016
   CML247	0.261138727
   CML254	0.330294057
   CML287	0.077330646
   CML300	0.098037893
   CML301	0.157140089
   CML302	0.213946167

3. 准备显著SNP标记列表，每行一个，不含标题。【[示例文件.txt](https://dataholdcn.cn/R/PVE/slist.txt)】

   S6_145529216	
   S3_140023577
   S5_86861333
   S1_283640778

4. 将下列代码复制到Rstudio，然后全选运行。期间会弹出对话框3次，第一次选择工作目录，第二次选择基因型文件（*.hmp.txt），第三次选择表型文件。

   表型文件格式：

   | Taxa                         | CombinedENV |
   | ---------------------------- | ----------- |
   | RCGS_1:C3KBGACXX:6:250279881 | 0.598732    |
   | RCGS_2:C3KBGACXX:6:250279882 | 0.117541    |
   | RCGS_4:C3KBGACXX:6:250279884 | 0.255027    |
   | RCGS_5:C3KBGACXX:6:250279885 | 0.449318    |
   | RCGS_9:C3KBGACXX:6:250279889 | 0.277038    |

海外：

```R
source("https://aozhangchina.github.io/R/PVE/PVE.r")   # 加载程序文件，需要联网
```

国内：

```R
source("https://dataholdcn.cn/R/PVE/PVE.r")   # 加载程序文件，需要联网
```

