#################### 头部信息开始 ####################
# 功能：提琴图制作工具
# 制作人：张敖（Ao Zhang）
# 更新日期： 2019-10-10 
# Violin Plot tool
# Made by Ao Zhang
# Updata:10/10/2019
#################### 头部信息结束 ####################

# install and load vioplot  package
if(!"vioplot" %in% installed.packages()) {
  install.packages("vioplot")
  library(vioplot)
}else{
  library(vioplot)
}

# transfor vector to string
stringAddQuotation <- function(theVector=theVector,sep=",",sepIn="'",pre_element="",bac_element="",first="'",last="'"){
  # 前、中、后，解向量
  theVector <- paste0(pre_element,theVector,bac_element,collapse = sep)
  # 内部,变成','，基于first
  sep_in <- paste0(sepIn,sep,sepIn)
  # 分割符号前，增加引号
  theVector <- gsub(sep,sep_in,theVector)
  # 前后增加引号
  theVector <- paste0(first,theVector,last)
  return(theVector)
}

cat("请选择用来画图的数据文件！\nPlease choose the file for making plot.")
datainfo <- file.choose()
usefiles <- read.table(datainfo,header = T, sep="\t")
str(usefiles)
setwd(dirname(datainfo))
xnames <- names(usefiles)

xnames1 <- stringAddQuotation(xnames,",","","usefiles$","","","")
xnames2 <- stringAddQuotation(xnames,",","'","","","c('","')")


runfile <- parse(text=paste0("vioplot(",xnames1,
                               ", xlab = '",Xlab,
                               "', ylab='",Ylab,
                               "', names=",xnames2,
                               ",col='",plotcolor,
                               "')"))

eval(runfile)
