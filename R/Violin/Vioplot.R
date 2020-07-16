#################### 头部信息开始 ####################
# 功能：提琴图制作工具
# 制作人：张敖（Ao Zhang）
# 更新日期： 2020-07-16
# Violin Plot tool
# Made by Ao Zhang
# Updata:7/16/2020
#################### 头部信息结束 ####################
if (!exists("Xlab")) {
   Xlab <- "Xlab"
}
if (!exists("Ylab")) {
   Ylab <- "Ylab"
}
if (!exists("plotcolor")) {
   tempPlotcolor <- paste0("'gray'")
}else if (length(plotcolor)>=1) {
  tempPlotcolor <- "c("
   for (plotcolor_n in 1:length(plotcolor)) {
      tempPlotcolor <- paste0(tempPlotcolor,"'",plotcolor[plotcolor_n],"',")
   }
   tempPlotcolor <- sub(",$","",tempPlotcolor)
   tempPlotcolor <- paste0(tempPlotcolor,")")
}
if (!exists("Ylim")) {
   Ylim <- NULL
}


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
  # 内部,变成‘，
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

for (i in 1:ncol(usefiles)){
  usefiles[,i] <- as.numeric(usefiles[,i])
}

xnames1 <- stringAddQuotation(xnames,",","","usefiles$","","","")
xnames2 <- stringAddQuotation(xnames,",","'","","","c('","')")

if(length(Ylim)==0){
  runfile <- parse(text=paste0("vioplot(",xnames1,
                               ", xlab = '",Xlab,
                               "', ylab='",Ylab,
                               "', names=",xnames2,
                                ",col=",tempPlotcolor,
                               ")"))
}else{
  runfile <- parse(text=paste0("vioplot(",xnames1,
                               ", xlab = '",Xlab,
                               "', ylab='",Ylab,
                               "', names=",xnames2,
                                ",col=",tempPlotcolor,
                                ",ylim=c(",Ylim[1],",",Ylim[2],
                               "))"))
}
eval(runfile)
