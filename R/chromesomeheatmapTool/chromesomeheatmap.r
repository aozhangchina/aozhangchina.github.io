source("https://dataholdcn.cn/R/GBIT/GBIT.R")
if (exists("chr_size")){chr_size <- "value"}
if(exists("worksp")==F){
  worksp_f <- choose.files()
  setwd(dirname(worksp_f))
}

GBIT.library("readr")
info_chr <- choose.files()
if (length(info_chr)>0){
  info_chr_f <- read.delim(info_chr)
}

MyGeno <- read.delim(worksp_f,nrows =1)
MyGenoName <- names(MyGeno[3:4])
MyGeno <-  read_delim(worksp_f,delim="\t")
MyGeno <- MyGeno[,MyGenoName]


myDf <- data.frame(Chr=NA,Start=0,End=c(307041717,
                                        244442276,
                                        235667834,
                                        246994605,
                                        223902240,
                                        174033170,
                                        182381542,
                                        181122637,
                                        159769782,
                                        150982314),CE_start=c(127050000,+
                                                                96900000,+
                                                                88400000,+
                                                                82400000,+
                                                                106600000,+
                                                                32250000,+
                                                                49500000,+
                                                                68250000,+
                                                                34150000,+
                                                                41550000))
myDf$CE_end <- myDf$CE_start + 5000000

if(length(which(is.na(MyGeno$chrom)==TRUE))>0){
  MyGeno <- MyGeno[-which(is.na(MyGeno$chrom)),]
}

MyGeno$chrom <- as.character(MyGeno$chrom)

if (identical(unique(MyGeno$chrom),as.character(1:10))){
  myDf$Chr <- unique(MyGeno$chrom)
}

if (length(info_chr)>0){
  for (i in 2:5){
    if (!any(is.na(info_chr_f[,i]))){
      for (j in 1:nrow(myDf)){
        myDf[which(myDf$Chr==info_chr_f[j,1]),i] <- info_chr_f[j,i]
      }
    }
    
  }
}

for (i in 1:nrow(myDf)){
  if (myDf$End[i] - max(MyGeno$pos[which(MyGeno$chrom==i)]) > 0){
    if (myDf$End[i] - max(MyGeno$pos[which(MyGeno$chrom==i)]) < 1000000){
      myDf$End[i] <- myDf$End[i] + 1000000
    }
  }else{
    cat("The chr ",i," position size is not enough.\n")
  }
}


myDensity <- NULL

for (i in 1:nrow(myDf)){
  temp <- subset(MyGeno,chrom==i)
  if (myDf$End[i] < max(temp$pos)){
    myDf$End[i] <- max(temp$pos)
  }
  breaknumber <- seq(1,((ceiling(max(temp$pos)/1000000)+1)*1000000),1000000)
  tempdata <- hist(temp$pos, breaks = breaknumber)
  getHistInfo <- data.frame(Chr=i,Start=breaknumber[-length(breaknumber)],End=(tempdata$breaks[-1]-1),Value=tempdata$counts)
  if (is.null(myDensity)){
    myDensity <- getHistInfo
  }else{
    myDensity <- rbind(myDensity,getHistInfo)
  }
}

if (chr_size == "auto"){
  for (i in 1:nrow(myDf)){
    temp <- subset(MyGeno,chrom==i)
    myDf$End[i] <- max(temp$pos) + 1000000
  }
}

labelmyDensity <- cbind(myDensity,Color="fc8d62")

GBIT.library("RIdeogram")

ideogram(karyotype = myDf, overlaid = myDensity, colorset1 = plotcolor)
convertSVG("chromosome.svg", device="png")
convertSVG("chromosome.svg", device = "tiff", dpi = 600)
