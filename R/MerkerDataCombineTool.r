# Read files
readFiles <- function(header=TRUE,choose=FALSE,fname){
  if(!"readr" %in% installed.packages()) {
    install.packages("readr")
    library(readr)
  }else{
    library(readr)
  }
  if (choose==TRUE){
    myFileName <- file.choose()
  }else{
    myFileName <- fname
  }
  if (grepl("\\.csv$",myFileName)){
    if (header==TRUE){
      myFile <- read_csv(myFileName, col_names = TRUE)
    }else{
      myFile <- read_csv(myFileName, col_names = FALSE)
    }
  }else{
    if (header==TRUE){
      myFile <- read_tsv(myFileName, col_names = TRUE)
    }else{
      myFile <- read_tsv(myFileName, col_names = FALSE)
    }
  }
  return(myFile)
}

# 获得文件信息 
getFileInfo <- function(){
myFileName <- choose.files()
f.name <- basename(myFileName)
f.dir <- dirname(myFileName)
f <- list(f.name=f.name,f.dir=f.dir)
return(f)
}

# 将Tassel不能打开的HMP文件转换成Tassel能够打开的HMP文件
writeHMP <- function(wData,fName){
  chromAndPos <- data.frame(chrom=wData$chrom,pos=wData$pos)
  chromAndPos <- as.matrix(chromAndPos)
  chromAndPosNA <- which(is.na(chromAndPos))
  chromAndPos[chromAndPosNA] <- "NA"   # 将NA转化为"NA"
  chromAndPos <- as.data.frame(chromAndPos)
  temp <- wData[order(chromAndPos$chrom,chromAndPos$pos),]   # 排列非数字内容
  sortResult <- temp[order(as.numeric(as.character(temp$chrom)),as.numeric(as.character(temp$pos))),] 
  write.table(sortResult,paste0(fName,".hmp.txt"),sep="\t", quote =F,row.names = F)   # 输出hmp文件
  print("directory:")
  print(getwd())
}


fileInfo <- getFileInfo()
setwd(fileInfo$f.dir[1])
cat(paste0("Work Space:\n",getwd()))

for (i in 1:length(fileInfo$f.name)){
  genotext <- paste0("geno",i," <- readFiles(header=T,fname='",fileInfo$f.name[i],"')")
  eval(parse(text=genotext))
  genotext <- paste0("geno",i,"<- as.data.frame(geno",i,")")
  eval(parse(text=genotext))
}


n1 <- as.character(geno1[,1])
length(n1)
n2 <- as.character(geno2[,1])
length(n2)
nintersect <- intersect(n1,n2)
length(nintersect)

index1 <- n1 %in% nintersect
index2 <- n2 %in% nintersect

geno1 <- geno1[index1,]
geno2 <- geno2[index2,]

geno <- cbind(geno1,geno2)

if (length(unique(names(geno)))<length(names(geno))){
  geno <- geno[,!duplicated(names(geno))]
}


if (HMP==TRUE){
  writeHMP(geno,"combineGenoData")
}else{
  write.table(geno,"combineGenoData.txt",sep="\t",row.names=F,quote=F)
}
