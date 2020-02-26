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
  chromAndPos[chromAndPosNA] <- "NA"   # 将NA转化�?"NA"
  chromAndPos <- as.data.frame(chromAndPos)
  temp <- wData[order(chromAndPos$chrom,chromAndPos$pos),]   # 排列非数字内�?
  sortResult <- temp[order(as.numeric(as.character(temp$chrom)),as.numeric(as.character(temp$pos))),] 
  write.table(sortResult,paste0(fName,".hmp.txt"),sep="\t", quote =F,row.names = F)   # 输出hmp文件
  print("directory:")
  print(getwd())
}

fileInfo <- getFileInfo()
setwd(fileInfo$f.dir[1])
cat(paste0("Work Space:\n",getwd(),"\n"))

# read files
for (i in 1:length(fileInfo$f.name)){
  genotext <- paste0("geno",i," <- readFiles(header=T,fname='",fileInfo$f.name[i],"')")
  eval(parse(text=genotext))
  genotext <- paste0("geno",i,"<- as.data.frame(geno",i,")")
  eval(parse(text=genotext))
}

cat("This tool coding by Zhang Ao. Update: 2020-02-25. \nhttps://datahold.cn/\n")
cat("The number of markers are:\n")
for (j in 1: length(fileInfo$f.name)){
  ntext <- paste0("print(length(n",j," <- names(geno", j,")))\n")
  eval(parse(text=ntext))
}

nintersect <- n1   # initialize ninterect
for (k in 2: length(fileInfo$f.name)){
  nintersecttext <- paste0("nintersect <- intersect(n",k,",nintersect)")
  eval(parse(text=nintersecttext))
}

cat("The number of final markers: \n")
cat(paste(length(nintersect),"\n"))

for (kk in 1:length(fileInfo$f.name)){
  indexMtext <- paste0("index",kk," <- n",kk," %in% nintersect")
  eval(parse(text=indexMtext))
  genoMtext <- paste0("genoM",kk," <- geno",kk,"[,index",kk,"]")
  eval(parse(text=genoMtext))
}


geno <- "" # Initialization
# combine all materials
for (l in 1:length(fileInfo$f.name)){
  combinetext <- paste0("geno <- rbind(genoM",l,",geno)")
  eval(parse(text=combinetext))
}
if (length(which(geno[,1]==""))>0){
geno <- geno[-which(geno[,1]==""),]
}

if (length(unique(geno[,1]))<length(geno[,1])){
  cat("There are some repeat materials:\n")
  repGeno <- which(duplicated(geno[,1]))
  print(geno[repGeno,1])
}

if (delRep==TRUE & length(which(duplicated(geno[,1])))){
  geno <- geno[-which(duplicated(geno[,1])),]
}else{
  cat("You can set 'delRep <- TRUE' to auto-delete duplicates.\n")
}

if (HMP==TRUE){
  writeHMP(geno,"combineGenoData")
}else{
  write.table(geno,"combineGenoData.txt",sep="\t",row.names=F,quote=F)
}
