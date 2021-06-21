source("https://dataholdcn.cn/R/GBIT/GBIT.R")
if (!exists('HMP')) {HMP <- TRUE} #? TRUE or FALSE, TRUE is generating a usable format for TASSEL
if (!exists('delRep')) {delRep <- TRUE} #? TRUE or FALSE, TRUE is Auto-delete duplicates
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

cat("This tool is made by Zhang Ao. Update: 2020-02-25. \nhttps://datahold.cn/\n")
cat("The number of individuals are:\n")
for (j in 1: length(fileInfo$f.name)){
  ntext <- paste0("print(length(n",j," <- names(geno", j,")))\n")
  eval(parse(text=ntext))
}

nintersect <- n1   # initialize ninterect
for (k in 2: length(fileInfo$f.name)){
  nintersecttext <- paste0("nintersect <- intersect(n",k,",nintersect)")
  eval(parse(text=nintersecttext))
}

if (length(nintersect)==11){
  cat("Pass!\n",nintersect)
}else{
  cat("Warning!\n",nintersect)
}

n2intersect <- geno1[,1]
for (k2 in 2: length(fileInfo$f.name)){
  n2intersecttext <- paste0("n2intersect <- intersect(geno",k2,"[,1],n2intersect)")
  eval(parse(text=n2intersecttext))
}

for (k3 in 1:length(fileInfo$f.name)){
  txt <- paste0("index <- geno",k3,"[,1] %in% n2intersect")
  eval(parse(text=txt))
  txt <- paste0("geno",k3," <- geno",k3,"[index,]")
  eval(parse(text=txt))
}

cat("The number of final markers: \n",length(n2intersect),"\n")

geno <- geno1 # Initialization
# combine all materials
for (l in 2:length(fileInfo$f.name)){
  txt <- paste0("index <- names(geno",l,") %in% nintersect")
  eval(parse(text=txt))
  txt <- paste0("geno",l," <- geno",l,"[,!index]")
  eval(parse(text=txt))
  combinetext <- paste0("geno <- cbind(geno,geno",l,")")
  eval(parse(text=combinetext))
  txt <- paste0("rm(geno",l,")")
  eval(parse(text=txt))
}
rm(geno1)
if (length(which(geno[,1]==""))>0){
geno <- geno[-which(geno[,1]==""),]
}

if (length(unique(geno[,1]))<length(geno[,1])){
  cat("There are some repeat markers:\n")
  repGeno <- which(duplicated(geno[,1]))
  print(geno[repGeno,1])
}

if (length(unique(names(geno)))<length(names(geno))){
  cat("There are some repeat material:\n")
  repGeno <- which(duplicated(names(geno)))
  print(names(geno)[repGeno])
}

if (delRep==TRUE & length(which(duplicated(names(geno))))){
  geno <- geno[,-which(duplicated(names(geno)))]
}else{
  cat("You can set 'delRep <- TRUE' to auto-delete duplicates.\n")
}

if (HMP==TRUE){
  GBIT.writeHMP(geno,"combineGenoData")
}else{
  write.table(geno,"combineGenoData.txt",sep="\t",row.names=F,quote=F)
}

