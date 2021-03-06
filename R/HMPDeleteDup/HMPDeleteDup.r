if(!exists("fName")){fName <- "newMyGeno"}
if(!exists("subtext")){subtext <- ":.*"}
source("https://dataholdcn.cn/R/GBIT/GBIT.R")
GBIT.setwd(choose.dir())
myGeno <- GBIT.readFile(choose= T,header = T)

myGenoHeaderClean <- GBIT.geno.simplifiedName(myGeno,subtext = subtext)

myGenoHeaderTwo <- cbind(names(myGeno),names(myGenoHeaderClean))
write.table(myGenoHeaderTwo,paste0("GBIT/",fName,".name.compare.csv"),sep=",",row.names=FALSE,col.names = FALSE,quote = FALSE)

myGenoHeaderRep <- GBIT.geno.delDup(myGenoHeaderClean)

names(myGeno) <- names(myGenoHeaderClean)

index <- which(duplicated(names(myGenoHeaderClean)))
if (length(index)>0){
  myGeno2 <- myGeno[,-index]
}

GBIT.writeHMP(myGeno,fName)
