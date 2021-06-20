if(!exists("fName")){fName <- "newMyGeno"}
source("https://dataholdcn.cn/R/GBIT/GBIT.R")
GBIT.setwd(choose.dir())
myGeno <- GBIT.readFile(choose= T,header = T)

myGenoHeaderClean <- GBIT.geno.simplifiedName(myGeno,subtext = ":.*|-.*")

myGenoHeaderTwo <- cbind(names(myGeno),names(myGenoHeaderClean))
write.csv(myGenoHeaderTwo,paste0("GBIT/",fName,".name.compare.csv"),row.names=FALSE,col.names = NULL)

myGenoHeaderRep <- GBIT.geno.delDup(myGenoHeaderClean)

index <- which(duplicated(names(myGenoHeaderClean)))
myGeno <- myGeno[,-index]

GBIT.writeHMP(myGeno,fName)