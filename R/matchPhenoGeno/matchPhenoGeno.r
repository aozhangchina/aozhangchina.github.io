source("https://dataholdcn.cn/R/GBIT/GBIT.R")
if(!exists("fName.geno")){fName.geno <- "myGeno"}
if(!exists("fName.pheno")){fName.pheno <- "myPheno"}
if(!exists("ID")){ID <- "GenoID"}
if(!exists("subtext")){subtext <- ":.*|-.*"}

GBIT.setwd(choose.dir())
cat("Geno\n")
myGeno.f <- choose.files()
cat("Pheno\n")
myPheno.f <- choose.files()
myGeno <- GBIT.readFile(choose = F,header = T,fname = myGeno.f)
myPheno <- GBIT.readFile(choose = F,header = T,fname = myPheno.f)

myGeno_head <- myGeno[,1:11]

myPhenoClean <- GBIT.pheno.simplifiedName(myPheno,which(names(myPheno)==ID),subtext = subtext)

myGenoHeaderClean <- GBIT.geno.simplifiedName(myGeno,subtext = subtext)

all_material <- intersect(names(myGenoHeaderClean),myPhenoClean$GenoID)

cat(length(all_material)," individuals were matched.\n")
cat("Geno: ",setdiff(names(myGenoHeaderClean)[-c(1:11)],all_material),".\n")
cat("Pheno: ",setdiff(myPhenoClean$GenoID,all_material),".\n")

myGeno <- cbind(myGeno_head,myGenoHeaderClean[,all_material])
myPheno <- myPhenoClean[-which(myPhenoClean$GenoID %in% setdiff(myPhenoClean$GenoID,all_material)),]

write.table(myPheno,paste0(fName.pheno,".csv"),sep = ",",quote = F,col.names = T,row.names = F)
GBIT.writeHMP(wData = myGeno,fName = paste0(fName.geno,"_",(ncol(myGeno)-11),"_",nrow(myGeno)))
