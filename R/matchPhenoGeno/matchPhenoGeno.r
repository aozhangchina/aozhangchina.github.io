source("https://dataholdcn.cn/R/GBIT/GBIT.R")
if(!exists("fName.geno")){fName.geno <- "myGeno"}
if(!exists("fName.pheno")){fName.Pheno <- "myPheno"}
if(!exists("ID")){ID <- "GenoID"}
if(!exists("subtext")){subtext <- ":.*|-.*"}

GBIT.setwd(choose.dir())
myGeno <- GBIT.readFile(choose = T,header = T)
myPheno <- GBIT.readFile(choose = T,header = T)

myGeno_head <- myGeno[,1:11]

myPhenoClean <- GBIT.pheno.simplifiedName(myPheno,which(names(myPheno)==ID),subtext = subtext)

myGenoHeaderClean <- GBIT.geno.simplifiedName(myGeno,subtext = subtext)

all_material <- intersect(names(myGenoHeaderClean),myPhenoClean$GeneID)

cat(length(all_material)," individuals were matched.\n")
cat("Geno: ",setdiff(names(myGenoHeaderClean)[-c(1:11)],all_material),".\n")
cat("Pheno: ",setdiff(myPhenoClean$GeneID,all_material),".\n")

myGeno <- cbind(myGeno_head,myGenoHeaderClean[,all_material])
myPheno <- myPhenoClean[-which(myPhenoClean$GeneID %in% setdiff(myPhenoClean$GeneID,all_material)),]

write.table(myPheno,paste0(fName.Pheno,".csv"),sep = ",",quote = F,col.names = T,row.names = F)
GBIT.writeHMP(wData = myGeno,fName = paste0(fName.geno,"_",ncol(myGeno),"_",nrow(myGeno)))
