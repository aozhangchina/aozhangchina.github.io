source("https://dataholdcn.cn/R/GBIT/GBIT.R")
cat("Work space. A folder.\n")
GBIT.setwd(choose.dir())
cat("Genotype. (*.hmp.txt)\n")
geno.info <- choose.files()
cat("Phenotype.\n")
pheno.info <- choose.files()
cat("A list of SNP.\n")
slist.info <- choose.files()

geno <- GBIT.readFile(header=T,fname = geno.info)
pheno <- GBIT.readFile(header=T,fname = pheno.info)
slist <- GBIT.readFile(header=F,fname = slist.info)

GBIT.geno.delDup(geno)

geno_name <- names(geno[12:ncol(geno)])
pheno_name <- pheno[,1]

ck <- GBIT.matching.test(geno_name,pheno_name)

if (!ck){
  geno_1 <- geno[,1:11]
  geno_2 <- geno[,12:ncol(geno)]
  p_i <- which(names(pheno) %in% c("Genotype","GenoID","ID","Taxa"))
  pheno_2 <- pheno[,p_i]
  index <- intersect(names(geno_2),pheno_2)
  geno_2 <- geno_2[,index]
  pheno <- pheno[which(pheno[,p_i] %in% index),]
  pheno <- pheno[order(pheno[,p_i]),]
  geno_2 <- geno_2[,order(names(geno_2))]
  geno  <- cbind(geno_1,geno_2)
}

myPVE <- GBIT.get.PVE(geno,pheno,slist)
cat("The results of PVEs are in your clipboard.\n")
write.table(myPVE,"GBIT\\PVE.txt",col.names=FALSE,sep="\t",quote=FALSE)
write.table(myPVE,"clipboard",col.names=FALSE,sep="\t")
