if (!exists('output_f')) {output_f = F} #? T 输出结果，F 不输出
source("https://dataholdcn.cn/R/GBIT/GBIT.R")
cat("Work space.\n")
GBIT.setwd(choose.dir())

cat("Pheno file.\n")
GCASCAFile <- GBIT.readFile(header = T, choose = T)

rownames(GCASCAFile) <- demoGCA_one_rows <- GCASCAFile[,1]
GCASCAFile <- GCASCAFile[,-1]

GCA_SCA <- GBIT.get.GCASCA(GCASCAFile,output_f = F)