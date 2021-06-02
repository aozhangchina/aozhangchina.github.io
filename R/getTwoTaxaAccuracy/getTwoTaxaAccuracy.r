fGeno <- GBIT.readFile(header=TRUE,choose=TRUE) # HMP
vlist <- read.csv(choose.files())
# *.csv
# Names	GID
# GMPGBS24-P12:MRG:3:250142093,275823
# GMPGBS24-P13:MRG:3:250142210,275823
# GMPGBS25-P12:MRG:3:250142094,34867
# GMPGBS25-P13:MRG:3:250142209,34867
# GMPGBS541:MRG:3:250142096,279158
# GMPGBS541:MRG:3:250142163,279158

same <- function(ss){
  if (ss[1]==ss[2]&(ss[1]!="N"|ss[2]!="N")){
    ind <- TRUE
  }else{
    ind <- FALSE
  }
}
noNA <- function(ss){
  if(length(which(ss=="N"))==0){
    ind <- TRUE
  }else{
    ind <- FALSE
  }
}

res <- numeric()
for (i in 1:length(unique(vlist$GID))){
  local1 <- vlist[which(vlist$GID==unique(vlist$GID)[i])[1],1]
  local2 <- vlist[which(vlist$GID==unique(vlist$GID)[i])[2],1]
  a <- cbind(fGeno[,which(names(fGeno)==local1)],fGeno[,which(names(fGeno)==local2)])
  b <- apply(a,1,same)
  bc <- apply(a,1,noNA)
  
  res <- c(res,length(which(b==TRUE))/length(which(bc==TRUE)))
}
names(res) <- unique(vlist$GID)
write.table(res,"clipboard",sep="\t",col.names=FALSE)
