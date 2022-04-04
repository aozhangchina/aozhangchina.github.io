if(!exists("fName")){fName <- "newMyGeno"}
if(!exists("subtext")){subtext <- ":.*"}
if(!exists("missingSign")){missingSign <- "N"}
if(!exists("chooseOne")){chooseOne <- TRUE}
cat("Updated: 20220331 10:50:27\n")
tempC1 <- NULL
source("https://dataholdcn.cn/R/GBIT/GBIT.R")
integratingMarkers <- function(myVector){
    myVector[which(myVector == missingSign)] <- NA
  if (any(is.na(myVector))){
    myVector <- myVector[-which(is.na(myVector))]
  }

  table_vector <- table(myVector)

  if(length(table(table_vector)) > 0){
    myVector <- names(sort(table_vector,decreasing=T))[1]
  }else{
    myVector <- "N"
  }
  return(myVector)
}
GBIT.setwd(choose.dir())
myGeno <- GBIT.readFile(choose= T,header = T)

myGenoHeaderClean <- GBIT.geno.simplifiedName(myGeno,subtext = subtext)

myGenoHeaderTwo <- cbind(names(myGeno),names(myGenoHeaderClean))
write.table(myGenoHeaderTwo,paste0("GBIT/",fName,".name.compare.csv"),sep=",",row.names=FALSE,col.names = FALSE,quote = FALSE)

myGenoHeaderRep <- unique(GBIT.geno.delDup(myGenoHeaderClean))

if (length(myGenoHeaderRep)>0){
  names(myGeno) <- names(myGenoHeaderClean)

vlist <- as.data.frame(cbind(myGenoHeaderRep,c(1:length(myGenoHeaderRep))))
names(vlist) <- c("Names","GID")

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
  a <- cbind(myGeno[,which(names(myGeno)==local1)],myGeno[,which(names(myGeno)==local2)])
  b <- apply(a,1,same)
  bc <- apply(a,1,noNA)
  
  res <- c(res,length(which(b==TRUE))/length(which(bc==TRUE)))
}
names(res) <- myGenoHeaderRep
write.table(res,"clipboard",sep="\t",col.names=FALSE)
cat("The two taxa accuracy has been exported in the clipboard.\n")

if (chooseOne==TRUE){
  index <- which(duplicated(names(myGenoHeaderClean)))
  if (length(index)>0){
  myGeno2 <- myGeno[,-index]
  }
}else{
  index <- which(names(myGenoHeaderClean) %in% myGenoHeaderRep)
  myGeno2 <- myGeno[,-index]
  for (i in myGenoHeaderRep){
      cat("Calculating ",i,"\n")
  tempA <- which(names(myGeno) == i)
  tempB <- as.matrix(cbind(myGeno[tempA]))
  tempB <- as.data.frame(tempB)
  NA_vector <- NULL

  system.time(tempC<- apply(tempB,1,integratingMarkers))
  tempC1 <- cbind(tempC1,tempC)
  }
  colnames(tempC1) <- myGenoHeaderRep
  
  myGeno2 <- cbind(myGeno2,tempC1)
}

  GBIT.writeHMP(myGeno2,fName)
}else{
  cat("No duplicate material was found.\n")
}