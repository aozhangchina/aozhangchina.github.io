source("https://dataholdcn.cn/R/GBIT/GBIT.R")
GBIT.WD(NULL)
mydata <- GBIT.readFile(choose = TRUE)
Ns <- names(mydata)
Tn <- GBIT.grep(Ns,"T\\d+","Time")
Vn <- GBIT.grep(Ns,"V\\d+","value")
if (Tn$lst-Tn$fst != Vn$lst-Vn$fst){
  stop(paste0("The number of time is not equal to veriable.\n"))
}
getmode <- function(x){
  xt <- NULL
  xt <- as.numeric(names(sort(table(x),decreasing = TRUE)))[1]
  return(xt)
}
getmean <- function(x){
  xt <- NULL
  xt <- mean(x,na.rm=TRUE)
  return(xt)
}
TnTable <- apply(mydata[Tn$fst:Tn$lst], 2, getmode)
VnTable <- apply(mydata[Vn$fst:Vn$lst], 2, getmean)


getMeanTable<- function(x,y){
  if (any(is.na(x))){
    x <- x[-which(is.na(x))]
  }
  if (any(is.na(y))){
    y <- y[-which(is.na(y))]
  }
  NnTable <- matrix(NA,ncol=(length(x)-1),nrow=2)
  for (i2 in 1:(length(x)-1)){
    NnTable[1,i2] <- as.numeric((x[i2+1]+x[i2])/2)
    NnTable[2,i2] <- as.numeric((y[i2+1]+y[i2])/2)
  }
  return(NnTable)
}

a2 <- getMeanTable(TnTable,VnTable)

AUDPC <- function(x){
  a <- 0
  for(i in 1:ncol(x)){
    a <- a + x[1,i]*x[2,i]
  }
  return(a)
}

myAUDPC <- AUDPC(a2)

# plot(TnTable,VnTable,xlim=c(0,max(TnTable)),ylim=c(0,max(VnTable)),type = "o")
plot(TnTable,VnTable,type = "o")
text((TnTable[1]/2)+(TnTable[2]/2),a2[2,ncol(a2)],c("AUDPC\n\n",myAUDPC))

for(i in 1:ncol(a2)){
  segments(x0=TnTable[i],y0=VnTable[1],x1=TnTable[i+1],y1=VnTable[1],col="red",lty=6)
  segments(x0=TnTable[i],y0=a2[2,i],x1=TnTable[i+1],y1=a2[2,i],col="red",lty=6)
  segments(x0=TnTable[i],y0=VnTable[1],x1=TnTable[i],y1=a2[2,i],col="red",lty=6)
  segments(x0=TnTable[i+1],y0=VnTable[1],x1=TnTable[i+1],y1=a2[2,i],col="red",lty=6)
}

eachT <- mydata[,Tn$fst:Tn$lst]
eachV <- mydata[,Vn$fst:Vn$lst]
myAUDPCeach <- rep(NA,nrow(eachT))
for (i in 1:nrow(eachT)){
  a2each <- getMeanTable(eachT[i,],eachV[i,])
  myAUDPCeach[i] <- AUDPC(a2each)
}
names(myAUDPCeach) <- mydata$Taxa
myAUDPCeach[nrow(eachT)+1] <- myAUDPC
names(myAUDPCeach)[nrow(eachT)+1] <- "mean"
write.table(myAUDPCeach,paste0("GBIT\\AUDPCs.txt"),col.names = F,quote = F)
