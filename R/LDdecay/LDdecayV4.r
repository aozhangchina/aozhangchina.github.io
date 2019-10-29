system.time({
cat("LD Decay Plot Tool by Zhang Ao\t LD衰减图工具由张敖制作\nMost recent update\t最后更新 2019.10.29\n-------\nPlease choose the data of LDdata\t 请选择连锁不平衡文件")
myFile <- file.choose()
file.name <- basename(myFile)
file.dir <- dirname(myFile)
setwd(file.dir)
if(!"readr" %in% installed.packages()) {
  install.packages("readr")
  library(readr)
}else{
  library(readr)
}

temp <- read_csv(file.name,col_name=TRUE)

colnames(temp)[which(colnames(temp)=="R^2")] <- "R.2"
myLD <- data.frame(locus=temp$Locus1,Dist_bp=(temp$Dist_bp),R.2=temp$R.2)   
#rm(temp)

coordinateScale <- sort(coordinateScale)
ac <- c(coordinateScale[1],coordinateScale[length(coordinateScale)],coordinateScale[c(-1,-length(coordinateScale))])

chrdataframe <- data.frame(chrlocs=unique(myLD$locus),
                           chrs=paste0("Chr",unique(myLD$locus)),
                           locs=c(1:length(unique(myLD$locus))))

chrs <- as.character(chrdataframe$chrs)
locs <- as.character(chrdataframe$chrlocs)


if (length(lineColor)==1){
  lineColor <- c(lineColor,rainbow(length(chrs)))
}

myScale <- as.data.frame(matrix(rep(NA,length(ac)*(length(chrs)+2)),length(ac),(length(chrs)+2)))
names(myScale) <- c("Scale",chrs,"Mean")
myScale$Scale <- coordinateScale
# The value of scale 0 is 1.
myLD[,1] <- as.character(myLD[,1])
myLD[,2] <- as.numeric(myLD[,2])
myLD[,3] <- as.numeric(myLD[,3])
str(myLD)
for (i in 1:dim(myScale)[1]){
  cat(paste("Calculating Scale ",myScale$Scale[i], "\n",sep=''))
  if (myScale$Scale[i]==0){
    myScale[which(myScale$Scale==0),c(2:dim(myScale)[2])] <- round(1.00,2)
  }else if(i==2){
    for (j in 1:length(chrs)){
      txtcode <- paste0("myScale$",chrs[j], "[2] <- mean(myLD$R.2[which(myLD$Dist_bp<=myScale$Scale[2]&myLD$locus=='",locs[j],"')],na.rm=T)")
      eval(parse(text=txtcode))
    }
    myScale$Mean[2] <- mean(myLD$R.2[which(myLD$Dist_bp<=myScale$Scale[2])],na.rm=T)
  }else{
    for (j in 1:length(chrs)){
      txtcode <- paste0("myScale$",chrs[j],"[",i,"]" ,"<-mean(myLD$R.2[which(myLD$Dist_bp<=myScale$Scale[",i,"]&myLD$Dist_bp>myScale$Scale[",i-1,"]&myLD$locus=='",locs[j],"')],na.rm=T)")
      eval(parse(text=txtcode))
    }
    myScale$Mean[i] <-mean(myLD$R.2[which(myLD$Dist_bp<=myScale$Scale[i]&myLD$Dist_bp>myScale$Scale[i-1])],na.rm=T)
  }
}  

if (unit=="Kb"){
  myScale$Scale <- coordinateScale/1000
}else{
  myScale$Scale <- coordinateScale
}
print(round(myScale,2))


op<-par(no.readonly=TRUE)
par(op)

if (unit=="Kb"){
  plot(myScale$Scale,myScale$Mean,xaxt="n",yaxt="n",type=dot,lwd=lineWidth, col=lineColor,main=paste(plotTitle,sep=' '), ylim = c(0,1) ,xlab=expression("Physical Distance (Kb)"), ylab=expression("r"^" 2"))
}else{
  plot(myScale$Scale,myScale$Mean,xaxt="n",yaxt="n",type=dot,lwd=lineWidth, col=lineColor,main=paste(plotTitle,sep=' '), ylim = c(0,1) ,xlab=expression("Physical Distance (bp)"), ylab=expression("r"^" 2"))  
}

#  acc <- ac[-c(which(ac>0&ac<5000))]
ad <- seq(1000,4000,1000)
acc <- c(ac[1:2],ac[6:9])

if(unit=="Kb"){
  axis(side=1,at=acc/1000)
  axis(side=1,at=ad/1000,labels=FALSE)
}else{
  axis(side=1,at=acc)
  axis(side=1,at=ad,labels=FALSE)
}

axis(side=2,at=seq(0.0,1.0,0.1))
if (dottedLine == TRUE){
  abline(h=dividingLine,lty=2,col="red")
}
mine<-function(x1,x2,y1,y2,y){
  a <- matrix(c(x1,x2,1,1),nrow=2,ncol=2)
  b <- matrix(c(y1,y2),nrow=2,ncol=1)
  result<-solve(a,b)
  (y-result[2])/result[1]
}
LDDecay <- data.frame(label=colnames(myScale)[2:dim(myScale)[2]],LDDecayValue=NA)

for (j in 2:(length(locs)+2)){
  
  if(min(which(myScale[,j]<=dividingLine[1])) ==Inf){
    stop(paste0("!!! The threshold is too low, please up to  ",min(myScale[,j]),"\n(dividingLine)\n"))
  }
  
  a <- max(which(myScale[,j]>dividingLine[1]),na.rm=T)
  b <- min(which(myScale[,j]<=dividingLine[1]),na.rm=T)
  ab <- data.frame(scale=myScale[c(a,b),1],aver=myScale[c(a,b),j])
  if (names(myScale[j])=="Mean"){
    d <- ab
    y3 <- myScale[a,j]
    y4 <- myScale[b,j]
  }
  x1 <- myScale$Scale[a]
  x2 <- myScale$Scale[b]
  y1 <- myScale[a,j]
  y2 <- myScale[b,j]
  
  LDPoint <- mine(x1,x2,y1,y2,dividingLine[1])
  LDDecay[j-1,2] <- LDPoint
}
if (unit=="Kb"){
  LDDecay$LDDecayValue <- round(LDDecay$LDDecayValue,2)
}else{
  LDDecay$LDDecayValue <- round(LDDecay$LDDecayValue)
}

LDPoint <- LDDecay$LDDecayValue[dim(LDDecay)[1]]
midax <- myScale$Scale[dim(myScale)[1]]/2-strwidth("r2=0.1, xxxx Kb")/2
text(midax,0.5, expression(r^2),adj = c(0,0))
srt<-strwidth("r2")
text(midax+srt,0.5,paste("= ",dividingLine[1],", ",LDPoint," ",unit,sep=""),adj=c(0,0))
cat(paste("The estimate value of LD decay is",round(LDDecay$LDDecayValue[dim(LDDecay)[1]],2),"\n"))

dev.new()
plot(myScale$Scale,myScale$Mean,xaxt="n",yaxt="n",type=dot,lwd=lineWidth, ylim=c(0,1),col=lineColor[1],main=paste(plotTitle,sep=' '), xlab=expression("Physical Distance (Kb)"), ylab=expression("r"^" 2"))
if(unit=="Kb"){
  axis(side=1,at=acc/1000)
  axis(side=1,at=ad/1000,labels=FALSE)
}else{
  axis(side=1,at=acc)
  axis(side=1,at=ad/1000,labels=FALSE)
}
axis(side=2,at=seq(0,1,0.1))
if (dottedLine == TRUE){
  abline(h=dividingLine,lty=2,col="red")
}
for (i in 2:(ncol(myScale)-1)){
  lines(myScale[,1],myScale[,i],type=dot,col=lineColor[i])
}
if(plotChrR2==TRUE){
  legend("topright", c(paste(LDDecay$label,", r^2 = ",dividingLine,",",LDDecay$LDDecayValue,unit)), lty = 1, merge=TRUE,col = c(lineColor[2:11],lineColor[1]))
}else{
  legend("topright", c(paste(LDDecay$label)), lty = 1, merge=TRUE,col = c(lineColor[2:length(chrs)],lineColor[1]))
}

write1 <- myScale
names(write1)[1] <- paste("Scale Distance (",unit,")",sep='')
if (unit=="Kb"){
  write.csv(round(write1,2),paste(pop.name,plotTitle,".csv",sep=' '),row.names=FALSE)
}else{
  write.csv(round(write1,2),paste(pop.name,plotTitle,".csv",sep=' '),row.names=FALSE)
}
write2 <- LDDecay
names(write2) <- c("Chromosome",paste("LD Decay Distance (",unit,")",sep=''))
write.csv(write2,paste(pop.name,"LD Decay Chr.csv"),row.names=FALSE)

cat("Mission completed.\t任务完成！\n")
})
