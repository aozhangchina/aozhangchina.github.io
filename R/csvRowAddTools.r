combineFile <- function(){
tempfile <- choose.files()
for (i in 1:length(tempfile)){
  ttt <- paste0("file",i,"<- read.csv(tempfile[i],header=T,stringsAsFactors = F)")
  eval(parse(text=ttt))
}

combineFile <- NULL
for (i in 1:length(tempfile)){
  ttt <- paste0("combineFile <- rbind(combineFile,file",i,")")
  eval(parse(text=ttt))
}
return(combineFile)
}