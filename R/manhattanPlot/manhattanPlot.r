#    ______                   ________  __                                     __        
#   /      \                 /        |/  |                                   /  |       
#  /$$$$$$  |  ______        $$$$$$$$/ $$ |____    ______   _______    ______ $$/_______ 
#  $$ |__$$ | /      \           /$$/  $$      \  /      \ /       \  /      \$//       |
#  $$    $$ |/$$$$$$  |         /$$/   $$$$$$$  | $$$$$$  |$$$$$$$  |/$$$$$$  |/$$$$$$$/ 
#  $$$$$$$$ |$$ |  $$ |        /$$/    $$ |  $$ | /    $$ |$$ |  $$ |$$ |  $$ |$$      \ 
#  $$ |  $$ |$$ \__$$ |       /$$/____ $$ |  $$ |/$$$$$$$ |$$ |  $$ |$$ \__$$ | $$$$$$  |
#  $$ |  $$ |$$    $$/       /$$      |$$ |  $$ |$$    $$ |$$ |  $$ |$$    $$ |/     $$/ 
#  $$/   $$/  $$$$$$/        $$$$$$$$/ $$/   $$/  $$$$$$$/ $$/   $$/  $$$$$$$ |$$$$$$$/  
#                                                                    /  \__$$ |          
#                                                                    $$    $$/           
#                                                                     $$$$$$/            
#   ________                   __ 
#  /        |                 /  |
#  $$$$$$$$/______    ______  $$ |              Tool name: manhattanPlot.r
#     $$ | /      \  /      \ $$ |
#     $$ |/$$$$$$  |/$$$$$$  |$$ |              Updated: 05/09/2023
#     $$ |$$ |  $$ |$$ |  $$ |$$ |
#     $$ |$$ \__$$ |$$ \__$$ |$$ |              Website: datahold.cn
#     $$ |$$    $$/ $$    $$/ $$ |
#     $$/  $$$$$$/   $$$$$$/  $$/ 
#    

if (!exists("threshold_line")){
  threshold_line <- NULL   # threshold = 6 LOD or NULL
}


if(!"CMplot" %in% installed.packages()) {
  install.packages("CMplot")
  library(CMplot)
}else{
  library(CMplot)
}

# The funtion of reading files
readFiles <- function(header=TRUE,choose=FALSE,fname){
  if(!"readr" %in% installed.packages()) {
    install.packages("readr")
    library(readr)
  }else{
    library(readr)
  }
  if (choose==TRUE){
    myFileName <- file.choose()
  }else{
    myFileName <- fname
  }
  if (grepl("\\.csv$",myFileName)){
    if (header==TRUE){
      myFile <- read_csv(myFileName, col_names = TRUE)
    }else{
      myFile <- read_csv(myFileName, col_names = FALSE)
    }
  }else{
    if (header==TRUE){
      myFile <- read_tsv(myFileName, col_names = TRUE)
    }else{
      myFile <- read_tsv(myFileName, col_names = FALSE)
    }
  }
  return(myFile)
}
# Reading File

cf <- choose.files()
fn <- basename(cf)
setwd(dirname(cf))
getwd()

a <- readFiles(header=TRUE,fname=fn) # reading the file

temp1 <- data.frame("SNP"=a$Marker,"Chromosome"=a$Chr,"Position"=a$Pos,"name"=a$Trait,"trait"=a$p)

temp2 <- temp1[order(temp1$name,temp1$Chromosome,temp1$Position),]

if (length(which(temp2$SNP=="None"))>0){
  temp2 <- temp2[-which(temp2$SNP=="None"),]   # delete the line of None
}

if (length(which(is.na(temp2$SNP)))>0){
  temp2 <- temp2[-which(is.na(temp2$SNP)),]
}

if (length(which(is.na(temp2$Chromosome)))>0){
  temp2 <- temp2[-which(is.na(temp2$Chromosome)),]
}

titles <- unique(temp2$name)

paintData <- temp2[,1:3]
paintData <- temp2[!duplicated(temp2[,1:3]),1:3]

cat("  ______                   ________  __                                     __        \n")
cat(" /      \\                 /        |/  |                                   /  |       \n")
cat("/$$$$$$  |  ______        $$$$$$$$/ $$ |____    ______   _______    ______ $$/_______ \n")
cat("$$ |__$$ | /      \\           /$$/  $$      \\  /      \\ /       \\  /      \\$//       |\n")
cat("$$    $$ |/$$$$$$  |         /$$/   $$$$$$$  | $$$$$$  |$$$$$$$  |/$$$$$$  |/$$$$$$$/ \n")
cat("$$$$$$$$ |$$ |  $$ |        /$$/    $$ |  $$ | /    $$ |$$ |  $$ |$$ |  $$ |$$      \\ \n")
cat("$$ |  $$ |$$ \\__$$ |       /$$/____ $$ |  $$ |/$$$$$$$ |$$ |  $$ |$$ \\__$$ | $$$$$$  |\n")
cat("$$ |  $$ |$$    $$/       /$$      |$$ |  $$ |$$    $$ |$$ |  $$ |$$    $$ |/     $$/ \n")
cat("$$/   $$/  $$$$$$/        $$$$$$$$/ $$/   $$/  $$$$$$$/ $$/   $$/  $$$$$$$ |$$$$$$$/  \n")
cat("                                                                  /  \\__$$ |          \n")
cat("                                                                  $$    $$/           \n")
cat("                                                                   $$$$$$/            \n")
cat(" ________                   __ \n")
cat("/        |                 /  |                                                       \n")
cat("$$$$$$$$/______    ______  $$ |              Tool name: manhattanPlot.r\n")
cat("   $$ | /      \\  /      \\ $$ |\n")
cat("   $$ |/$$$$$$  |/$$$$$$  |$$ |              Updated: 05/09/2023\n")
cat("   $$ |$$ |  $$ |$$ |  $$ |$$ |\n")
cat("   $$ |$$ \\__$$ |$$ \\__$$ |$$ |              Website: datahold.cn\n")
cat("   $$ |$$    $$/ $$    $$/ $$ |\n")
cat("   $$/  $$$$$$/   $$$$$$/  $$/ \n\n")

tp <- NULL

if (length(unique(temp2$name))>1){
  for (i in 1:length(unique(temp2$name))){
    cat(paste0("trait",i,"=",unique(temp2$name)[i],"\n"))
    tp <- c(tp,unique(temp2$name)[i])
    indexI <- which(temp2$name==unique(temp2$name)[i])
    txt <- paste0("paintData$trait",i,"<-temp2$trait[indexI]")
    eval(parse(text=txt))
    names(paintData)[3+i] <- titles[i]
  }
}else{
  txt <- paste0("paintData$trait","<-temp2$trait")
  eval(parse(text=txt))
  paintData$trait <- temp2$trait
  names(paintData)[4] <- titles[1]
  tp <- titles[1]
}

CMplot(paintData,plot.type="m",LOG10=TRUE,main=titles,threshold=threshold_line,chr.den.col=NULL,file="jpg",file.name=tp,file.output=TRUE,verbose=TRUE)

b <- dir(pattern="Rect_Manhtn")
for (i in 1:length(b)){
  file.rename(b[i],paste0(fn,"_",tp[i],".jpg"))
}

CMplot(paintData, plot.type='q',conf.int.col=NULL, box=TRUE, file="jpg", file.name=tp,dpi=300, file.output= TRUE, verbose= TRUE)

bq <- dir(pattern="QQplot")
for (i in 1:length(bq)){
  file.rename(bq[i],paste0(fn,"_",tp[i],"_QQplot.jpg"))
}
