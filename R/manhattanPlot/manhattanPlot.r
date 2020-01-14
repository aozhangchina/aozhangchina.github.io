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
#  $$$$$$$$/______    ______  $$ |              Tool name: manhattanPlotforTassel.r
#     $$ | /      \  /      \ $$ |
#     $$ |/$$$$$$  |/$$$$$$  |$$ |              Updated: 01/13/2020
#     $$ |$$ |  $$ |$$ |  $$ |$$ |
#     $$ |$$ \__$$ |$$ \__$$ |$$ |              Website: datahold.cn
#     $$ |$$    $$/ $$    $$/ $$ |
#     $$/  $$$$$$/   $$$$$$/  $$/ 
#    


rm(list=ls())

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

temp2 <- temp1[order(temp1$Chromosome,temp1$Position),]
temp2 <- temp2[-which(temp2$SNP=="None"),]   # delete the line of None

paintData <- temp2[,1:3]
paintData <- temp2[!duplicated(temp2[,1:3]),1:3]

if (length(unique(a$Trait))>1){
for (i in 1:length(unique(a$Trait))){
  cat(paste0("trait",i,"=",unique(a$Trait)[i],"\n"))
  indexI <- which(temp2$name==unique(temp2$name)[i])
  txt <- paste0("paintData$trait",i,"<-temp2$trait[indexI]")
  eval(parse(text=txt))
}
}else{
  txt <- paste0("paintData$trait","<-temp2$trait")
  eval(parse(text=txt))
  paintData$trait <- temp2$trait
}


CMplot(paintData,plot.type="m",LOG10=TRUE,threshold=NULL,chr.den.col=NULL,file="jpg",memo=fn,file.output=TRUE,verbose=TRUE)

