rm(list=ls())
# Read files
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

cat("Please choose a HMP file. \nPlease make sure your file dont't have NA.\n")
fDir <- choose.files()
setwd(dirname(fDir))
cat(paste0("Workspace:",getwd(),"\n"))
HMP_file <- readFiles(header=TRUE,choose=FALSE,basename(fDir))
# View(HMP_file)
HMP_file_name <- as.character(HMP_file[,1])
HMP_file_content <- as.data.frame(HMP_file[,12:ncol(HMP_file)])

# check NA
if (any(is.na(HMP_file_content))){
  stop("Have NA.\n")
}


# the function of get the same line number
myFuncGetSameLineNo <- function(y){
  x <- tempi
  length(which(x==y))/nrow(x)
}

# initiallize matrix
genetic_similarity <- NULL
for (i in 1:ncol(HMP_file_content)){
  cat(paste0(i,"/",ncol(HMP_file_content),"\n"))
  tempi <- HMP_file_content[i]
    tempii <- apply(HMP_file_content,2,myFuncGetSameLineNo)
    genetic_similarity <- rbind(genetic_similarity,tempii)
}
row.names(genetic_similarity) <- names(HMP_file_content)

write.table(genetic_similarity,paste0("genetic_similarity_",basename(fDir),".txt"),quote=FALSE,row.names = T,sep="\t")