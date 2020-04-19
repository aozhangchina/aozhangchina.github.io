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

compareTwoSets <- function(a,b){
  same_a_b <- length(a==b)
  percentage <- same_a_b/length(a)
  return(percentage)
}

cat("Please choose a HMP file. \nPlease make sure your file dont't have NA.")
fDir <- choose.files()
setwd(dirname(fDir))
cat(paste0("Workspace:",getwd(),"\n"))
HMP_file <- readFiles(header=TRUE,choose=FALSE,basename(fDir))
# View(HMP_file)
HMP_file_name <- as.character(HMP_file[,1])
HMP_file_content <- HMP_file[,12:ncol(HMP_file)]
# View(HMP_file_content)
# Initialize a matrix
genetic_similarity <- matrix(NA,ncol(HMP_file_content),ncol(HMP_file_content))
colnames(genetic_similarity) <- names(HMP_file_content)
rownames(genetic_similarity) <- names(HMP_file_content)
# View(genetic_similarity)
for (i in 1:ncol(genetic_similarity)){
    tempi <- HMP_file_content[i,]
    cat(paste0(i,"/",length(tempi),"\n"))
    apply(HMP_file_content,1,which(tempi==HMP_file_content))
    for (k in 1:ncol(genetic_similarity)){
        tempk <- HMP_file_content[k,]
        genetic_similarity_value <- length(which(tempi==tempk))/length(tempi)
        genetic_similarity[k,i] <- genetic_similarity_value
    }
}
write.table(genetic_similarity,paste0("genetic_similarity_",basename(fDir),".txt"),quote=FALSE,row.names = T,sep="\t")