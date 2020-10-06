library(LinkageMapView)
if(!"readr" %in% installed.packages()) {
  install.packages("LinkageMapView")
}

cat("Please choose a import file.\n")
f.info <- choose.files()
setwd(dirname(f.info))

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
f <- readFiles(header=TRUE,fname = basename(f.info)) # reading the file
fNameG <- paste0(basename(f.info),".map.pdf")
fNameD <- paste0(basename(f.info),".density.pdf")
lmv.linkage.plot(f,outfile = fNameG,autoconnadj = FALSE)
lmv.linkage.plot(f,outfile = fNameD,denmap=TRUE)
cat("The output directory is:",getwd(),"\n")
