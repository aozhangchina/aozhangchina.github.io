

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

cat("Please choose a HMP file. \nPlease make sure your file dont't have NA.")
fDir <- choose.files()
setwd(dirname(fDir))
cat(paste0("Workspace:",getwd(),"\n"))

# initialize stop sign
stopCal <- FALSE

# raed a hmp file
HMP_file <- readFiles(header=TRUE,choose=FALSE,basename(fDir))

# transform to a dataframe
HMP_file <- as.data.frame(HMP_file)

# get merket numbers
nMarkers <- nrow(HMP_file)

# transform to a character
HMP_file[,1] <- as.character(HMP_file[,1])

# get the name of markers
HMP_file_name <- c(HMP_file[,1])

# get Obtaining the genotype to be compared
x <- as.data.frame(HMP_file[,12])
standard_line <- x
  names(standard_line) <- names(HMP_file)[12]

# get data part of the HMP file
HMP_file_content <- HMP_file[,13:ncol(HMP_file)]

# the function of get the same line number
myFuncGetSameLineNo <- function(y){
  which(x==y)
}

# get the Same Line Number
SLM <- apply(HMP_file_content,2,myFuncGetSameLineNo)

# get the length of each SLM
SLML <- lapply(SLM,length)

# wheter 1 parent can satisfy
if (max(unlist(SLML))==nMarkers){
  cat("Just 1 parent is enough!\n")
  
  parent1 <- names(SLML[which(SLML==max(unlist(SLML)))])
  
  # creat a dataframe
  resultDataframe <- data.frame(parent1)
  
  # write file 
  write.table(resultDataframe,paste0("Possible_parents_",basename(fDir),".txt"),quote=FALSE,row.names = F,sep="\t")

  # get parents' genotype
  ParentsGeno <- cbind(standard_line,HMP_file_content[parent1])
  
  # write file
  write.table(ParentsGeno,paste0("Parents_Geno_",basename(fDir),".txt"),quote=FALSE,row.names = F,sep="\t")
  
  stopCal = TRUE   # No more counting.
}

# sort the file
orderHMP_file_content <- HMP_file_content[,order(unlist(SLML),decreasing = T)]
  
if (method == "fast" & stopCal != TRUE){
  
    # get the materials with the most line numbers as a new x
    x <- SLM[names(orderHMP_file_content[1])]

    # SLM
    y <- SLM

    # the function of union
    myFuncUnion <- function(y)
    {
      union(unlist(x),unlist(y))
    }

    # get line numbers's union
    SLMUnion <- lapply(y,myFuncUnion)

    # get totoal number of SLMUnion
    SLMUnionL <- unlist(lapply(SLMUnion,length))

    # parent 1 and parent 2
    parent1 <- names(orderHMP_file_content[1])
    parent2 <- names(SLMUnionL)[which(SLMUnionL==max(SLMUnionL))]
    percentage <- max(SLMUnionL)/nMarkers

    # creat a dataframe
    resultDataframe <- data.frame(parent1,parent2,percentage)
           
    # write file 
    write.table(resultDataframe,paste0("Possible_parents_",basename(fDir),".txt"),quote=FALSE,row.names = F,sep="\t")

    # get parents' genotype
    parentsGeno <- cbind(standard_line,HMP_file_content[c(parent1,parent2)])

    # write file
    write.table(parentsGeno,paste0("Parents_Geno_",basename(fDir),".txt"),quote=FALSE,row.names = F,sep="\t")
    
}

if (method == "fast" & parents >= 3 & stopCal != TRUE){
  
  formax <- parents-2
  parentsNames <- c(parent1,parent2[1])
  addParent <- NULL
  for (ii in 1:formax){
  
    parentsNames <- c(parentsNames,addParent)
    
    # combine union
    y <- as.list(parentsGeno[parentsNames])
  
    # get line numbers's union
    SLMUnion3 <- Reduce(union,unlist(y))
  
    # delete NA
    SLMUnion3 <- SLMUnion3[!is.na(SLMUnion3)]
  
    percentage <- length(SLMUnion3)/nMarkers
  
    if(percentage == 1){
      stopCal = TRUE   # No more counting.
      stop("Parents are enough!\n")
    }else{
      y <- SLM
      x <- SLMUnion3
      
      # multiple parents same marker numbers 
      SLM3Union <- lapply(y,myFuncUnion)
      
      # get totoal number of SLMUnion
      SLM3UnionL <- unlist(lapply(SLM3Union,length))
      
      percentage <- max(SLM3UnionL)/nMarkers
      
      parents3 <- names(which(SLM3UnionL==max(SLM3UnionL)))
      parents3n <- which(names(orderHMP_file_content) %in% parents3)
      parents3n <- names(orderHMP_file_content)[parents3n[1]]
      
      # adding parents names
      addParent <- parents3n
      
      # get parents' genotype
      parentsGeno <- cbind(standard_line,HMP_file_content[c(parentsNames,addParent)])
      
      # creat a dataframe
      resultDataframe <- c(parentsNames,addParent,percentage)
      
      #names of resultDataframe
      names(resultDataframe) <- rep("parent",parents)
      names(resultDataframe)[length(resultDataframe)] <- "percentage"
    }
  }
  
  # write file 
  write.table(t(resultDataframe),paste0("Possible_parents_",basename(fDir),".txt"),quote=FALSE,row.names = F,sep="\t")
  write.table(parentsGeno,paste0("Parents_Geno_",basename(fDir),".txt"),quote=FALSE,row.names = F,sep="\t")
}

if (method == "full" & stopCal != TRUE){
  # All combinations
  n <- ncol(orderHMP_file_content)
  k <- choose(n, parents)
  cat(paste0("In totol should calculation ",k," times.\n"))
  allTimes <- combn(n,parents)
  # initialize percentage
  percentage <- 0
  for (l in 1:k){
    temp <- SLM[allTimes[,l]]
    tempUnion <- Reduce(union,unlist(temp))
    percentage1 <- length(tempUnion)/nMarkers
    
    if (percentage != 1 & percentage1>percentage){
      percentage <- percentage1
      templ <- temp
      cat(paste0(percentage,"\n"))
    }else if(percentage==1){
      parentsNames <- names(templ)
      orderHMP_file_content[parentsNames]
      
      resultDataframe <- c(parentsNames,percentage)
      #names of resultDataframe
      names(resultDataframe) <- rep("parent",parents)
      names(resultDataframe)[length(resultDataframe)] <- "percentage"
      parentsGeno <- cbind(standard_line,HMP_file_content[c(parentsNames)])
      
      # write file 
      write.table(t(resultDataframe),paste0("Possible_parents_",basename(fDir),".txt"),quote=FALSE,row.names = F,sep="\t")
      write.table(parentsGeno,paste0("Parents_Geno_",basename(fDir),".txt"),quote=FALSE,row.names = F,sep="\t")
      
      stop("percentage is 1.\n")
    }
    
    parentsNames <- names(templ)
    orderHMP_file_content[parentsNames]
    
    resultDataframe <- c(parentsNames,percentage)
    #names of resultDataframe
    names(resultDataframe) <- rep("parent",parents)
    names(resultDataframe)[length(resultDataframe)] <- "percentage"
    parentsGeno <- cbind(standard_line,HMP_file_content[c(parentsNames)])
  }
  # write file 
  write.table(t(resultDataframe),paste0("Possible_parents_",basename(fDir),".txt"),quote=FALSE,row.names = F,sep="\t")
  write.table(parentsGeno,paste0("Parents_Geno_",basename(fDir),".txt"),quote=FALSE,row.names = F,sep="\t")
}


