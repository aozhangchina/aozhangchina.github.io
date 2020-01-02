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
#  $$$$$$$$/______    ______  $$ |              Tool name: GOMultiRowintoOneRow.r
#     $$ | /      \  /      \ $$ |
#     $$ |/$$$$$$  |/$$$$$$  |$$ |              Updated: 01/02/2020
#     $$ |$$ |  $$ |$$ |  $$ |$$ |
#     $$ |$$ \__$$ |$$ \__$$ |$$ |              Website: datahold.cn
#     $$ |$$    $$/ $$    $$/ $$ |
#     $$/  $$$$$$/   $$$$$$/  $$/ 
#                                 


rm(list=ls())
fn <- file.choose()
bn <- basename(fn)
nameOutput <- paste0(bn,"_transform.csv")
# set workspace
setwd(dirname(fn))

# read a CSV file to be converted
# The format of CSV file:
# RGZM2G010871	GO:0003700
#           		GO:0003701
#               GO:0003702
#               GO:0003703
#               GO:0003704
#               GO:0003705
#               GO:0003706
#RGZM2G010872	  GO:0003707
#               GO:0003708
#               GO:0003709
#RGZM2G010873	  GO:0003710
#               GO:0003711
#               GO:0003712
#               GO:0003713
test <- read.csv(bn,header=F)

# Complement the fist column
for (i in 1:nrow(test)){
  if(test$V1[i]==""){
    test$V1[i] <- test$V1[i-1]
  }
}

# Initialize the data frame of result
dataF <- data.frame(names=as.character(unique(test$V1)),value=NA)
dataF$names <- as.character(dataF$names)

for (j in 1:nrow(dataF)){
  
  dataF$value[j] <- paste(test$V2[which(test$V1==dataF$names[j])],collapse = ",")
}

cat("  ______                   ________  __                                     __        ")
cat(" /      \                 /        |/  |                                   /  |       ")
cat("/$$$$$$  |  ______        $$$$$$$$/ $$ |____    ______   _______    ______ $$/_______ ")
cat("$$ |__$$ | /      \           /$$/  $$      \  /      \ /       \  /      \ //       |")
cat("$$    $$ |/$$$$$$  |         /$$/   $$$$$$$  | $$$$$$  |$$$$$$$  |/$$$$$$  |/$$$$$$$/ ")
cat("$$$$$$$$ |$$ |  $$ |        /$$/    $$ |  $$ | /    $$ |$$ |  $$ |$$ |  $$ |$$      \ ")
cat("$$ |  $$ |$$ \ _$$ |       /$$/____ $$ |  $$ |/$$$$$$$ |$$ |  $$ |$$ \ _$$ | $$$$$$  |")
cat("$$ |  $$ |$$    $$/       /$$      |$$ |  $$ |$$    $$ |$$ |  $$ |$$    $$ |/     $$/ ")
cat("$$/   $$/  $$$$$$/        $$$$$$$$/ $$/   $$/  $$$$$$$/ $$/   $$/  $$$$$$$ |$$$$$$$/  ")
cat("                                                                  /  \ _$$ |          ")
cat("                                                                  $$    $$/           ")
cat("                                                                   $$$$$$/            ")
cat(" ________                   __ ")
cat("/        |                 /  |                                                       ")
cat("$$$$$$$$/______    ______  $$ |              Tool name: GOMultiRowintoOneRow.r")
cat("   $$ | /      \  /      \ $$ |")
cat("   $$ |/$$$$$$  |/$$$$$$  |$$ |              Updated: 01/02/2020")
cat("   $$ |$$ |  $$ |$$ |  $$ |$$ |")
cat("   $$ |$$ \ _$$ |$$ \ _$$ |$$ |              Website: datahold.cn")
cat("   $$ |$$    $$/ $$    $$/ $$ |")
cat("   $$/  $$$$$$/   $$$$$$/  $$/ ")
  

# Output file
write.table(dataF,nameOutput,row.names = F,col.names = F,sep = ",")
cat(paste("File directory:\n"),getwd(),"\n","File name:\n",nameOutput,"\n")