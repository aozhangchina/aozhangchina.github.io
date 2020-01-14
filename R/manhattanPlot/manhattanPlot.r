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

cf <- choose.files()
fn <- basename(cf)
setwd(dirname(cf))
getwd()

a <- read.table(cf,header=T)

temp1 <- data.frame("SNP"=a$Marker,"Chromosome"=a$Chr,"Position"=a$Pos,"name"=a$Trait,"trait"=a$p)

temp2 <- temp1[order(temp1$Chromosome,temp1$Position),]

paintData <- temp2[,1:3]
paintData <- temp2[!duplicated(temp2[,1:3]),1:3]



if (length(levels(a$Trait))!=1){
for (i in 1:length(levels(a$Trait))){
  cat(paste0("trait",i,"=",levels(a$Trait)[i],"\n"))
  indexI <- which(temp2$name==levels(temp2$name)[i])
  txt <- paste0("paintData$trait",i,"<-temp2$trait[indexI]")
  eval(parse(text=txt))
}
}


CMplot(paintData,plot.type="m",LOG10=TRUE,threshold=NULL,chr.den.col=NULL,file="jpg",file.output=TRUE,verbose=TRUE)

