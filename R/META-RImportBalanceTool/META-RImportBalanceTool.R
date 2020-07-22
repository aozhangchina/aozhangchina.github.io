cat("Please choose a CSV file.\n")
fn <- choose.files()
fnWithoutCSV <- sub(".csv","",basename(fn))

setwd(dirname(fn))
fileA <- read.csv(basename(fn))

txt <- paste0("varList <- unique((fileA$",var_env,"))")
eval(parse(text=txt))

finalMaterial <- NULL
for (i in 1:length(varList)){
  txt <- paste0("Env",varList[i]," <- unique(fileA$",var_material,"[which(fileA$",var_env,"=='",varList[i],"')])")
  eval(parse(text=txt))
  if (!i == 1){
    txt <- paste0("finalMaterial <- intersect(Env",varList[i],",Env",varList[(i-1)],")")
    eval(parse(text=txt))
  }
}

txt <- paste0("resultFile <- fileA[which(fileA$",var_material," %in% finalMaterial),]")
eval(parse(text=txt))
write.table(resultFile, paste0(dirname(fn),"/",fnWithoutCSV,"_balance.csv"), quote = F, row.names = F, sep = ",")
