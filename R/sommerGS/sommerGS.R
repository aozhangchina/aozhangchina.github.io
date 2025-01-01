if(!exists("traitName")){traitName <- "BLUE_Fv_"}
if(!exists("taxa_name_in_Pheno")){taxa_name_in_Pheno <- "Genotype"}
if(!exists("Env")){Env <- NULL}   #"Environment"
if(!exists("GBLUP")){GBLUP <- TRUE}
if(!exists("rrBLUP")){rrBLUP <- FALSE}
if(!exists("rrBLUPr")){rrBLUPr <- TRUE}
if(!exists("GbyEMain")){GbyEMain <- FALSE}
if(!exists("GbyEDG")){GbyEDG <- FALSE}
if(!exists("cycles")){cycles <- 100}
if(!exists("get.train_test")){get.train_test <- FALSE}
if(!exists("A_pred_B")){A_pred_B <- FALSE}   # TRUE means a population predict b population
if(!exists("group")){group <- NULL}   # c("SYNDH","DTMA")
if(!exists("train_pop")){train_pop <- "SYNDH"}
if(!exists("test_pop")){test_pop <- "SYNDH"}
if(!exists("K_fold")){K_fold <- 5}
if(!exists("optimization")){optimization <- NULL}   # NULL,"CDmean"
if(!exists("NTrn.Optmize")){NTrn.Optmize <- 20}   # opt trn size
if(!exists("nIter")){nIter <- 10000}   # iterative number for CDmean
if(!exists("WD")){WD <- NULL}
source("https://dataholdcn.cn/R/GBIT/GBIT.R")

if (is.null(WD)){
  cat("Choose a work direction.\n")
  GBIT.setwd(choose.dir())
}else{
  GBIT.setwd(WD)
}

if("Pheno.data" %in% dir()&"markers.imputed.data" %in% dir()){
  cat("Discover handled files, skip several steps.\n")
  load(file="Pheno.data")
  load("markers.imputed.data")
}else{
  if (is.null(sommerGenof)){
    cat("Choose a genotype file.\n")
    sommerGenof <- choose.files()
  }
  
  if(is.null(sommerPhenof)){
    cat("Choose a phenotypic file.\n")
    sommerPhenof <- choose.files()
  }
  
  
  sommerGeno <- GBIT.readFile(choose = F,header = T,fname = sommerGenof)
  
  sommerPheno <- GBIT.readFile(choose = F,header = T,fname = sommerPhenof)
  
  # geno part
  if(is.character(sommerGeno[2,12])&grep(".hmp|.HMP",basename(sommerGenof))==1){
    cat("This is a HMP file.\nWe are starting a numerical conversion process.\n")
    mySommerGeno <- GBIT.initiate.geno(sommerGeno)
    genoG <- mySommerGeno$genoG
    row.names(genoG) <- genoG[,1]
    markers <- apply(genoG[,-1], 1, GBIT.geno.numerical)
    row.names(markers) <- names(genoG[,-1])
  }else{
    markers <- sommerGeno
    rownames(markers) <- markers[,1]
    if (!is.numeric(markers[1,1])){
      markers <- markers[,-1]
    }
    markers <- as.matrix(markers)
  }
  # Remove big file
  if(exists("mySommerGeno")){
    rm(mySommerGeno)
  }
  
  
  if(any(is.na(markers))){
    cat("Missing values exist.\nWe are starting an imputing process.\n")
    markers_imputed <- apply(markers,2,GBIT.impute.mean)
    row.names(markers_imputed) <- names(genoG[,-1])
  }else{
    markers_imputed <- markers
  }
  # Remove big file
  if(exists("markers")){
    rm(markers)
  }
  
  # pheno part
  Name_trait <- grep(traitName,names(sommerPheno))
  if (length(Name_trait)==0){
    stop("Please check the traitName variable.\n")
  }
  if (any(is.na(sommerPheno[,1]))){
    sommerPheno <- sommerPheno[-which(is.na(sommerPheno[,1])),]
  }
  
  # merge Pheno and Geno
  if(!is.null(group)&length(group)>0){
    sommerPheno <- sommerPheno[sommerPheno$group %in% group,]
  }
  cat("Pheno:", length(unique(sommerPheno[, taxa_name_in_Pheno])), "\nGeno:", length(row.names(markers_imputed)),"\n")
  index <- intersect(unique(sommerPheno[, taxa_name_in_Pheno]), row.names(markers_imputed))
  cat("Intersection:", length(index),"\n")
  sommerPheno <-  sommerPheno[which(sommerPheno[, taxa_name_in_Pheno] %in% index), ]
  markers_imputed <-  markers_imputed[which(row.names(markers_imputed) %in% index), ]
  cat("Pheno:", length(unique(sommerPheno[, taxa_name_in_Pheno])), "\nGeno:", length(row.names(markers_imputed)),"\n")
  if (!identical(unique(sommerPheno[, taxa_name_in_Pheno]), row.names(markers_imputed)) &
      is.null(Env)) {
    sommerPheno <-
      sommerPheno[order(unique(sommerPheno[, taxa_name_in_Pheno])), ]
    markers_imputed <-
      markers_imputed[order(row.names(markers_imputed)), ]
  }
  if (identical(unique(sommerPheno[, taxa_name_in_Pheno]), row.names(markers_imputed))) {
    cat("The Geno and Pheno are identical. Pass!\n")
  } else{
    cat("Please check the Geno and Pheno files.\n")
  }
  save(markers_imputed,file = "markers.imputed.data")
  save(sommerPheno,file = "Pheno.data")
}


## CDmean
if (optimization=="CDmean"){
  cat("CDmean process.\n")
  Reference <- sommerPheno[,taxa_name_in_Pheno]
  CDmean_tst <- which(sommerPheno$group==test_pop)
  tst <- sommerPheno[CDmean_tst,taxa_name_in_Pheno]
  CDmean_trn <- which(!sommerPheno$group==test_pop)
  ##! CDMean
  
  trn_tst_Raw <- sommerPheno[,taxa_name_in_Pheno]   # names
  
  if(is.unsorted(trn_tst_Raw)){
    cat("The genotypic names are not sored.\n")
    trn_tst_Raw <- trn_tst_Raw[order(trn_tst_Raw)]  # sort
  }
  
  Gmatrix <- GBIT.get.kinship.VanRaden(markers_imputed)
  
  tst_size <- length(CDmean_tst)
  
  if (!is.null(NTrn.Optmize)) {
    trn_size = NTrn.Optmize
  }
  
  
  trn_new <- sample(CDmean_trn, NTrn.Optmize, replace = F)   # to sample a first trn
  trn_tst <- c(CDmean_tst, trn_new)   # 濞戞挸瀛╃粊瀵告嫚閺囥垺鑲犻柛姘墕閹酣鐛�?
  
  trn_tst <- Reference[trn_tst]   # 鐎电増顨呴崺宀勫春閸濆嫭绀堥柛銊ヮ儏閹洜绮�?
  
  Gmatrix_new <- Gmatrix[trn_tst, trn_tst]   # 闁圭粯鍔曡ぐ鍥矗椤忓懏绠掓繛鏉戭儓閻︻垶姊块崱娆愮暠kinship
  
  invAt_new <- solve(Gmatrix_new)   # 婵炴潙顑堥惁顖炴⒖閸℃鍊ゅù婊嗗皺缁鳖參宕楀畷鍥厙闁活厸鏅犲Ο鈧柛娆愮墵閳ь剙妫涢悡鈺呮⒓?
  
  # TT闁哄嫷鍨伴顔夹掗弮鍌滃彁闂傚啰顣槐婵喰掗弮鍥╃獩
  TT <- matrix(0,length(rownames(Gmatrix_new)),tst_size)  # 鐎点倛娅ｉ悵娑欑▔閳ь剚绋�?0闁活厸鏅犲Ο鈧柕鍡楀€块弳閬嶅及椤栨稓銈撮悹鍥ㄦ礋濞夘洃inship闁汇劌瀚、鎴︽晬鐏炶棄鐏欓柡鍕靛灠鐢偄霉鐎ｎ厾妲搁梻鍡楁閵囧洨浜�?
  TT[1:tst_size,] <- -1/tst_size   # 婵絽绻楅、?= -1/闁告ḿ鍠愮粊瀵告嫚閺囥垺鑲犲鍫嗗啰姣�
  # 閻庣數顢婇～妤冪棯鐠恒劍鏆�1-闁活厸鏅犲Ο鈧�
  for (i in 1:tst_size) { 
    TT[i,i]=1-1/tst_size 
  }     
  dim(TT) 
  
  X_trn <- rep(1,trn_size)  # 闂佹彃绉撮ˇ鍙夊濡搫顕х€点倗鍎よ啯闂傚棗妫楅幃搴♀枎?1
  Ident_trn <- diag(trn_size)   # 鐎点倛娅ｉ悵娑氣偓鐢殿攰椤鐥幐搴⑿�1闁汇劌瀚悡鈺呮⒓?
  M_trn <- Ident_trn - (X_trn %*% solve(t(X_trn) %*% X_trn) %*% t(X_trn) )  # 闁兼儳鍢茬欢杈ㄥ濡搫顕х€点倗鍎よ啯闂傚棗妫楅幃搴ㄦ儍閸曨噮鍟庨悹渚婄磿閻撯晠姊�?
  dim(M_trn)
  
  # 閻炴凹鍋勭紓鎾澄熼敍鍕暠閻犱焦宕橀鎼佹儗閳哄懏鈻�
  Z_trn <- matrix(0,trn_size,trn_size+tst_size) 
  for (i in 1:trn_size){
    Z_trn[i,tst_size+i]=1 
  }
  dim(Z_trn)
  
  VarE <- 0.5
  VarG <- 0.5
  
  lambda= (VarE/VarG)
  
  RawCD <- (t(TT)%*%(Gmatrix_new-lambda*solve(t(Z_trn)%*%M_trn%*%Z_trn + lambda*invAt_new))%*%TT)/(t(TT)%*%Gmatrix_new%*%TT) 
  
  CD <- diag(RawCD)
  # 闁告帗绻傞～鎬孌闁秆冩搐閳�??
  CDmeanSave <- mean(CD)
  CDmeanSaveCopy <- mean(CD)
  
  CDmeanExchange <- rep(NA,nIter)   # CDMean閺夆晩鍘洪崬顒勬儗閳哄懏鈻�
  CDmeanExchange[1] <- CDmeanSave   # 1闁告瑨娓圭紞鍛磾椤旂厧鐏ュ┑顔碱儏閳ь剟妫跨拹鐑燚Mean闁秆冩搐閳�??
  
  cpt2 <- 1
  cpt <- 0
  
  while (cpt2<=nIter) {
    cpt2 <- cpt2 + 1
    # 闁兼儳鍢茬欢杈ㄧ▔瀹ュ懏韬梻鍡楁閹孩绋夐鐘崇暠闁哄鍔栭弸?
    PartOfTrainingSet <- CDmean_trn[-match(trn_new,CDmean_trn)] 
    
    # 闁革负鍔岀紓鎾澄熼敓鐘宠偁闁告艾鐗呴懙鎴︽⒕韫囨梹绨氶梺顐㈩槹鐎氥劍绋夐埀顒佺▔椤忓懏缍忛柡??
    Sample2 <- sample(trn_new,1)
    
    # 闁革负鍔戝ḿ顏勵嚈閻戞◥渚€姊块崱妤佸€ら梻鍛箲濠р偓闂侇偄顦扮€氥劍绋夐埀顒佺▔椤忓懏缍忛柡??
    Sample3 <- sample(PartOfTrainingSet,1)
    
    # 闁告ḿ绮敮鈧梻鍛箲濠р偓闂侇偄顦扮€氥劑鎯冮崟顐ょ处婵☆垽绻濆▔锕傚触閸垺鐣遍柡浣哄瀹撲線鏁嶇仦钘夋櫃闁告梻濮撮崣鍡樼▔閳ь剚绋夐鍛厐闁汇劌瀚弳鐔煎箲椤曞棛绀夌紓浣稿閸ㄦ岸寮幍顔界暠鐎点倗鍎よ啯闂傚棗妫楅幃?
    Sample4 <- c(Sample3,trn_new[trn_new!=Sample2])
    
    # 缂備礁瀚紓鎾诲棘閹殿喗鐣辨繛鏉戭儓閻︻垶姊块崱顓犵闁告ḿ鍠栧▔锕傚触?+闁哄倹澹嗗▓鎴澝圭€ｎ厾妲搁梻鍡楁閹酣鏁�?
    Sample5 <- c(CDmean_tst, Sample4)
    
    Sample5 <- Reference[Sample5]   # 鐎电増顨呴崺灞矫圭€ｎ厾妲搁梻鍡楁閹酣鎯冮崙顤廌闁�??
    
    ## 鐎电増顨呴崺瀛弇atrix, 闁告鐗癷nship闁活厸鏅犲Ο鈧�
    GmAt <- Gmatrix[Sample5, Sample5]
    
    ## 闁活厸鏅犲Ο鈧柛娆愮墵閳�??
    invAt1 <- solve(GmAt)
    
    ## 閻犱緤绱曢悾璇残掕箛搴ｎ伇閺夌儐鍠氬▓鎱岲闁秆冩搐閳�??
    RawCD2 <- (t(TT)%*%(GmAt-lambda*solve(t(Z_trn)%*%M_trn%*%Z_trn + lambda*invAt1))%*%TT)/(t(TT)%*%GmAt%*%TT)
    
    CD2 <- diag(RawCD2)
    CDmeanSave2 <- mean(CD2)	
    
    if (mean(CD2) > CDmeanSave ) { 
      trn_new <- Sample4 # Accept the new Calibration set if CDmean is increased, reject otherwise.
      CDmeanSave <- mean(CD2)  
    }
    CDmeanExchange[cpt2]=CDmeanSave
  } 
  
  jpeg("CDmeanExchange.jpg",width = 8, height = 5, units = "in", res = 300)
  
  
  plot(CDmeanExchange) # CDMean闁稿﹨鍋愰垾妯荤┍濠靛洦绠掗悺鎺戝暱椤у嫰鎯冮崟顔煎殹濞�??
  
  dev.off()
  
  SampleOptimiz <- trn_new # 濞村吋锚鐎垫煡姊块崱妤佸€�
  
  ### Subset customized training set from ref
  
  Result <- Reference[SampleOptimiz]
  
  Result <- data.frame(genotype= Result, stringsAsFactors = F)
  
  write.csv(Result, paste("GBIT\\OptimizedTRN", NTrn.Optmize,optimization, ".csv", sep = "_"), row.names = F)
  # write.csv(COR, paste("PredACC",NTrn.Optmize, optimization, ".csv", sep = "_"))
  
  ## End
}else if(optimization=="avgGRM"&!paste0("avgGRM_",test_pop,".csv") %in% dir("optimization")){
  cat("avgGRM process.\n")
  GBIT.library("reshape2")
  Reference <- sommerPheno[,taxa_name_in_Pheno]
  avgGRM_tst <- which(sommerPheno$group==test_pop)
  avgGRM_trn <- which(!sommerPheno$group==test_pop)
  ##! avgGRM
  
  TSTsample <- sommerPheno[avgGRM_tst,taxa_name_in_Pheno]
  TRNsample <- sommerPheno[avgGRM_trn,taxa_name_in_Pheno]
  
  if(class(Reference)!="data.frame") {
    Reference <- data.frame(GID=Reference, stringsAsFactors = F)
  }
  
  Result <- data.frame(GID=TRNsample, AVG=NA, stringsAsFactors = F)
  
  trn_tst_Raw <- as.character(sommerPheno[,taxa_name_in_Pheno])   # names 
  
  if(is.unsorted(trn_tst_Raw)){
    cat("The genotypic names are not sored.\n")
    trn_tst_Raw <- trn_tst_Raw[order(trn_tst_Raw)]  # sort
  }
  
  Gmatrix <- GBIT.get.kinship.VanRaden(markers_imputed)   # get kinship matrix using VanRaden
  if(is.unsorted(rownames(Gmatrix))){
    cat("The genotypes are not sored.\n")
    Gmatrix <- Gmatrix[trn_tst_Raw,trn_tst_Raw]  # sort
  }
  
  Result[,2] <- apply(Result[,1, drop=F], 1, function(x){
    
    trn_tst <- c(TSTsample, x)   # 瀵ょ儤膩闂嗗棗鎮庢晶鐐插娑撯偓娑擃亝娼楅弬?
    
    trn_tst <- trn_tst[order(trn_tst)]   # 閹烘帒绨�
    
    Reference <- subset(Reference, GID %in% trn_tst)   # 閹绘劕褰囬崣鍌濃偓鍐ц厬閻ㄥ嚕ID
    
    Reference <- Reference[order(Reference$GID), ,drop =F]   # 閹绘劕褰囬幋鎰殶閹诡喗顢�
    
    
    ### 閼惧嘲绶卞鐑樐侀梿鍡楁値閻ㄥ嫬鍙х化鑽ょ叐闂�?
    Gmatrix <- Gmatrix[trn_tst, trn_tst]
    
    ## 婢х偛濮濺ef閻ㄥ嫬鍨敍灞藉斧瀵ょ儤膩缂囥倓缍�1閿涘苯鍙炬禒鏍﹁礋2
    Reference$Value <- with(Reference, ifelse(GID %in%TSTsample , 1,2))   # 閺嶅洦鏁�1閿�?2
    refm <- data.frame(GID= Reference$GID, GRM=Reference$Value)
    refm$GID <- as.character(refm$GID)   # GID閸欐ɑ鍨氱€涙顑佹稉?
    ## 鐏忓棛鐓╅梼浣冩祮娑撹桨琚辨稉銈囩矋閸氬牏娈戣ぐ銏犵础 
    GmatrixL <- reshape2::melt(as.matrix(Gmatrix))   # reshape
    GmatrixL$Var1 <- as.character(GmatrixL$Var1)   # 鏉烆剚鍨氱€涙顑佹稉?
    GmatrixL$Var2 <- as.character(GmatrixL$Var2)   # 鏉烆剚鍨氱€涙顑佹稉?
    names(GmatrixL) <- c('rowL', 'colL', 'GRM')   # 婢х偛濮為崥宥囆�
    GmatrixLR <- merge(refm, GmatrixL, by.x = 'GID', by.y = 'rowL', all = TRUE)   # 閸栧綊鍘ょ粭顑跨閸掓绱濋惄绋跨秼娴滃孩绁寸拠鏇犲參娴�?
    GmatrixLR$GID <- as.character(GmatrixLR$GID)   # GID鏉烆剚鍨氱€涙顑侀崹?
    GmatrixLRC <- merge(refm, GmatrixLR, by.x = 'GID', by.y = 'colL', all = TRUE)   # 閸栧綊鍘ょ粭顑跨癌閸掓绱濋惄绋跨秼娴滃骸缂撳Ο锛勫參娴�?
    GmatrixLRC$GID <- as.character(GmatrixLRC$GID)   # GID鏉烆剚宕查幋鎰摟缁楋箑鐎�
    names(GmatrixLRC) <- c('TST', 'TSTpop', 'TRN', 'TRNpop', 'GRM')   # 缂佹瑦鐖ｆ０?
    ##Subset only the pairwise combination of the Testing set(TSTpop) and the Training set
    tmp <- subset(GmatrixLRC, TSTpop == 1 & TRN%in%x)   # 閹绘劕褰囧ù瀣槸缂囥倓缍嬫稉?1閿涘苯缂撳Ο锛勫參娴ｆ挷璐�2閻ㄥ嫨鈧�?
    
    AVG <- mean(tmp$GRM)   # GRM閻ㄥ嫬娼庨崐?
    
  })
  
  Result <- Result[order(-Result$AVG), , drop = F]
  GBIT.write.clipboard(Result)
  write.table(Result,paste0("optimization\\avgGRM_",test_pop,".csv"),sep=",",row.names = F,col.names = F)
  
  
}

if(exists("Result")&optimization=="CDmean"){
  # merge Pheno and Geno
  row.names(sommerPheno) <- sommerPheno[,taxa_name_in_Pheno]
  sommerPheno <- sommerPheno[c(tst,Result[,1]),]
  cat("Pheno:", length(unique(sommerPheno[, taxa_name_in_Pheno])), "\nGeno:", length(row.names(markers_imputed)),"\n")
  index <- intersect(unique(sommerPheno[, taxa_name_in_Pheno]), row.names(markers_imputed))
  cat("Intersection:", length(index),"\n")
  sommerPheno <-  sommerPheno[which(sommerPheno[, taxa_name_in_Pheno] %in% index), ]
  markers_imputed <-  markers_imputed[which(row.names(markers_imputed) %in% index), ]
  cat("Pheno:", length(unique(sommerPheno[, taxa_name_in_Pheno])), "\nGeno:", length(row.names(markers_imputed)),"\n")
  if (!identical(unique(sommerPheno[, taxa_name_in_Pheno]), row.names(markers_imputed)) & is.null(Env)) {
    sommerPheno <-
      sommerPheno[order(unique(sommerPheno[, taxa_name_in_Pheno])), ]
    markers_imputed <-
      markers_imputed[order(row.names(markers_imputed)), ]
  }
  if (identical(unique(sommerPheno[, taxa_name_in_Pheno]), row.names(markers_imputed))) {
    cat("The Geno and Pheno are identical. Pass!\n")
  } else{
    cat("Please check the Geno and Pheno files.\n")
  }
}else if (optimization=="CDmean"){
  stop("CDmean algorithm failed.\n")
}

if(optimization=="avgGRM"&paste0("avgGRM_",test_pop,".csv") %in% dir("optimization")){
  Result <- GBIT.readFile(header = F,choose = F,fname = paste0("optimization\\avgGRM_",test_pop,".csv"))
  names(Result) <- c("id","X1")
  Result <- Result[1:NTrn.Optmize,]
  tst <- sommerPheno[which(sommerPheno$group==test_pop),taxa_name_in_Pheno]
  sommerPheno2 <- sommerPheno
  # sommerPheno2 <- sommerPheno[,c(taxa_name_in_Pheno,traitName)]
  # sommerPheno2 <- sommerPheno2[which(tst %in% sommerPheno2[,taxa_name_in_Pheno]),]
  # 
  # sommerPheno3 <- sommerPheno2[,c(taxa_name_in_Pheno,traitName)]
  # row.names(sommerPheno3) <- sommerPheno2[,taxa_name_in_Pheno]
  cat("Pheno:", length(unique(sommerPheno2[, taxa_name_in_Pheno])), "\nGeno:", length(row.names(markers_imputed)),"\n")
  index <- intersect(unique(sommerPheno2[, taxa_name_in_Pheno]), row.names(markers_imputed))
  cat("Intersection:", length(index),"\n")
  sommerPheno <-  sommerPheno[which(sommerPheno[, taxa_name_in_Pheno] %in% index), ]
  markers_imputed <-  markers_imputed[which(row.names(markers_imputed) %in% index), ]
  cat("Pheno:", length(unique(sommerPheno[, taxa_name_in_Pheno])), "\nGeno:", length(row.names(markers_imputed)),"\n")
  if (!identical(unique(sommerPheno[, taxa_name_in_Pheno]), row.names(markers_imputed)) & is.null(Env)) {
    sommerPheno <-
      sommerPheno[order(unique(sommerPheno[, taxa_name_in_Pheno])), ]
    markers_imputed <-
      markers_imputed[order(row.names(markers_imputed)), ]
  }
  if (identical(unique(sommerPheno[, taxa_name_in_Pheno]), row.names(markers_imputed))) {
    cat("The Geno and Pheno are identical. Pass!\n")
  } else{
    cat("Please check the Geno and Pheno files.\n")
  }
}

# GS
GBIT.library("sommer")
K <- A.mat(markers_imputed)
Name_trait <- grep(traitName,names(sommerPheno))
sommerPheno$X1 <- sommerPheno[,Name_trait]

sommerPheno$id <- as.factor(sommerPheno[,taxa_name_in_Pheno])
if(is.null(Env)){
  row.names(sommerPheno) <- sommerPheno$id
}else{
  names(sommerPheno)[names(sommerPheno)==Env] <- "Env"
}

set.seed(12345)
GS_cor_GBLUP <- NULL
GS_cor_GBLUPr <- NULL
GS_cor_rrBLUP <- NULL
GS_cor_GbyEMean <- NULL
GS_cor_GbyEDG <- NULL
GS_cor_rrBLUPr <- NULL
for(i in 1:cycles){
  if (exists("Result")&!is.null(train_pop)&!is.null(test_pop)&(optimization=="CDmean"|optimization=="avgGRM")){
    # CDmean or avgGRM
    vv <- sommerPheno[sommerPheno$group == test_pop,taxa_name_in_Pheno]   # test_pop for CDmean
    y.trn <- sommerPheno[sommerPheno$id %in% c(vv,Result[,1]),]
    row.names(y.trn) <- y.trn$id
    y.trn$id <- as.factor(as.character(y.trn$id))
    y.trn[vv,"X1"] <- NA
    head(y.trn)
    # if (is.unsorted(K))
    #   K <- K[sort(c(vv,Result[,1])),sort(c(vv,Result[,1]))]
  }else if (!is.null(group)&!is.null(train_pop)&!is.null(test_pop)){
    # A pre B
    # 寤烘ā缇や綋鎶芥牱
    vv <- sommerPheno[sommerPheno$group==test_pop,taxa_name_in_Pheno]   # test_pop as a testing population
    # 棰勬祴缇や綋
    y.trn <- sommerPheno
    row.names(y.trn) <- y.trn$id
    y.trn$id <- as.factor(as.character(y.trn$id))
    y.trn[vv,"X1"] <- NA
    head(y.trn)
  }else if(is.null(Env)){
    # NoGbyE
    # 寤烘ā缇や綋鎶芥牱
    vv <- sample(rownames(sommerPheno),round(nrow(sommerPheno)/K_fold))   # 鎶藉彇20%浣滀负楠岃瘉缇や綋锛堣棰勬祴锛�
    # 棰勬祴缇や綋
    y.trn <- sommerPheno
    y.trn[vv,"X1"] <- NA
    head(y.trn)
  }else{
    # GbyE
    sommerPheno_mean <- data.frame('X1'=tapply(sommerPheno$X1,sommerPheno$id,mean))
    sommerPheno_mean$id <- row.names(sommerPheno_mean)
    y.trn <- sommerPheno_mean
    vv <- sample(rownames(sommerPheno_mean),round(nrow(sommerPheno_mean)/K_fold))
    y.trn[vv,"X1"] <- NA
    y.trn2 <- sommerPheno   # 澶氱幆澧冩暟鎹�
    y.trn2[which(sommerPheno$id %in% vv),"X1"] <- NA
    head(y.trn2)
  }
  
  
  # Write testing population
  if (get.train_test==TRUE){
    write.table(t(vv),"GBIT/GS_testPop.csv",append = T,row.names = F,col.names = F,sep = ",",quote=F)
  }
  ## GBLUP
  if(GBLUP&is.null(Env)){
    ans_GBLUP <- mmer(X1~1,
                      random=~vs(id,Gu=K),
                      rcov=~units,
                      data=y.trn, verbose = FALSE) # kinship based
    ans_GBLUP$U$`u:id`$X1 <- as.data.frame(ans_GBLUP$U$`u:id`$X1)
    GS_cor_GBLUP <- c(GS_cor_GBLUP,cor(ans_GBLUP$U$`u:id`$X1[vv,],sommerPheno[vv,"X1"], use="complete"))
    
    if (length(optimization)>0){
      write.table(GS_cor_GBLUP,paste0("GBIT/GBLUP_",optimization,"_",NTrn.Optmize,"_",test_pop,".cor.csv"),row.names = F,col.names = F,sep = ",",quote=F)
    }else{
      write.table(GS_cor_GBLUP,paste0("GBIT/GBLUP_",optimization,"_",NTrn.Optmize,".cor.csv"),row.names = F,col.names = F,sep = ",",quote=F)
    }
    
    if(i==cycles){
      cat("GBLUP:",GS_cor_GBLUP,"\n")
    }
  }
  
  
  
  ## rrBLUP
  if(rrBLUP&is.null(Env)){
    ans_rrBLUP <- mmer(X1~1,
                       random=~vs(list(markers_imputed)),
                       rcov=~units,
                       data=y.trn, verbose = FALSE) # kinship based
    u <- markers_imputed %*% as.matrix(ans_rrBLUP$U$`u:markers_imputed`$X1) # BLUPs for individuals
    rownames(u) <- rownames(markers_imputed)
    GS_cor_rrBLUP <- c(GS_cor_rrBLUP,cor(u[vv,],sommerPheno[vv,"X1"])) # same correlation
    if (length(optimization)>0){
      write.table(ans_rrBLUP,paste0("GBIT/rrBLUP_",optimization,"_",NTrn.Optmize,"_",test_pop,".cor.csv"),row.names = F,col.names = F,sep = ",",quote=F)
    }else{
      write.table(ans_rrBLUP,paste0("GBIT/rrBLUP_",optimization,"_",NTrn.Optmize,".cor.csv"),row.names = F,col.names = F,sep = ",",quote=F)
    }
    if(i==cycles){
      cat("rrBLUP:",GS_cor_rrBLUP,"\n")
    }
  }
  
  ## rrBLUP using rrBLUP package
  if(rrBLUPr&is.null(Env)){
    GBIT.library("rrBLUP")
    ww <- setdiff(y.trn$id,vv)
    y <- y.trn$X1
    names(y) <- y.trn$id
    Pheno_train <- as.matrix(y[ww])
    m_train <- markers_imputed[ww,]
    Pheno_valid <- as.matrix(sommerPheno[vv,"X1"])
    m_valid <- markers_imputed[vv,]
    ans_rrBLUPr <- mixed.solve(Pheno_train, Z=NULL, K=m_train, SE = FALSE, return.Hinv=FALSE)
    e <- as.matrix(ans_rrBLUPr$u)
    rrBLUPr_valid <- m_valid %*% e
    GS_cor_rrBLUPr <- c(GS_cor_rrBLUPr,cor(rrBLUPr_valid, Pheno_valid, use="complete"))
    if (length(optimization)>0){
      write.table(GS_cor_rrBLUPr,paste0("GBIT/rrBLUPr_",optimization,"_",NTrn.Optmize,"_",test_pop,".cor.csv"),row.names = F,col.names = F,sep = ",",quote=F)
    }else{
      write.table(GS_cor_rrBLUPr,paste0("GBIT/rrBLUPr_",optimization,"_",NTrn.Optmize,".cor.csv"),row.names = F,col.names = F,sep = ",",quote=F)
    }
    if(i==cycles){
      cat("rrBLUPr:",GS_cor_rrBLUPr,"\n")
    }
  }

  ## GBLUP using rrBLUP package
  if(GBLUPr&is.null(Env)){
    GBIT.library("rrBLUP")
    ww <- setdiff(y.trn$id,vv)
    y <- y.trn$X1
    names(y) <- y.trn$id
    Pheno_train <- as.matrix(y[ww])
    m_train <- markers_imputed[ww,]
    Pheno_valid <- as.matrix(sommerPheno[vv,"X1"])
    m_valid <- markers_imputed[vv,]
    ans_GBLUPr <- mixed.solve(Pheno_train, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
    e <- as.matrix(ans_GBLUPr$u)
    GBLUPr_valid <- m_valid %*% e
    GS_cor_GBLUPr <- c(GS_cor_GBLUPr,cor(GBLUPr_valid, Pheno_valid, use="complete"))
    if (length(optimization)>0){
      write.table(GS_cor_GBLUPr,paste0("GBIT/GBLUPr_",optimization,"_",NTrn.Optmize,"_",test_pop,".cor.csv"),row.names = F,col.names = F,sep = ",",quote=F)
    }else{
      write.table(GS_cor_GBLUPr,paste0("GBIT/GBLUPr_",optimization,"_",NTrn.Optmize,".cor.csv"),row.names = F,col.names = F,sep = ",",quote=F)
    }
    if(i==cycles){
      cat("GBLUPr:",GS_cor_GBLUPr,"\n")
    }
  }

  ## GbyEMain
  if(!is.null(Env)&GbyEMain){
    ansMain <- mmer(X1~Env,
                    random= ~ vs(id, Gu=K),
                    rcov= ~ units,
                    data=y.trn2, verbose = FALSE)
    ansMain$U$`u:id`$X1 <- as.data.frame(ansMain$U$`u:id`$X1)
    write.table(data.frame(ansMain$U$`u:id`$X1,i=i),"GBIT/GbyE_GEBV.csv",append=T,row.names = T,col.names = F,sep = ",",quote=F)
    GS_cor_GbyEMean <- c(GS_cor_GbyEMean,cor(ansMain$U$`u:id`$X1[vv,],sommerPheno_mean[vv,"X1"], use="complete"))
    write.table(GS_cor_GbyEMean,"GBIT/GbyEMean.cor.csv",row.names = F,col.names = F,sep = ",",quote=F)
    if (i ==cycles){
      cat("GbyEMean:",GS_cor_GbyEMean)
    }
  }
  
  # GbyEDG
  if(!is.null(Env)&GbyEDG){
    ansDG <- mmer(X1~Env,
                  random= ~ vs(ds(Env),id, Gu=K),
                  rcov= ~ units,
                  data=y.trn2, verbose = FALSE)
    summary(ansDG)
    as.data.frame(ansDG$U$`u:id`$X1)
  }
  
  if (A_pred_B==TRUE){
    break
  }
}
cat("Done!\n")

