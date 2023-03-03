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
  trn_tst <- c(CDmean_tst, trn_new)   # 娑撳孩绁寸拠鏇㈡肠閸氬牆鎮庨獮?
  
  trn_tst <- Reference[trn_tst]   # 瀵版鍩岄崺鍝勬礈閸ㄥ鎮曠粔?
  
  Gmatrix_new <- Gmatrix[trn_tst, trn_tst]   # 閹绘劕褰囬崣顏呮箒濞村鐦梿鍡欐畱kinship
  
  invAt_new <- solve(Gmatrix_new)   # 濞村鐦梿鍡楁値娴滆尙绱崗宕囬兇閻晠妯€閸欐牠鈧棛鐓╅梼?
  
  # TT閺勵垰顕В鏃傜叐闂冪绱濆В鏃囩窛
  TT <- matrix(0,length(rownames(Gmatrix_new)),tst_size)  # 瀵よ櫣鐝涙稉鈧稉?0閻晠妯€閵嗗倿鏆遍弰顖涚ゴ鐠囨洟娉inship閻ㄥ嫯顢戦敍灞藉灙閺勵垰甯ù瀣槸闂嗗棗銇囩亸?
  TT[1:tst_size,] <- -1/tst_size   # 濮ｅ繗顢?= -1/閸樼喐绁寸拠鏇㈡肠婢堆冪毈
  # 鐎电顫楃痪璺ㄦ暏1-閻晠妯€
  for (i in 1:tst_size) { 
    TT[i,i]=1-1/tst_size 
  }     
  dim(TT) 
  
  X_trn <- rep(1,trn_size)  # 闁插秴顦叉导妯哄瀵ょ儤膩闂嗗棗鎮庡▎?1
  Ident_trn <- diag(trn_size)   # 瀵よ櫣鐝涚€电顫楃痪鎸庢Ц1閻ㄥ嫮鐓╅梼?
  M_trn <- Ident_trn - (X_trn %*% solve(t(X_trn) %*% X_trn) %*% t(X_trn) )  # 閼惧嘲绶辨导妯哄瀵ょ儤膩闂嗗棗鎮庨惃鍕啎鐠侊紕鐓╅梼?
  dim(M_trn)
  
  # 鐞氼偄缂撳Ο锛勬畱鐠佹崘顓搁惌鈺呮█
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
  # 閸掓繂顫怌D閸у洤鈧??
  CDmeanSave <- mean(CD)
  CDmeanSaveCopy <- mean(CD)
  
  CDmeanExchange <- rep(NA,nIter)   # CDMean鏉╊厺鍞惌鈺呮█
  CDmeanExchange[1] <- CDmeanSave   # 1閸欒渹缍呯純顔煎灥婵鈧棿璐烠DMean閸у洤鈧??
  
  cpt2 <- 1
  cpt <- 0
  
  while (cpt2<=nIter) {
    cpt2 <- cpt2 + 1
    # 閼惧嘲绶辨稉宥呮躬闂嗗棗鎮庢稉顓犳畱閺夋劖鏋?
    PartOfTrainingSet <- CDmean_trn[-match(trn_new,CDmean_trn)] 
    
    # 閸︺劌缂撳Ο锟犳肠閸氬牅鑵戦梾蹇旀簚闁瀚ㄦ稉鈧稉顏呮綏閺??
    Sample2 <- sample(trn_new,1)
    
    # 閸︺劑娼鐑樐侀梿鍡楁値闂呭繑婧€闁瀚ㄦ稉鈧稉顏呮綏閺??
    Sample3 <- sample(PartOfTrainingSet,1)
    
    # 閸樼粯甯€闂呭繑婧€闁瀚ㄩ惃鍕紦濡繝娉﹂崥鍫㈡畱閺佺増宓侀敍灞藉晙閸旂姴鍙嗘稉鈧稉顏呮煀閻ㄥ嫭鏆熼幑顕嗙礉缂佸嫭鍨氶弬鎵畱瀵ょ儤膩闂嗗棗鎮?
    Sample4 <- c(Sample3,trn_new[trn_new!=Sample2])
    
    # 缂佸嫬缂撻弬鎵畱濞村鐦梿鍡礄閸樼喖娉﹂崥?+閺傛壆娈戝ù瀣槸闂嗗棗鎮庨敍?
    Sample5 <- c(CDmean_tst, Sample4)
    
    Sample5 <- Reference[Sample5]   # 瀵版鍩屽ù瀣槸闂嗗棗鎮庨惃鍑D閸??
    
    ## 瀵版鍩孏matrix, 閸楃牰inship閻晠妯€
    GmAt <- Gmatrix[Sample5, Sample5]
    
    ## 閻晠妯€閸欐牠鈧??
    invAt1 <- solve(GmAt)
    
    ## 鐠侊紕鐣诲В蹇庣鏉烆喚娈慍D閸у洤鈧??
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
  
  
  plot(CDmeanExchange) # CDMean閸婅偐鈥樻穱婵囨箒鐡掑啿顧勯惃鍕嚡娴??
  
  dev.off()
  
  SampleOptimiz <- trn_new # 娴兼ê瀵查梿鍡楁値
  
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
    
    trn_tst <- c(TSTsample, x)   # 寤烘ā闆嗗悎澧炲姞涓€涓潗鏂?
    
    trn_tst <- trn_tst[order(trn_tst)]   # 鎺掑簭
    
    Reference <- subset(Reference, GID %in% trn_tst)   # 鎻愬彇鍙傝€冧腑鐨凣ID
    
    Reference <- Reference[order(Reference$GID), ,drop =F]   # 鎻愬彇鎴愭暟鎹
    
    
    ### 鑾峰緱寤烘ā闆嗗悎鐨勫叧绯荤煩闃?
    Gmatrix <- Gmatrix[trn_tst, trn_tst]
    
    ## 澧炲姞Ref鐨勫垪锛屽師寤烘ā缇や綋1锛屽叾浠栦负2
    Reference$Value <- with(Reference, ifelse(GID %in%TSTsample , 1,2))   # 鏍囨敞1锛?2
    refm <- data.frame(GID= Reference$GID, GRM=Reference$Value)
    refm$GID <- as.character(refm$GID)   # GID鍙樻垚瀛楃涓?
    ## 灏嗙煩闃佃浆涓轰袱涓ょ粍鍚堢殑褰㈠紡 
    GmatrixL <- reshape2::melt(as.matrix(Gmatrix))   # reshape
    GmatrixL$Var1 <- as.character(GmatrixL$Var1)   # 杞垚瀛楃涓?
    GmatrixL$Var2 <- as.character(GmatrixL$Var2)   # 杞垚瀛楃涓?
    names(GmatrixL) <- c('rowL', 'colL', 'GRM')   # 澧炲姞鍚嶇О
    GmatrixLR <- merge(refm, GmatrixL, by.x = 'GID', by.y = 'rowL', all = TRUE)   # 鍖归厤绗竴鍒楋紝鐩稿綋浜庢祴璇曠兢浣?
    GmatrixLR$GID <- as.character(GmatrixLR$GID)   # GID杞垚瀛楃鍨?
    GmatrixLRC <- merge(refm, GmatrixLR, by.x = 'GID', by.y = 'colL', all = TRUE)   # 鍖归厤绗簩鍒楋紝鐩稿綋浜庡缓妯＄兢浣?
    GmatrixLRC$GID <- as.character(GmatrixLRC$GID)   # GID杞崲鎴愬瓧绗﹀瀷
    names(GmatrixLRC) <- c('TST', 'TSTpop', 'TRN', 'TRNpop', 'GRM')   # 缁欐爣棰?
    ##Subset only the pairwise combination of the Testing set(TSTpop) and the Training set
    tmp <- subset(GmatrixLRC, TSTpop == 1 & TRN%in%x)   # 鎻愬彇娴嬭瘯缇や綋涓?1锛屽缓妯＄兢浣撲负2鐨勩€?
    
    AVG <- mean(tmp$GRM)   # GRM鐨勫潎鍊?
    
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
  sommerPheno <- sommerPheno[c(tst,Result[,1]),]
  row.names(sommerPheno) <- sommerPheno[,taxa_name_in_Pheno]
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
      if (is.unsorted(K))
      K <- K[sort(c(vv,Result[,1])),sort(c(vv,Result[,1]))]
    }else if (!is.null(group)&!is.null(train_pop)&!is.null(test_pop)){
      # A pre B
      # 建模群体抽样
      vv <- sommerPheno[sommerPheno$group==test_pop,taxa_name_in_Pheno]   # test_pop as a testing population
      # 预测群体
      y.trn <- sommerPheno
      row.names(y.trn) <- y.trn$id
      y.trn$id <- as.factor(as.character(y.trn$id))
      y.trn[vv,"X1"] <- NA
      head(y.trn)
    }else if(is.null(Env)){
      # NoGbyE
      # 建模群体抽样
      vv <- sample(rownames(sommerPheno),round(nrow(sommerPheno)/K_fold))   # 抽取20%作为验证群体（被预测）
      # 预测群体
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
      y.trn2 <- sommerPheno   # 多环境数据
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
      ans_rrBLUPr <- mixed.solve(Pheno_train, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
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