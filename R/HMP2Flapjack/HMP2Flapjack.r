# Some information for users.
cat("Please select a file in HMP format. This tool made by Zhang Ao. 6/24/2020 https://datahold.cn \n")

# Choose a file and get some information from it.
fBaseInfo <- choose.files()
fBaseName <- basename(fBaseInfo)
fDirName <- dirname(fBaseInfo)

# Set workspace
setwd(fDirName)

# load a package.
if(!"data.table" %in% installed.packages()){
    install.packages("data.table")
    library(data.table)
}else{
    library(data.table)
}

# Read a genotypic file.
geno <- fread(fBaseInfo)
geno <- as.data.frame(geno)

# Get header of files for MAP and GENOTYPE.
fjFile_map_first <- c("# fjFile = MAP",NA,NA)
fjFile_geno_first <- "# fjFile = GENOTYPE"

# Get names of all of chromosomes.
all_chromosome <- unique(geno$chrom)

# delet NA chromosome.
if (any(is.na(all_chromosome))){
    all_chromosome <- all_chromosome[-which(is.na(all_chromosome))]
    cat("There is the NA chromosome, it has been deleted. \n")
}

# delet NA pos.
geno <- geno[-which(is.na(geno$chrom)),]

# modify the unit of positions.
geno$pos <- round(geno$pos/1e6,1)

MAP_temp <- NULL   # initialization of MAP
# the max of position in each chromosome as a chromosome length.
for (i in 1:length(all_chromosome)){
    chromosome_length <- max(geno$pos[geno$chrom==all_chromosome[i]],na.rm=TRUE)
    MAP_tmp <- geno[geno$chrom==all_chromosome[i],c("rs#","chrom","pos")]
    chromosome_length_first <- c(all_chromosome[i],chromosome_length,NA)
    MAP_tmp <- rbind(chromosome_length_first,MAP_tmp)
    MAP_temp <- rbind(MAP_temp,MAP_tmp)
}
MAP <- rbind(fjFile_map_first,MAP_temp)

# Output MAP file for flapjack.
write.table(MAP,paste0(fBaseName,"_MAP.txt"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,na="")

# delet some colume from geno.
useless_column <- which(names(geno) %in% c("alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panelLSID","QCcode"))
geno <- geno[,-useless_column]

# Transpose geno
geno_t <- t(geno)
rownames_geno_t <- row.names(geno_t)
names_geno_t[1] <- NA
colnames_geno_t <- geno_t[1,]
geno_t <- geno_t[-1,]

geno_t <- sub("N","-",geno_t)
geno_t <- sub("NA","-",geno_t)

# change R,Y,S,W,K,M to A/G,C/T,C/G,A/T,G/T,A/C.
geno_t <- sub("R","A/G",geno_t)
geno_t <- sub("Y","C/T",geno_t)
geno_t <- sub("S","C/G",geno_t)
geno_t <- sub("W","A/T",geno_t)
geno_t <- sub("K","G/T",geno_t)
geno_t <- sub("M","A/C",geno_t)

geno_t <- rbind(colnames_geno_t,geno_t)
geno_t <- cbind(rownames_geno_t,geno_t)
geno_t[1,1] <- NA
fjFile_geno_first <- c(fjFile_geno_first,rep(NA,ncol(geno_t)-1))
GENOTYPE <- rbind(fjFile_geno_first,geno_t)

# Output GENOTYPE file for flapjack.
write.table(GENOTYPE,paste0(fBaseName,"_GENOTYPE.txt"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,na="")