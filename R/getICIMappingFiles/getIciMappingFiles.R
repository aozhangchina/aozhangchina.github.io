# Get a snp file for WinIciMapping from a HMP file
# Made by Zhang Ao
# Updated 20220401 18:12:34
if(!exists("f_name")){f_name <- "my_icimapping_data"}
if(!exists("first_change")){first_change <- FALSE}
if(!exists("Eliminate_heterozygous_bases")){Eliminate_heterozygous_bases <- TRUE}
if(!exists("pheno_colname")){pheno_colname <- BLUE}
source("https://dataholdcn.cn/R/GBIT/GBIT.R")
# The part of import file.
GBIT.setwd(choose.dir())   # Set a directory for work
myGeno <- GBIT.readFile(choose = T,header = T)
myPheno <- GBIT.readFile(choose = T,header = T)

# Match Pheno and Geno
if (exists("myPheno")){
    # View(myPheno)
    myPheno2 <- myPheno[c("Genotype",pheno_colname)]
    myPheno2 <- myPheno2[-which(is.na(myPheno2[,1])),]

    geno_part1 <- myGeno[,1:11]
    geno_part2 <- myGeno[,12:dim(myGeno)[2]]
    if (!identical(names(geno_part2),myPheno2[,1])){
        all_material <- intersect(names(geno_part2),myPheno2[,1])
        which(myPheno2[,1] %in% all_material)
        myPheno2 <- myPheno2[which(myPheno2[,1] %in% all_material),]
        geno_part2 <- geno_part2[,which(names(geno_part2) %in% all_material)]
        pheno_geno_index <- match(names(geno_part2),myPheno2[,1])
        myPheno2 <- myPheno2[pheno_geno_index,]
        myGeno <- cbind(geno_part1,geno_part2)
    }
}

# The part for snp file.
## Recode marker names
taxa <- myGeno[,1]
taxa <- cbind(taxa,paste0("M",1:length(taxa)))
write.table(taxa, "GBIT/marker_rename.txt", sep = "\t", col.names = F, row.names = F)

# myGeno[1:15, 1:15]
## snp_GeneralInfo
dim(myGeno)
base_information <- c(4, 1, 2, dim(myGeno)[1], dim(myGeno)[2]-11)
write.table(base_information, "GBIT/snp_GeneralInfo.txt", sep = "\t", col.names = F, row.names = F)
write.table(base_information, paste0("GBIT/",f_name,".snp"), sep = "\t", col.names = F, row.names = F)
write.table("", paste0("GBIT/",f_name,".snp"), sep = "\t", col.names = F, row.names = F, append = T, quote = F)
base_parents <- myGeno$alleles
getParents <- function(my_vector){
    # my_vector <- unlist(base_parents2[1])
    my_vector <- strrep(as.character(my_vector), 2)
    return(my_vector[1:2])
}
base_parents2 <- strsplit(myGeno$alleles, split = "/")
if (first_change==TRUE){
        worng_index <- which(lengths(base_parents2)>2)
    if (length(worng_index)>0){
        IUPAC_matrix <- matrix(c("AG","CT","CG","AT","GT","AC","R","Y","S","W","K","M"), ncol = 2)
        change_letter <- c("CG","CT","GT")
            for (i in worng_index){
                table_1 <- names(table(as.character((myGeno[i,12:dim(myGeno)[2]]))))
                correct_names <- table_1[1:2]
                if (correct_names %in% change_letter){
                    correct_names <- paste0(sort(unlist(strsplit(correct_names,split="")),decreasing=T),collapse="")
                }
                correct_names <- paste0(sort(correct_names), collapse = "")
                index <- !myGeno[i,12:dim(myGeno)[2]] %in% c(as.character(table_1[1:2]),IUPAC_matrix[which(IUPAC_matrix[,1]==correct_names),][2],"N")
                myGeno[i,12:dim(myGeno)[2]][index] <- "N"
            }
    }

    base_parents3 <- t(sapply(base_parents2, getParents))
}else{
    base_parents3 <- t(sapply(base_parents2, getParents))
}

base_parents3 <- cbind(taxa[,2] , base_parents3)
# result_list <- list(matrix(base_information,ncol=1),base_parents3)
write.table(base_parents3,"GBIT/snp_Parents.txt",sep="\t",col.names =F, row.names = F, quote = F)
write.table(base_parents3, paste0("GBIT/",f_name,".snp"), sep = "\t", col.names = F, row.names = F, append = T, quote = F)
write.table("", paste0("GBIT/",f_name,".snp"), sep = "\t", col.names = F, row.names = F, append = T, quote = F)
# capture.output(result_list,file="filename.txt")

base_genotype <- myGeno[,12:dim(myGeno)[2]]
# base_genotype[1:15,1:15]
getGenotypes <- function(my_vector){
    IUPAC2_matrix <- matrix(c("AG","CT","CG","AT","GT","AC","RR","YY","SS","WW","KK","MM"), ncol = 2)
    
    # my_vector <- unlist(base_parents2[1])
    # my_vector <- base_genotype[1314,]
    table2 <- names(sort(table(as.character(my_vector)),decreasing = T))
    table2 <- table2[which(table2 %in% c("A","T","C","G"))]
    main_names <- paste0(sort(table2[1:2]),collapse = "")
    heb_names <- IUPAC2_matrix[which(IUPAC2_matrix[,1]==main_names),2]
    my_vector <- strrep(as.character(my_vector), 2)
    # table(my_vector)
    # unique(my_vector)
                #####$ Heterozygosity is not allowed
                    change_letter <- c("CG","CT","GT")
                    if(main_names %in% change_letter){
                        main_names <- paste0(sort(unlist(strsplit(main_names,split="")),decreasing=T),collapse="")
                        index <- which(my_vector==heb_names)
                        my_vector[index] <- main_names
                        # index <- my_vector %in% c(main_names)
                        # my_vector[index] <- "NN"
                    }
                    index <- !my_vector %in% c(names(sort(table(my_vector),decreasing = T))[1:2],main_names)
                    my_vector[index] <- "NN"
                        if (Eliminate_heterozygous_bases==TRUE){
                            index <- !my_vector %in% c("AA","TT","CC","GG","NN")
                            my_vector[index] <- "NN"
                        }
                


    return(my_vector)
}
base_genotype2 <- apply(base_genotype,1,getGenotypes)
base_genotype2 <- t(base_genotype2)
# View(base_genotype)
# base_genotype[1:15,1:15]
dim(base_genotype2)
dim(myGeno)
base_genotype2 <- cbind(taxa[,2],base_genotype2)
write.table(base_genotype2,"GBIT/snp_Genotype.txt",sep="\t",col.names =F, row.names = F, quote = F)
write.table(base_genotype2, paste0("GBIT/",f_name,".snp"), sep = "\t", col.names = F, row.names = F, append = T, quote = F)
write.table("", paste0("GBIT/",f_name,".snp"), sep = "\t", col.names = F, row.names = F, append = T, quote = F)
base_anchor <- cbind(taxa[,2],myGeno[,3])
write.table(base_anchor,"GBIT/snp_Anchor.txt",sep="\t",col.names =F, row.names = F, quote = F)
write.table(base_anchor, paste0("GBIT/",f_name,".snp"), sep = "\t", col.names = F, row.names = F, append = T, quote = F)
write.table("", paste0("GBIT/",f_name,".snp"), sep = "\t", col.names = F, row.names = F, append = T, quote = F)

# The part for bip file.
## bip_GeneralInfo
if (exists("myPheno")){
    trait_number <- length(pheno_colname)
    bip_base_information <- c(1, 4, 1, 2, 1, length(unique(myGeno$chrom)), dim(myGeno)[2]-11, trait_number)
    write.table(bip_base_information, "GBIT/bip_GeneralInfo.txt", sep = "\t", col.names = F, row.names = F)
    write.table(bip_base_information, paste0("GBIT/",f_name,".bip"), sep = "\t", col.names = F, row.names = F)
    write.table("", paste0("GBIT/",f_name,".bip"), sep = "\t", col.names = F, row.names = F, append = T, quote = F)

    ch <- 1:length(unique(myGeno$chrom))
    ch <- matrix(paste0("Ch",ch, sep = ""),ncol=1)
    ch <- cbind(ch,"number"=NA)
    for (i in 1:length(unique(myGeno$chrom))){
        ch[i,2] <- length(which(myGeno$chrom==i))
    }
    write.table(ch, "GBIT/bip_Chromosome.txt", sep = "\t", col.names = F, row.names = F, quote = F)
    write.table(ch, paste0("GBIT/",f_name,".bip"), sep = "\t", col.names = F, row.names = F, append =T, quote = F)

    bip_base_LinkageMap <- cbind(taxa[,2],myGeno[,c("chrom","pos")])
    bip_base_LinkageMap$pos <- bip_base_LinkageMap$pos/100000
    write.table(bip_base_LinkageMap, "GBIT/bip_LinkageMap.txt", sep = "\t", col.names = F, row.names = F, quote = F)
    write.table(bip_base_LinkageMap, paste0("GBIT/",f_name,".bip"), sep = "\t", col.names = F, row.names = F, quote = F, append = T)
    write.table("", paste0("GBIT/",f_name,".bip"), sep = "\t", col.names = F, row.names = F, append = T, quote = F)

    bip_base_genotype <- base_genotype2
    for (i in 1:dim(bip_base_genotype)[1]){
        bip_base_genotype[i,] <- gsub("NN",-1,bip_base_genotype[i,])
        bip_base_genotype[i,] <- gsub(base_parents3[i,2],2,bip_base_genotype[i,])
        bip_base_genotype[i,] <- gsub(base_parents3[i,3],0,bip_base_genotype[i,])
    }
    write.table(bip_base_genotype,"GBIT/bip_Genotype.txt",sep="\t",col.names =F, row.names = F, quote = F)
    write.table(bip_base_genotype, paste0("GBIT/",f_name,".bip"), sep = "\t", col.names = F, row.names = F, append = T, quote = F)
    write.table("", paste0("GBIT/",f_name,".bip"), sep = "\t", col.names = F, row.names = F, append = T, quote = F)

    bip_base_Phenotype <- round(t(myPheno2[,2]),3)
    row.names(bip_base_Phenotype) <- pheno_colname
    write.table(bip_base_Phenotype,"GBIT/bip_Phenotype.txt",sep="\t",col.names =F, row.names = T, quote = F)
    write.table(bip_base_Phenotype, paste0("GBIT/",f_name,".bip"), sep = "\t", col.names = F, row.names = T, append = T, quote = F)
    write.table("", paste0("GBIT/",f_name,".bip"), sep = "\t", col.names = F, row.names = F, append = T, quote = F)
}
