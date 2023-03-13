
rm(list=ls())


#look at relatedness
relthreshold <- 0.35
allout <- list.files("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/plink/766INDIV_3356SNPS/")
outwd <- "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/plink/766INDIV_3356SNPS/"
outputwd <-"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/plink/outputdirectory/"
relfiles <- allout[grepl(".rel", allout)]
files <- unique(gsub("\\..*", "", relfiles))
library(reshape2)
library(tidyr)
mat <- list()
ids <- list()
relmat <- list()
randreldf <- list()
sibdf <- list()
sibmissdf <- list()
chk <- list()
indkeep <- list()
for(x in 1:length(files)){
  nam <- files[x]
  mat[[nam]] <- read.table(paste(outwd, files[x], ".rel", sep=""))
  ids[[nam]] <- read.table(paste(outwd, files[x], ".rel.id", sep=""))
  relmat[[nam]] <- as.matrix(mat[[nam]])
  colnames(relmat[[nam]]) <- ids[[nam]]$V2
  row.names(relmat[[nam]]) <- ids[[nam]]$V2
  relmat[[nam]][upper.tri(relmat[[nam]])] <- NA
  randreldf[[nam]] <- reshape2::melt(relmat[[nam]])
  randreldf[[nam]] <- randreldf[[nam]][-which(is.na(randreldf[[nam]]$value)),]
  sibdf[[nam]] <- randreldf[[nam]][which(randreldf[[nam]]$value>=relthreshold),]
  miss <- read.table(paste(outwd, files[grepl(gsub("_relatedness", files[x]), files)], ".imiss", sep=""))
  sibmissdf[[nam]] <- merge(sibdf[[nam]], miss, by.x="Var1", by.y="ID")
  sibmissdf[[nam]] <- merge(sibmissdf[[nam]], miss, by.x="Var2", by.y="ID")
  sibmissdf[[nam]]$Cut <- submissdf[[nam]]$Var1
  submissdf[[nam]][which(submissdf[[nam]]$miss.x<subsmidddf[[nam]]$miss.y),]$Cut <- submissdf[[nam]]$Var2
  test <- randreldf[[nam]][-which(randreldf[[nam]]$Var1 %in% submissdf[[nam]]$Cut | randreldf[[nam]]$Var2 %in% submissdf[[nam]]$Cut),]
  chk[[nam]] <- data.frame("Pond"=nam, "Max.Relatedness"=max(test$value), "Num.Ind"=length(unique(c(test$Var1, test$Var2))))
  indkeep[[nam]] <- data.frame("Ind"=unique(c(test$Var1, test$Var2)))
  write.table(indkeep[[nam]], file=paste(outputwd, "/", nam, "_NoSibs.txt", sep=""), row.names=F, col.names = F, quote=F, sep="\t")
}


# relatedness over 0.35 
sibrelover035<-sibdf[[nam]] 
dim(sibrelover035) #24899     3


# delete match of individuals
sibrelover035_onlynomatchindiv<-subset(sibrelover035, Var1!=Var2)
dim(sibrelover035_onlynomatchindiv) #24133     3

### just within ponds
sub_sibrelover035_onlynomatchindiv = sibrelover035_onlynomatchindiv[substr(sibrelover035_onlynomatchindiv$Var1,1,9)==substr(sibrelover035_onlynomatchindiv$Var2,1,9),]
dim(sub_sibrelover035_onlynomatchindiv) #3645    3


### just isolate one row of names
sub_sibrelover035_names<-sub_sibrelover035_onlynomatchindiv$Var1
sub_sibrelover035_names

str(sub_sibrelover035_names)

sub_sibrelover035_names.char<-as.character(sub_sibrelover035_names)
str(sub_sibrelover035_names.char) #chr [1:3645]


#deletes the repeats --- gives list of individuals to remove that I can feed back into R to remove
sub_sibrelover035_names_unique<-unique(sub_sibrelover035_names)
str(sub_sibrelover035_names_unique) 

sub_sibrelover035_names_unique.char<-as.character(sub_sibrelover035_names_unique)
str(sub_sibrelover035_names_unique.char) #chr [1:477] 
(477/766)*100 #62.27154
766-477 # 289



### subset data that is being kept
pop.data_289<-pop.data[which(!pop.data$sample.id %in% sub_sibrelover035_names_unique.char),]
dim(pop.data_289) # 289   7
levels(pop.data_289$pop)

#check number of ponds
pop.data_289_2<-pop.data_289
pop.data_289_2$pop<-as.factor(pop.data_289_2$pop)

levels(pop.data_289_2$pop) # 44

# lost these populations
# "R01-SI-LL"
# "R02-SC-HA"
# "R02-SS-KH"









####### 0.4
relthreshold <- 0.4
allout <- list.files("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/plink/766INDIV_3356SNPS/")
outwd <- "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/plink/766INDIV_3356SNPS/"
outputwd <-"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/plink/outputdirectory/"
relfiles <- allout[grepl(".rel", allout)]
files <- unique(gsub("\\..*", "", relfiles))
library(reshape2)
library(tidyr)
mat <- list()
ids <- list()
relmat <- list()
randreldf <- list()
sibdf <- list()
sibmissdf <- list()
chk <- list()
indkeep <- list()
for(x in 1:length(files)){
  nam <- files[x]
  mat[[nam]] <- read.table(paste(outwd, files[x], ".rel", sep=""))
  ids[[nam]] <- read.table(paste(outwd, files[x], ".rel.id", sep=""))
  relmat[[nam]] <- as.matrix(mat[[nam]])
  colnames(relmat[[nam]]) <- ids[[nam]]$V2
  row.names(relmat[[nam]]) <- ids[[nam]]$V2
  relmat[[nam]][upper.tri(relmat[[nam]])] <- NA
  randreldf[[nam]] <- reshape2::melt(relmat[[nam]])
  randreldf[[nam]] <- randreldf[[nam]][-which(is.na(randreldf[[nam]]$value)),]
  sibdf[[nam]] <- randreldf[[nam]][which(randreldf[[nam]]$value>=relthreshold),]
  miss <- read.table(paste(outwd, files[grepl(gsub("_relatedness", files[x]), files)], ".imiss", sep=""))
  sibmissdf[[nam]] <- merge(sibdf[[nam]], miss, by.x="Var1", by.y="ID")
  sibmissdf[[nam]] <- merge(sibmissdf[[nam]], miss, by.x="Var2", by.y="ID")
  sibmissdf[[nam]]$Cut <- submissdf[[nam]]$Var1
  submissdf[[nam]][which(submissdf[[nam]]$miss.x<subsmidddf[[nam]]$miss.y),]$Cut <- submissdf[[nam]]$Var2
  test <- randreldf[[nam]][-which(randreldf[[nam]]$Var1 %in% submissdf[[nam]]$Cut | randreldf[[nam]]$Var2 %in% submissdf[[nam]]$Cut),]
  chk[[nam]] <- data.frame("Pond"=nam, "Max.Relatedness"=max(test$value), "Num.Ind"=length(unique(c(test$Var1, test$Var2))))
  indkeep[[nam]] <- data.frame("Ind"=unique(c(test$Var1, test$Var2)))
  write.table(indkeep[[nam]], file=paste(outputwd, "/", nam, "_NoSibs.txt", sep=""), row.names=F, col.names = F, quote=F, sep="\t")
}



  dim(  sibdf[[nam]]) #1809
  
  
  # relatedness over 0.4 
  sibrelover04<-sibdf[[nam]] 
  dim(sibrelover04) #1809
  
  
  # delete match of individuals
  sibrelover04_onlynomatchindiv<-subset(sibrelover04, Var1!=Var2)
dim(sibrelover04_onlynomatchindiv)  #1031

  ### just within ponds
  sub_sibrelover04_onlynomatchindiv = sibrelover04_onlynomatchindiv[substr(sibrelover04_onlynomatchindiv$Var1,1,9)==substr(sibrelover04_onlynomatchindiv$Var2,1,9),]
  dim(sub_sibrelover04_onlynomatchindiv) #1031
  
  
  ### just isolate one row of names
  sub_sibrelover04_names<-sub_sibrelover04_onlynomatchindiv$Var1
  sub_sibrelover04_names
  
  str(sub_sibrelover04_names)
  
  sub_sibrelover04_names.char<-as.character(sub_sibrelover04_names)
  str(sub_sibrelover04_names.char) #chr [1:1031]
  
  #deletes the repeats --- gives list of individuals to remove that I can feed back into into R to remove
  sub_sibrelover04_names_unique<-unique(sub_sibrelover04_names)
  str(sub_sibrelover04_names_unique) 
  
  sub_sibrelover04_names_unique.char<-as.character(sub_sibrelover04_names_unique)
  str(sub_sibrelover04_names_unique.char) #chr [1:256]
  
  (256/766)*100 # 33.42037
  766-256 # 510
  
  dim(pop.data.766)
  
  ### subset data that is being kept
  pop.data_510<-pop.data.766[which(!pop.data.766$sample.id %in% sub_sibrelover04_names_unique.char),]
  dim(pop.data_510) # 510   7

  write.table(sub_sibrelover04_names_unique.char, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/766INDIV_3356SNPS_listindiv_plink_sub_sib_rel_over_04_names_unique.char.txt", row.names = FALSE, col.names = FALSE,quote=FALSE )
  
  
  
  