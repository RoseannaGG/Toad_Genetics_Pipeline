#look at relatedness
libary(ggplot2)
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
#create histogram of relatedness values
  histprefilt <- ggplot(randreldf[[nam]], aes(x=value))+geom_histogram(color="black", fill="white")
  ggsave(file=paste(outputwd, "/", nam, "_presibfilt_hist.pdf", sep=""), histprefilt, width=10, height=8, units="in")
  sibdf[[nam]] <- randreldf[[nam]][which(randreldf[[nam]]$value>=relthreshold),]
#get rough idea of size of sibling groups
#base this on the loose idea that the number of times an individual is in Var1 represents the number of siblings
#but that might not be quite how the matrix works, and possibility that member of sibling group in there twice
#i.e. ToadA is siblings to ToadB and ToadC, and both ToadA and ToadC are in this plot
#COLONY would be better to look at this, but that will be up to time at the moment
  famsize <- aggregate(sibdf[[nam]]$Var2, by=list(sibdf[[nam]]$Var1), length)
  write.csv(famsize, file=paste(outputwd, "/", nam, "_sibling_groups.csv", sep=""))
  miss <- read.table(paste(outwd, files[grepl(gsub("_relatedness", files[x]), files)], ".imiss", sep=""))
  sibmissdf[[nam]] <- merge(sibdf[[nam]], miss, by.x="Var1", by.y="ID")
  sibmissdf[[nam]] <- merge(sibmissdf[[nam]], miss, by.x="Var2", by.y="ID")
  sibmissdf[[nam]]$Cut <- submissdf[[nam]]$Var1
  submissdf[[nam]][which(submissdf[[nam]]$miss.x<subsmidddf[[nam]]$miss.y),]$Cut <- submissdf[[nam]]$Var2
  test <- randreldf[[nam]][-which(randreldf[[nam]]$Var1 %in% submissdf[[nam]]$Cut | randreldf[[nam]]$Var2 %in% submissdf[[nam]]$Cut),]
  #create histogram of relatedness values after filtering
  histpostfilt <- ggplot(test, aes(x=value))+geom_histogram(color="black", fill="white")
  ggsave(file=paste(outputwd, "/", nam, "_postsibfilt_hist.pdf", sep=""), histpostfilt, width=10, height=8, units="in")
  chk[[nam]] <- data.frame("Pond"=nam, "Max.Relatedness"=max(test$value), "Num.Ind"=length(unique(c(test$Var1, test$Var2))))
  indkeep[[nam]] <- data.frame("Ind"=unique(c(test$Var1, test$Var2)))
  write.table(indkeep[[nam]], file=paste(outputwd, "/", nam, "_NoSibs.txt", sep=""), row.names=F, col.names = F, quote=F, sep="\t")
}

chkfilt <- do.call("rbind", chk)

histprefilt <- ggplot(randreldf[[nam]], aes(x=value))+geom_histogram(color="black", fill="white")
ggsave(file=paste(outputwd, "/", nam, "_presibfilt_hist.png", sep=""), histprefilt, width=10, height=8, units="in")

## my code

# make dataframe
sibdataframe_766_3356<-randreldf[[nam]]
str(sibdataframe_766_3356)
dim(sibdataframe_766_3356) # 293761      3

# delete match of individuals
sibdataframe_766_3356_onlynomatchindiv<-subset(sibdataframe_766_3356, Var1!=Var2)
dim(sibdataframe_766_3356_onlynomatchindiv) #292995      3

# subset to just ponds
sub_sibdataframe_766_3356_onlynomatchindiv = sibdataframe_766_3356_onlynomatchindiv[substr(sibdataframe_766_3356_onlynomatchindiv$Var1,1,9)==substr(sibdataframe_766_3356_onlynomatchindiv$Var2,1,9),]
dim(sub_sibdataframe_766_3356_onlynomatchindiv) #7240    3
head(sub_sibdataframe_766_3356_onlynomatchindiv)

# histprefilt_withinponds
histprefilt_withinponds<-ggplot(sub_sibdataframe_766_3356_onlynomatchindiv, aes(x=value))+geom_histogram(color="black", fill="white")
dev.off()
ggsave(file=paste(outputwd, "/", nam, "_presibfilt_WITHINPONDS_hist.png", sep=""), histprefilt_withinponds, width=10, height=8, units="in")

# basic R plotting
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/sibHist_wihtinponds_766INDIV_3356SNPS.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(sub_sibdataframe_766_3356_onlynomatchindiv$value,breaks=seq(0,1.3,l=50), main="Relatedness within ponds only",
    # ylim = c(0, 1500),
     xlab="Relatedness - higher => more related (over 0.4 -  full sibs or more) ",
     ylab="Frequency")
dev.off()

# one column
sub_sibdataframe_766_3356_onlynomatchindiv_onecol<-sub_sibdataframe_766_3356_onlynomatchindiv[,2:3]


#### MY CODE

# relatedness over 0.35 
sibrelover035<-sibdf[[nam]] 
dim(sibrelover035) #24899     3


# delete match of individuals
sibrelover035_onlynomatchindiv<-subset(sibrelover035, Var1!=Var2)
dim(sibrelover035_onlynomatchindiv) #24133     3

### just within ponds
sub_sibrelover035_onlynomatchindiv = sibrelover035_onlynomatchindiv[substr(sibrelover035_onlynomatchindiv$Var1,1,9)==substr(sibrelover035_onlynomatchindiv$Var2,1,9),]
dim(sub_sibrelover035_onlynomatchindiv) #3645    3


head(sub_sibrelover035_onlynomatchindiv)

### just isolate one row of names
sub_sibrelover035_names<-sub_sibrelover035_onlynomatchindiv$Var1
sub_sibrelover035_names

str(sub_sibrelover035_names)
head(sub_sibrelover035_names)

sub_sibrelover035_names.char<-as.character(sub_sibrelover035_names)
str(sub_sibrelover035_names.char) #chr [1:3645]


#deletes the repeats --- gives list of individuals to remove that I can feed back into vcf tools
sub_sibrelover035_names_unique<-unique(sub_sibrelover035_names)
str(sub_sibrelover035_names_unique) 

sub_sibrelover035_names_unique.char<-as.character(sub_sibrelover035_names_unique)
str(sub_sibrelover035_names_unique.char) #chr [1:477] 


head(sub_sibrelover035_names_unique.char)

# HIST
histprefilt2 <- ggplot(sub_sibrelover035_onlynomatchindiv, aes(x=value))+geom_histogram(color="black", fill="white")
ggsave(file=paste(outputwd, "/", nam, "_take035_presibfilt_hist.png", sep=""), histprefilt2, width=10, height=8, units="in")





