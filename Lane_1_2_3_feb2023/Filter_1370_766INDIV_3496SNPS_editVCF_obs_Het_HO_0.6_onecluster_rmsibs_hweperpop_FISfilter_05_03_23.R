rm(list=ls())

library(vcfR)
library(adegenet)
library(hierfstat)
library(dartR)
library(ggplot2)
library(spaa)
library(fields)
library(vegan)
library(StAMPP)
library(ggmap)
library(car)
library(VennDiagram)
library(radiator)
library(pegas)
library(RColorBrewer)
library(SNPRelate)
library(reshape2)
library(svglite)


setwd("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/")

sink("R_ref_Aboreas_PCA_RGGoutput_766INDIV_3496SNP.txt", split = TRUE)

####################### LD Filtering - ONLY RESULTING IN LOSING 1 SNP SO DIDN'T BOTHER REMOVING ################


vcfF22 <- "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS.vcf"
snpgdsVCF2GDS(vcfF22, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS.gds", method="biallelic.only")

snpgdsSummary("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS.gds")
""

snpgdsClose("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS.gds")

genofile <- SNPRelate::snpgdsOpen("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS.gds")

# LD Pruning
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, slide.max.bp = 100000, autosome.only = FALSE, ld.threshold = 0.2)

#3,390 markers are selected in total.


#increase to 500,000 - composite and others!!
set.seed(1000)
snpset1 <- snpgdsLDpruning(genofile, slide.max.bp = 500000, method="corr", autosome.only = FALSE, ld.threshold = 0.2)

# 3,213 markers are selected in total.


## try method = r 500k
set.seed(1000)
snpset2 <- snpgdsLDpruning(genofile, slide.max.bp = 500000, autosome.only = FALSE, method="r",ld.threshold = 0.2)

# 3,270 markers are selected in total.


#dprime 500k
set.seed(1000)
snpset3 <- snpgdsLDpruning(genofile, slide.max.bp = 500000, autosome.only = FALSE, method="dprime",ld.threshold = 0.2)

# 2,747 markers are selected in total.





#increase to 500,000 - composite and others!!
set.seed(1000)
snpset1 <- snpgdsLDpruning(genofile, slide.max.bp = 500000, method="corr", autosome.only = FALSE, ld.threshold = 0.2)

# 3,213 markers are selected in total. (from 3496) = removed 283 SNPs in high LD

str(snpset1)


"List of 65
 $ chrScaffold_1__2_contigs__length_508792519 : int [1:342] 1 2 3 5 6 7 8 9 10 11 ...
 $ chrScaffold_2__3_contigs__length_804548876 : int [1:519] 373 374 375 376 377 378 379 380 381 382 ...
 $ chrScaffold_3__2_contigs__length_655571293 : int [1:430] 934 935 936 937 938 939 940 941 942 943 ...
 $ chrScaffold_4__2_contigs__length_776011713 : int [1:539] 1398 1399 1400 1401 1402 1403 1404 1405 1407 1408 ...
 $ chrScaffold_5__2_contigs__length_593928783 : int [1:416] 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 ...
 $ chrScaffold_6__2_contigs__length_328854778 : int [1:230] 2437 2438 2439 2440 2441 2442 2443 2444 2445 2446 ...
 $ chrScaffold_7__1_contigs__length_206032987 : int [1:184] 2689 2690 2691 2694 2695 2696 2697 2698 2699 2700 ...
 $ chrScaffold_8__1_contigs__length_198584680 : int [1:159] 2897 2898 2899 2900 2901 2903 2904 2905 2906 2907 ...
 $ chrScaffold_9__1_contigs__length_196030596 : int [1:180] 3065 3066 3067 3068 3069 3071 3072 3073 3074 3075 ...
 $ chrScaffold_10__1_contigs__length_179080748: int [1:154] 3260 3262 3263 3264 3266 3267 3268 3269 3270 3271 ...
 $ chrScaffold_15__1_contigs__length_349370   : int 3435
 $ chrScaffold_20__1_contigs__length_332182   : int 3436
 $ chrScaffold_23__1_contigs__length_311584   : int 3437
 $ chrScaffold_33__1_contigs__length_253271   : int 3438
 $ chrScaffold_35__1_contigs__length_248379   : int 3439
 $ chrScaffold_37__1_contigs__length_241972   : int 3440
 $ chrScaffold_43__1_contigs__length_235388   : int 3441
 $ chrScaffold_53__1_contigs__length_217242   : int 3442
 $ chrScaffold_59__1_contigs__length_206510   : int 3443
 $ chrScaffold_61__1_contigs__length_203537   : int [1:2] 3445 3446
 $ chrScaffold_87__1_contigs__length_164672   : int 3447
 $ chrScaffold_93__1_contigs__length_161056   : int 3448
 $ chrScaffold_103__1_contigs__length_146291  : int 3449
 $ chrScaffold_116__1_contigs__length_134014  : int 3450
 $ chrScaffold_128__1_contigs__length_128476  : int [1:2] 3451 3452
 $ chrScaffold_131__1_contigs__length_127423  : int 3453
 $ chrScaffold_136__1_contigs__length_125746  : int 3454
 $ chrScaffold_141__1_contigs__length_124535  : int 3455
 $ chrScaffold_146__1_contigs__length_121756  : int 3457
 $ chrScaffold_156__1_contigs__length_117762  : int 3458
 $ chrScaffold_162__1_contigs__length_115834  : int 3459
 $ chrScaffold_164__1_contigs__length_115762  : int 3460
 $ chrScaffold_182__1_contigs__length_108084  : int 3461
 $ chrScaffold_211__1_contigs__length_100455  : int 3462
 $ chrScaffold_220__1_contigs__length_98461   : int 3463
 $ chrScaffold_257__1_contigs__length_89402   : int [1:3] 3464 3465 3466
 $ chrScaffold_263__1_contigs__length_87712   : int 3467
 $ chrScaffold_264__1_contigs__length_87568   : int 3468
 $ chrScaffold_286__1_contigs__length_84388   : int 3469
 $ chrScaffold_317__1_contigs__length_78170   : int 3470
 $ chrScaffold_397__1_contigs__length_67266   : int 3471
 $ chrScaffold_543__1_contigs__length_53344   : int 3472
 $ chrScaffold_559__1_contigs__length_52607   : int 3473
 $ chrScaffold_570__1_contigs__length_51877   : int 3474
 $ chrScaffold_638__1_contigs__length_47705   : int 3475
 $ chrScaffold_754__1_contigs__length_41843   : int [1:2] 3476 3477
 $ chrScaffold_835__1_contigs__length_38384   : int 3478
 $ chrScaffold_883__1_contigs__length_36667   : int 3479
 $ chrScaffold_916__1_contigs__length_35273   : int 3480
 $ chrScaffold_933__1_contigs__length_34700   : int 3481
 $ chrScaffold_967__1_contigs__length_33574   : int 3482
 $ chrScaffold_1008__1_contigs__length_32738  : int 3483
 $ chrScaffold_1011__1_contigs__length_32652  : int 3484
 $ chrScaffold_1138__1_contigs__length_29745  : int 3485
 $ chrScaffold_1166__1_contigs__length_29151  : int 3486
 $ chrScaffold_1194__1_contigs__length_28640  : int 3487
 $ chrScaffold_1271__1_contigs__length_27132  : int 3488
 $ chrScaffold_1447__1_contigs__length_24630  : int 3489
 $ chrScaffold_1914__1_contigs__length_19482  : int 3490
 $ chrScaffold_2127__1_contigs__length_17606  : int 3491
 $ chrScaffold_2129__1_contigs__length_17586  : int 3492
 $ chrScaffold_2462__1_contigs__length_15147  : int 3493
 $ chrScaffold_2540__1_contigs__length_14532  : int 3494
 $ chrScaffold_3333__1_contigs__length_9447   : int 3495
 $ chrScaffold_3438__1_contigs__length_8665   : int 3496"



snpset1.df<-as.data.frame(snpset1)

"
Error in (function (..., row.names = NULL, check.rows = FALSE, check.names = TRUE,  : 
                      arguments imply differing number of rows: 342, 519, 430, 539, 416, 230, 184, 159, 180, 154, 1, 2, 3 
                      "

LD.keep.loci.names <- rownames(snpset1.df)

#subset
gl.toad.LD <- gl.toad[ , LD.keep.loci.names]
gl.toad.LD


#######
############ CHANGE LOCI NAMES IN VCF #######

VCFfile<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS.vcf")
head(VCFfile@fix[,'ID'])

VCFfiletest<-VCFfile

#loci names
VCFfiletest@fix[,'ID']

#remove :+ at end of name
VCFfiletest@fix[,'ID']<-gsub("\\:\\+", "", VCFfiletest@fix[,'ID'])
VCFfiletest@fix[,'ID']

#remove :- at end of name
VCFfiletest@fix[,'ID']<-gsub("\\:\\-", "", VCFfiletest@fix[,'ID'])
VCFfiletest@fix[,'ID']

# add ANBO as a prefix
VCFfiletest@fix[,'ID']<-paste("ANBO", sep = '_', VCFfiletest@fix[,'ID'])

VCFfiletest@fix[,'ID']
head(VCFfiletest@fix[,'ID'])

filtered.VCF<-VCFfiletest

"***** Object of Class vcfR *****
766 samples
65 CHROMs
3,496 variants
Object size: 56.3 Mb
0 percent missing data
*****        *****         *****"


vcfR::write.vcf(filtered.VCF, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_editlocinamesinR.vcf.gz")

# then uncompress in folder


####################################
################# MAF ###################
##############################
#loading the data file
#filtered.VCF<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS.vcf")


## keep only fitlered loci and export 
## loci names I want to keep
gl.toad.nosibs.hwe.FIS@loc.names
locinamnes_766_1503<-gl.toad.nosibs.hwe.FIS@loc.names

## column 3 is ID with loci names
filtered.VCF_766_1503<-filtered.VCF[filtered.VCF@fix[,3]%in%locinamnes_766_1503,]

# write
vcfR::write.vcf(filtered.VCF_766_1503, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/419NDIV_1367SNPS_filteredinR.vcf.gz")

# then uncompress in folder



#minor allele freq
mafvcf<-maf(filtered.VCF,element=2)



mafvcf.df<-as.data.frame(mafvcf)

head(mafvcf.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(2,2,2,2))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_766INDIV_3496SNPS.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(mafvcf.df$Frequency,breaks=seq(0,0.5,l=50), main="MAF frequency",
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()


### NB ploting the number of SNPs/loci that have a minor allele freq of x, for every loci across all individuals

113/1082

# shows that this is the number of alleles across all individuals per locus - which cannot be more than the number of individuals*2 (diploid)
# in this case 601 indiv *2 = 1242 is the max number of alleles for a particular locus
summary(mafvcf.df$nAllele)

#major allele freq
mafvcfmajor<-maf(filtered.VCF,element=1)
str(mafvcf)



##### below 0.1

mafvcf.df_below0.1<-mafvcf.df[mafvcf.df$Frequency<=0.1,]
dim(mafvcf.df_below0.1)
dim(mafvcf.df)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_766INDIV_3496SNPS_below0.1_50bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(mafvcf.df_below0.1$Frequency,breaks=seq(0,0.1,l=50), main="MAF frequency",
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()


##### below 0.05

mafvcf.df_below0.05<-mafvcf.df[mafvcf.df$Frequency<=0.05,]
dim(mafvcf.df_below0.05)
dim(mafvcf.df)

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_766INDIV_3496SNPS_below0.05_50bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(mafvcf.df_below0.05$Frequency,breaks=seq(0,0.05,l=50), main="MAF frequency",
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()




###############################################
############# INPUT DATA FAST CHEAT ################
###################################################

#loading the data file
#filtered.VCF<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS.vcf")





#open shortened file
pop.data <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/popmap.766_samples_header_region_rmINdiv_sorted.txt", sep = "\t", header = TRUE)

#pop.data$onecluster<-rep("one")

head(pop.data)

#all(colnames(filtered.VCF@gt)[-1] == pop.data$seID)

### genlight from vcf #####
gl.toad <- vcfR2genlight(filtered.VCF)

## add region pop as pop to genlight object
pop(gl.toad) <- pop.data$fourclusters

#setting ploidy
ploidy(gl.toad) <-2
pop(gl.toad)


gl.toad

" /// GENLIGHT OBJECT /////////

 // 766 genotypes,  3,496 binary SNPs, size: 3.4 Mb
 349205 (13.04 %) missing data

 // Basic content
   @gen: list of 766 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  766 individual labels
   @loc.names:  3496 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 195-287)
   @other: a list containing: elements without names "


#######################
## HWE filtering per pond p=0.01 ####
######################
#https://grunwaldlab.github.io/Population_Genetics_in_R/Locus_Stats.html

gl.toad.nosibs<-gl.toad

## check new pop.data file
dim(pop.data)

#### genind dartR from genligght ####
genid.toad.nosibs<-dartR::gl2gi(gl.toad.nosibs, probar = FALSE, verbose = NULL)


library("pegas")
#(nanhwe.full <- hw.test(genid.toad.nosibs, B = 1000)) # performs 1000 permuatations

## hwe per pop

pop(genid.toad.nosibs)<-pop.data$pop

### check how many pops
pop(genid.toad.nosibs) # 44

#run hwe test for each pop seprately - and only focus on the analytical P value (B=0)
(nanhwe.pop <- seppop(genid.toad.nosibs) %>% lapply(hw.test, B = 0))

# Take the third column with all rows - p value column
(nanhwe.mat <- sapply(nanhwe.pop, "[", i = TRUE, j = 3)) 

#alpha  <- 0.05
#newmat <- nanhwe.mat
#newmat[newmat > alpha] <- 1


#par(mar = c(1, 1, 1, 1))
library("lattice")
#levelplot(t(newmat))

#summary(newmat)

#plot(newmat)

############## my code #####

#turn to dataframe
nanhwe.mat.df<-as.data.frame(nanhwe.mat)

str(nanhwe.mat.df) #1868

#View(nanhwe.mat.df)

### subset loci to remove - ie anything under 0.05 #####
library(dplyr)

nanhwe.mat.df3<-nanhwe.mat.df %>%  
  filter_all(any_vars(. < 0.01))

str(nanhwe.mat.df3) #
#View(nanhwe.mat.df3)
summary(nanhwe.mat.df3)


###### loci to keep! remember exclude NAs
nanhwe.mat.df6<-nanhwe.mat.df %>%  
  filter_all(all_vars(. >= 0.01|is.na(.)))

str(nanhwe.mat.df6)#

# double check no minimum values less than 0.05
summary(nanhwe.mat.df6)



###### MAKE SURE CHANGE NUMBER OF BREEDING SITES TO pop(genid.toad.nosibs) 
# -------------------- variations ########
nanhwe.mat.df$INhwe_count <- apply( nanhwe.mat.df[ ,1:44], 1, function(x) sum( x > 0.01 ))

#keep.loci1 <- nanhwe.mat.df[ which( nanhwe.mat.df$INhwe_count==3 ), ]
#nrow(keep.loci1) #988


#pick how many sites out of 26 have to be in hwe to keep snp
2/3*44 = 29.3

keep.loci2 <- nanhwe.mat.df[ which( nanhwe.mat.df$INhwe_count >=29 ), ]
nrow(keep.loci2) #  1526

#keep.loci3 <- nanhwe.mat.df[ which( nanhwe.mat.df$INhwe_count >=1 ), ]
#nrow(keep.loci3) #1853


#keep.loci4 <- nanhwe.mat.df[ which( nanhwe.mat.df$INhwe_count >=1 ), ]
#nrow(keep.loci4) #


####################################################
### select just the row names ###
hwe.keep.loci.names <- rownames(keep.loci2)
str(hwe.keep.loci.names)

#subset
gl.toad.nosibs.hwe <- gl.toad.nosibs[ , hwe.keep.loci.names]
gl.toad.nosibs.hwe

"
 /// GENLIGHT OBJECT /////////

 // 766 genotypes,  1,526 binary SNPs, size: 2.1 Mb
 145074 (12.41 %) missing data

 // Basic content
   @gen: list of 766 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  766 individual labels
   @loc.names:  1526 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 195-287)
   @other: a list containing: elements without names 

"

##############################
####### CALL RATE #######
#################################

genid.toad.nosibs.hwe<-dartR::gl2gi(gl.toad.nosibs.hwe, probar = FALSE, verbose = NULL)

"/// GENIND OBJECT /////////

 // 766 individuals; 1,526 loci; 3,052 alleles; size: 9.8 Mb

 // Basic content
   @tab:  766 x 3052 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 3052 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: df2genind(X = xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
    pop = gl@pop, NA.char = "-", ploidy = 2)

 // Optional content
   @pop: population of each individual (group size range: 195-287)"

genid.toad.nosibs.hwe.df<-genid.toad.nosibs.hwe@tab
str(genid.toad.nosibs.hwe.df)
genid.toad.nosibs.hwe.df<-as.data.frame(genid.toad.nosibs.hwe.df)
head(genid.toad.nosibs.hwe.df)

SNPscalled <- data.frame("Internal_ID"= row.names(genid.toad.nosibs.hwe.df), "Called"=rowSums(genid.toad.nosibs.hwe.df != "00"), "Call.Rate"=(rowSums(genid.toad.nosibs.hwe.df != "00")/ncol(genid.toad.nosibs.hwe.df)))
SNPscalled
hist(SNPscalled$Call.Rate)
hist(SNPscalled$Call.Rate,breaks=seq(0.60,1,l=100))
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/CallRatePerPond_plot_766INDIV_1526SNPS_hwe0.01_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(SNPscalled$Call.Rate,breaks=seq(0.60,1,l=100))
dev.off()


#################################
###        FIS filter            ###################
#########################


#### genind dartR from genligght ####
genidtoad.nosibs.hwe<-dartR::gl2gi(gl.toad.nosibs.hwe, probar = FALSE, verbose = NULL)
genidtoad.nosibs.hwe
""

pop(genidtoad.nosibs.hwe)<-pop.data$pop
pop(genidtoad.nosibs.hwe)

### check ###


#### heirfstat conversion from genind ########
hierftoad.nosibs.hwe<-genind2hierfstat(genidtoad.nosibs.hwe)

#hierftoad.nosibs.hwe<-genind2hierfstat(genidtoad.nosibs.hwe,pop=TRUE)
#hierftoad.nosibs.hwe$pop<-pop.data$pop

hierftoad.nosibs.hwebasic<-basic.stats(hierftoad.nosibs.hwe)

hierftoad.nosibs.hwebasic$Fis

head(hierftoad.nosibs.hwebasic$Fis)



### add new row names column
hierftoad.nosibs.hwebasic$SNP<-xxx

FISdataperpond<-hierftoad.nosibs.hwebasic$Fis

FISdataperpondavg <- data.frame("SNP"=rownames(FISdataperpond), "Average"=rowMeans(FISdataperpond, na.rm=T))
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))

#png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/FISperpond_plot_766INDIV_1526SNPS_hwe0.01_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))
dev.off()

#png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/FISperpond_plot_766INDIV_1526SNPS_hwe0.01_50bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=50))
dev.off()




######### filter

# snps I want to keep
nrow(FISdataperpondavg[which(FISdataperpondavg$Average<0.4 & FISdataperpondavg$Average>-0.35),]) 
#1503

# snps I want to delete
nrow(FISdataperpondavg[which(FISdataperpondavg$Average>0.4 & FISdataperpondavg$Average>(-0.35)),]) 
#149

str(FISdataperpondavg)

SNPkeepaveFIS<-FISdataperpondavg[which(FISdataperpondavg$Average<0.4 & FISdataperpondavg$Average>-0.35),]

SNPkeep<-SNPkeepaveFIS$SNP
str(SNPkeep) #chr [1:1503]

# edit names

SNPkeep2<- sub( '\\.', ':', SNPkeep)



gl.toad.nosibs.hwe.FIS <- gl.toad.nosibs[ , SNPkeep2]
gl.toad.nosibs.hwe.FIS

"/// GENLIGHT OBJECT /////////

 // 766 genotypes,  1,503 binary SNPs, size: 2.1 Mb
 143392 (12.45 %) missing data

 // Basic content
   @gen: list of 766 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  766 individual labels
   @loc.names:  1503 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 195-287)
   @other: a list containing: elements without names 
"

#######################
###### SFS SITE FREQUENCY SPECTRUM ##########
#########################


###########################################
#### SPLIT INTO THREE REGIONS ####
###########################################


pop(gl.toad.nosibs.hwe.FIS) <- pop.data$fourclusters

gl.toad.nosibs.hwe.FIS.bypop<-seppop(gl.toad.nosibs.hwe.FIS)


gl.toad.nosibs.hwe.FIS.VanIsland<-gl.toad.nosibs.hwe.FIS.bypop$VanIsland
gl.toad.nosibs.hwe.FIS.VanIsland

" /// GENLIGHT OBJECT /////////

 // 195 genotypes,  1,503 binary SNPs, size: 641.3 Kb
 37256 (12.71 %) missing data

 // Basic content
   @gen: list of 195 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  195 individual labels
   @loc.names:  1503 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 195-195)
   @other: a list containing: elements without names  "



gl.toad.nosibs.hwe.FIS.LowerMain<-gl.toad.nosibs.hwe.FIS.bypop$LowerMain
gl.toad.nosibs.hwe.FIS.LowerMain
" /// GENLIGHT OBJECT /////////

 // 284 genotypes,  1,503 binary SNPs, size: 863.5 Kb
 51431 (12.05 %) missing data

 // Basic content
   @gen: list of 284 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  284 individual labels
   @loc.names:  1503 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 284-284)
   @other: a list containing: elements without names 
"

gl.toad.nosibs.hwe.FIS.HaidaGwai<-gl.toad.nosibs.hwe.FIS.bypop$HaidaGwai
gl.toad.nosibs.hwe.FIS.HaidaGwai
"  /// GENLIGHT OBJECT /////////

 // 239 genotypes,  1,503 binary SNPs, size: 766.2 Kb
 48068 (13.38 %) missing data

 // Basic content
   @gen: list of 239 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  239 individual labels
   @loc.names:  1503 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 239-239)
   @other: a list containing: elements without names  "



gl.toad.nosibs.hwe.FIS.Northwest<-gl.toad.nosibs.hwe.FIS.bypop$Northwest
gl.toad.nosibs.hwe.FIS.Northwest

"/// GENLIGHT OBJECT /////////

 // 48 genotypes,  1,503 binary SNPs, size: 245.7 Kb
 6637 (9.2 %) missing data

 // Basic content
   @gen: list of 48 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  48 individual labels
   @loc.names:  1503 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 48-48)
   @other: a list containing: elements without names "

dartR::gl.report.maf(gl.toad.nosibs.hwe.FIS.HaidaGwai, maf.limit = 0.5, ind.limit = 0, loc.limit = 0,
              v = 2)

### convert genlight to genind objects 
genind.toad.nosibs.hwe.FIS.VanIsland<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS.VanIsland, probar = FALSE, verbose = NULL)


""


genind.toad.nosibs.hwe.FIS.LowerMain<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS.LowerMain, probar = FALSE, verbose = NULL)

""


genind.toad.nosibs.hwe.FIS.HaidaGwai<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS.HaidaGwai, probar = FALSE, verbose = NULL)

''

genind.toad.nosibs.hwe.FIS.Northwest<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS.Northwest, probar = FALSE, verbose = NULL)


###########################################
#### convert to structure ####
###########################################

########### genind2structureL ######
#https://github.com/lvclark/R_genetics_conv/blob/master/genind2structure.R

source("F:/GBS_data_03_02_21/genind2structureL_function.R")



genind2structureL(genind.toad.nosibs.hwe.FIS.HaidaGwai, file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/HG_genind_structure_419NDIV_1367SNPS_regionID.str", pops=TRUE)


genind2structureL(genind.toad.nosibs.hwe.FIS.VanIsland, file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/VI_VanIsland_genind_structure_419NDIV_1367SNPS_regionID.str", pops=TRUE)


genind2structureL(genind.toad.nosibs.hwe.FIS.LowerMain, file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/ML_LowerMain_genind_structure_419NDIV_1367SNPS_regionID.str", pops=TRUE)

###########################################
#### convert STR TO vcf USING PDF SPIDER ####
###########################################





###########################################
#### UPLOAD three VCFs ####
###########################################


###########################################
#### plot folded MAF - as count of alleles, for each region ####
###########################################


### HG #####
##############################
#loading the data file
HG_filtered.VCF<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/HaidaGwai_267indiv_766_2392.vcf")


### subset to just 1503 loci
subset(HG_filtered.VCF, HG_filtered.VCF@fix[,3]%in%locinamnes_766_1503)

### alternative subset
HG_filtered.VCF_766_1503<-HG_filtered.VCF[HG_filtered.VCF@fix[,3]%in%locinamnes_766_1503,]

#minor allele freq
HG_mafvcf<-maf(HG_filtered.VCF_766_1503,element=2)



HG_mafvcf.df<-as.data.frame(HG_mafvcf)
summary(HG_mafvcf.df)

head(HG_mafvcf.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(2,2,2,2))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_HG_419NDIV_1367SNPS_freq.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(HG_mafvcf.df$Frequency,breaks=seq(0,0.5,l=50), main=" ",
     ylim = c(0, 400),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_HG_419NDIV_1367SNPS_count.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(HG_mafvcf.df$Count, main=" ",
     ylim = c(0, 600),
     xlab="Minor allele count",
     ylab="Frequency")
dev.off()

HGcounthist<-hist(HG_mafvcf.df$Count, main=" ",
                  ylim = c(0, 2000),
                  xlab="Minor allele count",
                  ylab="Frequency")




### VI #####
##############################
#loading the data file
VI_filtered.VCF<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/VanIsland_33indiv_766_2392.vcf")



### subset to just 1503 loci
subset(VI_filtered.VCF, VI_filtered.VCF@fix[,3]%in%locinamnes_766_1503)


VI_filtered.VCF_766_1503<-VI_filtered.VCF[VI_filtered.VCF@fix[,3]%in%locinamnes_766_1503,]

#minor allele freq
VI_mafvcf<-maf(VI_filtered.VCF_766_1503,element=2)




VI_mafvcf.df<-as.data.frame(VI_mafvcf)

head(VI_mafvcf.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(2,2,2,2))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_VI_419NDIV_1367SNPS_freq.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(VI_mafvcf.df$Frequency,breaks=seq(0,0.5,l=50), main=" ",
     ylim = c(0, 400),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_VI_419NDIV_1367SNPS_count.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(VI_mafvcf.df$Count, main=" ",
     ylim = c(0, 600),
     xlab="Minor allele count",
     ylab="Frequency")
dev.off()


VIcounthist<-hist(VI_mafvcf.df$Count, main=" ",
                  ylim = c(0, 2000),
                  xlab="Minor allele count",
                  ylab="Frequency")

### LM #####
##############################
#loading the data file
LM_filtered.VCF<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/LowerMain_71indiv_766_2392.vcf")

### subset to just 1503 loci
subset(LM_filtered.VCF, LM_filtered.VCF@fix[,3]%in%locinamnes_766_1503)


LM_filtered.VCF_766_1503<-LM_filtered.VCF[LM_filtered.VCF@fix[,3]%in%locinamnes_766_1503,]

#minor allele freq
LM_mafvcf<-maf(LM_filtered.VCF_766_1503,element=2)




LM_mafvcf.df<-as.data.frame(LM_mafvcf)

head(LM_mafvcf.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(2,2,2,2))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_LM_419NDIV_1367SNPS_freq.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(LM_mafvcf.df$Frequency,breaks=seq(0,0.5,l=50), main=" ",
     ylim = c(0, 400),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_LM_419NDIV_1367SNPS_count.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(LM_mafvcf.df$Count, main=" ",
     ylim = c(0, 600),
     xlab="Minor allele count",
     ylab="Frequency")
dev.off()



LMcounthist<-hist(LM_mafvcf.df$Count, main=" ",
     ylim = c(0, 2000),
     xlab="Minor allele count",
     ylab="Number of SNPs")






######################################
############### make r facet plot of SFS ###########
####################################
HG_mafvcf.df$Region<-rep("HaidaGwai")
VI_mafvcf.df$Region<-rep("VanIsland")
LM_mafvcf.df$Region<-rep("LowerMain")

allthreeregions<-rbind(HG_mafvcf.df,VI_mafvcf.df)
allthreeregions<-rbind(allthreeregions,LM_mafvcf.df)

head(allthreeregions)
dim(allthreeregions) # 2196    5
## loci name vs freq

### maf - freq no bin
labels4clustersgenetics <- c(VanIsland = "Vancouver Island", LowerMain = "Lower mainland", HaidaGwai = "Haida Gwaii")

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/SFS_4clusters_MAF_vs_locicount_766_1503.png", width = 3.5, height = 6.5, units = 'in', res = 300)
ggplot(data=allthreeregions, aes(x=Frequency,fill=Region)) +
  geom_histogram()+
  scale_fill_manual(values=c("#ab4e03", "#fca45d","#009E73"),
                    name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  facet_wrap(~Region,ncol=1,labeller = labeller(Region = labels4clustersgenetics))+
  ylab("Number of loci") + xlab("Minor allele frequency")+
  theme_bw()+
  xlim(0,0.5)+
  theme(legend.position="none",panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()

ggsave("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/SFS_4clusters_MAF_vs_locicount_766_1503.svg", width = 3.5, height = 6.5, units = 'in')
ggplot(data=allthreeregions, aes(x=Frequency,fill=Region)) +
  geom_histogram()+
  scale_fill_manual(values=c("#ab4e03", "#fca45d","#009E73"),
                    name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  facet_wrap(~Region,ncol=1,labeller = labeller(Region = labels4clustersgenetics))+
  ylab("Number of loci") + xlab("Minor allele frequency")+
  theme_bw()+
  xlim(0,0.5)+
  theme(legend.position="none",panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()






######### weird other stuff ######


ggplot(data=allthreeregions, aes(x=lociname, y=Frequency,fill=Region)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#ab4e03", "#fca45d","#009E73"),
                    name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  facet_wrap(~Region,ncol=1)+
  theme_bw()+
  theme(legend.position="none",panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

### freq vs count
ggplot(data=allthreeregions, aes(x=Frequency, y=Count,fill=Region)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#ab4e03", "#fca45d","#009E73"),
                    name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  facet_wrap(~Region,ncol=1)+
  theme_bw()+
  theme(legend.position="none",panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())



#### maf - freq - bin 0.01
ggplot(data=allthreeregions, aes(x=Frequency,fill=Region)) +
  geom_histogram(binwidth = 0.01)+
  scale_fill_manual(values=c("#ab4e03", "#fca45d","#009E73"),
                    name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  facet_wrap(~Region,ncol=1)+
  theme_bw()+
  theme(legend.position="none",panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())



### subset to just minor allele
sum(LM_mafvcf.df$Count)
sum(HG_mafvcf.df$Count)
sum(VI_mafvcf.df$Count)

LM_mafvcf.df$count_div_sumcount<-LM_mafvcf.df$Count/sum(LM_mafvcf.df$Count)
VI_mafvcf.df$count_div_sumcount<-VI_mafvcf.df$Count/sum(VI_mafvcf.df$Count)
HG_mafvcf.df$count_div_sumcount<-HG_mafvcf.df$Count/sum(HG_mafvcf.df$Count)


ggplot(data=LM_mafvcf.df, aes(x=Frequency, y=Count)) +
  geom_bar(stat="identity")

ggplot(data=VI_mafvcf.df, aes(x=Frequency, y=Count)) +
  geom_bar(stat="identity")

ggplot(data=HG_mafvcf.df, aes(x=Frequency, y=Count)) +
  geom_bar(stat="identity")#+
#theme_bw()

ggplot(data=HG_mafvcf.df, aes(x=Frequency, y=Count)) +
  geom_bar(stat="identity")#+
#theme_bw()


HG_mafvcf.df$lociname<-row.names(HG_mafvcf.df)
ggplot(data=HG_mafvcf.df, aes(x=lociname, y=Frequency)) +
  geom_bar(stat="identity")


LM_mafvcf.df$lociname<-row.names(LM_mafvcf.df)
ggplot(data=LM_mafvcf.df, aes(x=lociname, y=Frequency)) +
  geom_bar(stat="identity")

VI_mafvcf.df$lociname<-row.names(VI_mafvcf.df)
ggplot(data=VI_mafvcf.df, aes(x=lociname, y=Frequency)) +
  geom_bar(stat="identity")


###############################################
############# INPUT DATA THIS IS THE ONE ################
###################################################

#loading the data file
filtered.VCF<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS.vcf")





#open shortened file
pop.data <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/popmap.766_samples_header_region_rmINdiv_sorted.txt", sep = "\t", header = TRUE)

#pop.data$onecluster<-rep("one")

head(pop.data)

#all(colnames(filtered.VCF@gt)[-1] == pop.data$seID)

### genlight from vcf #####
gl.toad <- vcfR2genlight(filtered.VCF)

## add region pop as pop to genlight object
pop(gl.toad) <- pop.data$fourclusters

#setting ploidy
ploidy(gl.toad) <-2
pop(gl.toad)


gl.toad

" 
/// GENLIGHT OBJECT /////////

 // 766 genotypes,  3,496 binary SNPs, size: 3.4 Mb
 349205 (13.04 %) missing data

 // Basic content
   @gen: list of 766 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  766 individual labels
   @loc.names:  3496 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 195-287)
   @other: a list containing: elements without names 


"




########################## BASIC STATS no filter ###################

#### genind dartR from genligght ####
genid.toad<-dartR::gl2gi(gl.toad, probar = FALSE, verbose = NULL)

#set pop as region
pop(genid.toad)<-pop.data$fourclusters


#has pop
genid.toad@pop

#### heirfstat conversion from genind ########
hierf.toad<-genind2hierfstat(genid.toad,pop=NULL)

#has pop - as region
hierf.toad$pop



dataIN<-hierf.toad
pop.stats <- function(dataIN = dataIN, boots = 10000){
  
  # dataIN = population genomic data in hierfstat format
  
  #- calculate stats and fis bootstrapping
  tmp.stats <- hierfstat::basic.stats( dataIN )
  tmp.boot <- hierfstat::boot.ppfis( dataIN, nboot = boots )
  
  #- create output df
  
  tmp.summary <- matrix( nrow = 8, ncol = length( levels(dataIN$pop)),
                         dimnames = list( c( "Hs_mean", "Hs_SD", "Ho_mean", "Ho_SD", "P", "Fis", "Fis_ll", "Fis_hl" ),
                                          levels(dataIN$pop) )
  )
  
  #- calculate mean and sd for Hs and Ho from individual loci estimates
  
  tmp.summary[ "Hs_mean", ] <- apply( tmp.stats$Hs, 2, function(x) mean( x, na.rm = T ) )
  tmp.summary[ "Hs_SD", ] <- apply( tmp.stats$Hs, 2, function(x) sd( x, na.rm = T) )
  
  tmp.summary[ "Ho_mean", ] <- apply( tmp.stats$Ho, 2, function(x) mean( x, na.rm = T ) )
  tmp.summary[ "Ho_SD", ] <- apply( tmp.stats$Ho, 2, function(x) sd( x, na.rm = T) )
  
  #- calculate percentage polymorphic loci from 'pop.freq' data
  
  tmp.freq <- do.call( rbind, lapply( tmp.stats$pop.freq, function(x) x[1,] ))
  tmp.P <- apply( tmp.freq, 2, function(x) 1 - ( length( x[ x == 0 | x == 1 ] ) / length(x) ) )
  tmp.summary[ "P", ] <- round( tmp.P*100, 1)
  
  #- calculate per population Fis as per hierfstat's formula: Fis = 1 - Ho/Hs.
  #- get 95% CI for Fist from bootstrap calculations
  
  tmp.summary[ "Fis", ] <- 1 - ( tmp.summary[ 'Ho_mean', ] / tmp.summary[ 'Hs_mean', ] )
  tmp.summary[ "Fis_ll", ] <- tmp.boot$fis.ci$ll
  tmp.summary[ "Fis_hl", ] <- tmp.boot$fis.ci$hl
  
  #- add overall stats
  
  tmp.overall <- c( tmp.stats$overall[[ 'Hs' ]],
                    sd( tmp.stats$Hs, na.rm = T ),
                    tmp.stats$overall[[ 'Ho' ]],
                    sd( tmp.stats$Ho, na.rm = T),
                    NA, # polymorphic loci
                    tmp.stats$overall[[ 'Fis' ]],
                    NA, NA) # overall Fis 95% CI
  
  tmp.summary <- cbind( tmp.summary, "OVERALL" = tmp.overall )
  
  return( as.data.frame( t( round( tmp.summary, 4 ) )))
  
}

forb2.pop.stats <- pop.stats(hierf.toad)
forb2.pop.stats

### filter max het!! no negative FIS

"    "

### pop stats output #####

write.table(forb2.pop.stats, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/hierfstat_popstats_region_766INDIV_3496SNPS.txt",sep = "\t",row.names = TRUE,col.names = TRUE)


########################################
############ Fis per region CI no filter ################
#########################################
str(forb2.pop.stats)

forb2.pop.stats$Fis

Fis_95CI_region<-cbind(forb2.pop.stats$Fis,forb2.pop.stats$Fis_ll,forb2.pop.stats$Fis_hl)
dfFis_95CI_region<-as.data.frame(Fis_95CI_region)

#delete last row
dfFis_95CI_region<-dfFis_95CI_region[-c(4), ] 


## add region column
dfFis_95CI_region$Region<-c("Vancouver Island","Lower Mainland", "Haida Gwaii")

#rename column
colnames(dfFis_95CI_region) <- c("Fis","Fis_ll","Fis_hl","Region")    # Applying colnames

## reorder
dfFis_95CI_region$Region<- factor(dfFis_95CI_region$Region, levels = c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
dfFis_95CI_region<- droplevels(dfFis_95CI_region)

## plot - no legend
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Inbreeding (Fis ± 95% CIs)")+theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p

### save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_95CI_per_region_766INDIV_3496SNPS_nofilter.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Inbreeding (Fis ± 95% CIs)")+theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p
dev.off()




#############################
######### PCA no filter ##########
###############################

##per region
pop(gl.toad) <- pop.data$fourclusters


my_pca <- glPca(gl.toad, nf=6, parallel = require("parallel"))


my_scores <- as.data.frame(my_pca$scores)
my_scores$pop <- pop(gl.toad)

#storing eigenvalues on a file


#Checking how many eigenvalues to keep
#plotting PCA
toad.pca.scores <- as.data.frame(my_pca$scores)
toad.pca.scores$pop <- pop(gl.toad)

#gl.toad.subset
"
"

library(ggplot2)

# get % pca for 1st two axes ????? it created 402 PCA axes, why????



my_pca$eig[c(1,2)]

#not a percent!!!!!!!!!
sum(my_pca$eig)

#### percentage variance explained for PC axes
100*my_pca$eig/sum(my_pca$eig)


100*my_pca$eig[c(1,2)]/sum(my_pca$eig)
#  22.849915  3.357097

##per region
pop(gl.toad) <- pop.data$fourclusters
toad.pca.scores$pop <- pop(gl.toad)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_4clusters_4axes_766INDIV_3496SNPS.png", width = 8, height = 6.5, units = 'in', res = 300)
set.seed(9)
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                            name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='Region') 
p<-p + ylab("PC2 (3.4% explained variance)") + xlab("PC1 (22.8% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="top",legend.text=element_text(size=16))
p<-p + guides(colour=guide_legend(nrow=2))
p
dev.off()

##per breeding site
pop(gl.toad) <- pop.data$pop
toad.pca.scores$pop <- pop(gl.toad)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_site_4axes_766INDIV_3496SNPS.png", width = 8, height = 6.5, units = 'in', res = 300)
set.seed(9)
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p<-p + labs(colour='pop') 
p<-p + ylab("PC2 (3.4% explained variance)") + xlab("PC1 (22.8% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="top",legend.text=element_text(size=16))
p<-p + guides(colour=guide_legend(nrow=12))
p
dev.off()

######################################
########## FST no filter ##############
########################################

#https://www.researchgate.net/post/Is_there_a_consensus_on_what_constitutes_a_high_and_low_value_of_FST_in_population_genetics

genid.toad<-dartR::gl2gi(gl.toad, probar = FALSE, verbose = NULL)


#### heirfstat conversion from genind ########
hierf.toad<-genind2hierfstat(genid.toad,pop=NULL)





#has pop - as moe regions
hierf.toad$pop

hierf.toad$pop <- pop.data$fourclusters
hierf.toad$pop <- factor(pop.data$fourclusters)

trial<-pairwise.WCfst(hierf.toad,diploid=T)
trial

"         HaidaGwai  LowerMain  Northwest  VanIsland
HaidaGwai        NA 0.38006037 0.54621109 0.39127661
LowerMain 0.3800604         NA 0.08791682 0.05543831
Northwest 0.5462111 0.08791682         NA 0.10370834
VanIsland 0.3912766 0.05543831 0.10370834         NA "


write.table(trial, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_fourclusters_766INDIV_3496SNPS.txt",sep = "\t",row.names = TRUE,col.names = TRUE )





##########################################
####### Rm sibs  ###########
#####################################

#list sibs to remove
listindiv_plink_sub_sib_rel_over_04_names_unique.char<-read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/listindiv_plink_sub_sib_rel_over_04_names_unique.char.txt")

listindiv_plink_sub_sib_rel_over_04_names_unique.char$V1

listsibs<-listindiv_plink_sub_sib_rel_over_04_names_unique.char$V1


# function alba made to remove siblings
source("F:/GBS_data_03_02_21/Function_only_remove_siblings_genlight_Alba_nov_8th_2021.R")

# run function
gl.toad.nosibs<- remove_sibs_genlight(genlight_object = gl.toad,
                                      names = listsibs)
gl.toad.nosibs


"/// GENLIGHT OBJECT /////////

 // 419 genotypes,  3,496 binary SNPs, size: 2.4 Mb
 236613 (16.15 %) missing data

 // Basic content
   @gen: list of 419 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  419 individual labels
   @loc.names:  3496 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 32-160)
   @other: a list containing: elements without names 
 
"

### edit pop file

pop.data<-pop.data[which(!pop.data$sample.id %in% listsibs),]
dim(pop.data) # 419   7

write.table(pop.data,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/popmap.419_samples_header_region_rmINdiv_rmsibs_sorted.txt",quote=FALSE,sep = "\t",row.names=FALSE)

# count pops and get list
str(pop.data$pop)

pop.data2<-pop.data

pop.data2$pop<-as.factor(pop.data2$pop)

levels(pop.data2$pop)

ponds_kept_419INDIV_3496SNPS<-levels(pop.data2$pop)




#######################
## HWE filtering per pond p=0.01 ####
######################
#https://grunwaldlab.github.io/Population_Genetics_in_R/Locus_Stats.html

#gl.toad.nosibs<-gl.toad

## check new pop.data file
dim(pop.data)

#### genind dartR from genligght ####
genid.toad.nosibs<-dartR::gl2gi(gl.toad.nosibs, probar = FALSE, verbose = NULL)


library("pegas")
#(nanhwe.full <- hw.test(genid.toad.nosibs, B = 1000)) # performs 1000 permuatations

## hwe per pop

pop(genid.toad.nosibs)<-pop.data$pop

### check how many pops
pop(genid.toad.nosibs) # 26

#run hwe test for each pop seprately - and only focus on the analytical P value (B=0)
(nanhwe.pop <- seppop(genid.toad.nosibs) %>% lapply(hw.test, B = 0))

# Take the third column with all rows - p value column
(nanhwe.mat <- sapply(nanhwe.pop, "[", i = TRUE, j = 3)) 

#alpha  <- 0.05
#newmat <- nanhwe.mat
#newmat[newmat > alpha] <- 1


#par(mar = c(1, 1, 1, 1))
library("lattice")
#levelplot(t(newmat))

#summary(newmat)

#plot(newmat)

############## my code #####

#turn to dataframe
nanhwe.mat.df<-as.data.frame(nanhwe.mat)

str(nanhwe.mat.df) #1868

#View(nanhwe.mat.df)

### subset loci to remove - ie anything under 0.05 #####
library(dplyr)

nanhwe.mat.df3<-nanhwe.mat.df %>%  
  filter_all(any_vars(. < 0.01))

str(nanhwe.mat.df3) #
#View(nanhwe.mat.df3)
summary(nanhwe.mat.df3)


###### loci to keep! remember exclude NAs
nanhwe.mat.df6<-nanhwe.mat.df %>%  
  filter_all(all_vars(. >= 0.01|is.na(.)))

str(nanhwe.mat.df6)#

# double check no minimum values less than 0.05
summary(nanhwe.mat.df6)


# -------------------- variations ########
nanhwe.mat.df$INhwe_count <- apply( nanhwe.mat.df[ ,1:44], 1, function(x) sum( x > 0.01 ))

#keep.loci1 <- nanhwe.mat.df[ which( nanhwe.mat.df$INhwe_count==3 ), ]
#nrow(keep.loci1) #988


#pick how many sites out of 26 have to be in hwe to keep snp
2/3*44 = 29.3

keep.loci2 <- nanhwe.mat.df[ which( nanhwe.mat.df$INhwe_count >=29 ), ]
nrow(keep.loci2) #1399


#keep.loci3 <- nanhwe.mat.df[ which( nanhwe.mat.df$INhwe_count >=1 ), ]
#nrow(keep.loci3) #1853


#keep.loci4 <- nanhwe.mat.df[ which( nanhwe.mat.df$INhwe_count >=1 ), ]
#nrow(keep.loci4) #349


####################################################
### select just the row names ###
hwe.keep.loci.names <- rownames(keep.loci2)

#subset
gl.toad.nosibs.hwe <- gl.toad.nosibs[ , hwe.keep.loci.names]
gl.toad.nosibs.hwe

"
 /// GENLIGHT OBJECT /////////

 // 419 genotypes,  1,399 binary SNPs, size: 1.2 Mb
 88746 (15.14 %) missing data

 // Basic content
   @gen: list of 419 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  419 individual labels
   @loc.names:  1399 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 32-160)
   @other: a list containing: elements without names 

"


#############################
######### PCA HWE ##########
###############################

##per region
pop(gl.toad.nosibs.hwe) <- pop.data$fourclusters


my_pca <- glPca(gl.toad.nosibs.hwe, nf=6, parallel = require("parallel"))


my_scores <- as.data.frame(my_pca$scores)
my_scores$pop <- pop(gl.toad.nosibs.hwe)

#storing eigenvalues on a file


#Checking how many eigenvalues to keep
#plotting PCA
toad.pca.scores <- as.data.frame(my_pca$scores)
toad.pca.scores$pop <- pop(gl.toad.nosibs.hwe)

#gl.toad.nosibs.hwe.subset
"
"

library(ggplot2)

# get % pca for 1st two axes ????? it created 402 PCA axes, why????



my_pca$eig[c(1,2)]

#not a percent!!!!!!!!!
sum(my_pca$eig)

#### percentage variance explained for PC axes
100*my_pca$eig/sum(my_pca$eig)


100*my_pca$eig[c(1,2)]/sum(my_pca$eig)
#  23.920595  1.906154

##per region
pop(gl.toad.nosibs.hwe) <- pop.data$fourclusters
toad.pca.scores$pop <- pop(gl.toad.nosibs.hwe)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_MOEpop_4axes_766INDIV_1211SNPS.png", width = 8, height = 6.5, units = 'in', res = 300)
set.seed(9)
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73"),
                            name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='moe region') 
p<-p + ylab("PC2 (1.9% explained variance)") + xlab("PC1 (23.9% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="top",legend.text=element_text(size=16))
p
dev.off()


######################################
########## FST HWE ##############
########################################

#https://www.researchgate.net/post/Is_there_a_consensus_on_what_constitutes_a_high_and_low_value_of_FST_in_population_genetics

genid.toad.nosibs.hwe<-dartR::gl2gi(gl.toad.nosibs.hwe, probar = FALSE, verbose = NULL)


#### heirfstat conversion from genind ########
hierf.toad.nosibs.hwe<-genind2hierfstat(genid.toad.nosibs.hwe,pop=NULL)





#has pop - as moe regions
hierf.toad.nosibs.hwe$pop

hierf.toad.nosibs.hwe$pop <- pop.data$fourclusters
hierf.toad.nosibs.hwe$pop <- factor(pop.data$fourclusters)

trial<-pairwise.WCfst(hierf.toad.nosibs.hwe,diploid=T)
trial

"      "

write.table(trial, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_4clusters_766INDIV_1211SNPS_hwe44pops0.01.txt",sep = "\t",row.names = TRUE,col.names = TRUE )



### FST per loc overall pops ########

basic.nosibs.hwe<-basic.stats(hierf.toad.nosibs.hwe)

perloc.basic.nosibs.hwe<-basic.nosibs.hwe$perloc

par(mar = c(3,3,3,3))

#get max and min for axes
summary(perloc.basic.nosibs.hwe$Fst)

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Hist_Fst_perlocus_metapop_601INDIV_577SNPS.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(perloc.basic.nosibs.hwe$Fst,breaks=seq(-0.00310,0.25260,l=50),
     main="Fst per locus - metapopulation",
     xlab="Fst",
     ylab="Frequency")
dev.off()



########################## BASIC STATS HWE ###################

#### genind dartR from genligght ####
genid.toad.nosibs.hwe<-dartR::gl2gi(gl.toad.nosibs.hwe, probar = FALSE, verbose = NULL)

#set pop as region
pop(genid.toad.nosibs.hwe)<-pop.data$fourclusters


#has pop
genid.toad.nosibs.hwe@pop

#### heirfstat conversion from genind ########
hierf.toad.nosibs.hwe<-genind2hierfstat(genid.toad.nosibs.hwe,pop=NULL)

#has pop - as region
hierf.toad.nosibs.hwe$pop



dataIN<-hierf.toad.nosibs.hwe
pop.stats <- function(dataIN = dataIN, boots = 10000){
  
  # dataIN = population genomic data in hierfstat format
  
  #- calculate stats and fis bootstrapping
  tmp.stats <- hierfstat::basic.stats( dataIN )
  tmp.boot <- hierfstat::boot.ppfis( dataIN, nboot = boots )
  
  #- create output df
  
  tmp.summary <- matrix( nrow = 8, ncol = length( levels(dataIN$pop)),
                         dimnames = list( c( "Hs_mean", "Hs_SD", "Ho_mean", "Ho_SD", "P", "Fis", "Fis_ll", "Fis_hl" ),
                                          levels(dataIN$pop) )
  )
  
  #- calculate mean and sd for Hs and Ho from individual loci estimates
  
  tmp.summary[ "Hs_mean", ] <- apply( tmp.stats$Hs, 2, function(x) mean( x, na.rm = T ) )
  tmp.summary[ "Hs_SD", ] <- apply( tmp.stats$Hs, 2, function(x) sd( x, na.rm = T) )
  
  tmp.summary[ "Ho_mean", ] <- apply( tmp.stats$Ho, 2, function(x) mean( x, na.rm = T ) )
  tmp.summary[ "Ho_SD", ] <- apply( tmp.stats$Ho, 2, function(x) sd( x, na.rm = T) )
  
  #- calculate percentage polymorphic loci from 'pop.freq' data
  
  tmp.freq <- do.call( rbind, lapply( tmp.stats$pop.freq, function(x) x[1,] ))
  tmp.P <- apply( tmp.freq, 2, function(x) 1 - ( length( x[ x == 0 | x == 1 ] ) / length(x) ) )
  tmp.summary[ "P", ] <- round( tmp.P*100, 1)
  
  #- calculate per population Fis as per hierfstat's formula: Fis = 1 - Ho/Hs.
  #- get 95% CI for Fist from bootstrap calculations
  
  tmp.summary[ "Fis", ] <- 1 - ( tmp.summary[ 'Ho_mean', ] / tmp.summary[ 'Hs_mean', ] )
  tmp.summary[ "Fis_ll", ] <- tmp.boot$fis.ci$ll
  tmp.summary[ "Fis_hl", ] <- tmp.boot$fis.ci$hl
  
  #- add overall stats
  
  tmp.overall <- c( tmp.stats$overall[[ 'Hs' ]],
                    sd( tmp.stats$Hs, na.rm = T ),
                    tmp.stats$overall[[ 'Ho' ]],
                    sd( tmp.stats$Ho, na.rm = T),
                    NA, # polymorphic loci
                    tmp.stats$overall[[ 'Fis' ]],
                    NA, NA) # overall Fis 95% CI
  
  tmp.summary <- cbind( tmp.summary, "OVERALL" = tmp.overall )
  
  return( as.data.frame( t( round( tmp.summary, 4 ) )))
  
}

forb2.pop.stats.hwe <- pop.stats(hierf.toad.nosibs.hwe)
forb2.pop.stats.hwe

### filter max het!! no negative FIS

"         Hs_mean  Hs_SD Ho_mean  Ho_SD    P     Fis  Fis_ll  Fis_hl
VanIsland      0.0735 0.0827  0.0769 0.1036 85.6 -0.0466 -0.0918  0.0000
LowerMain      0.0878 0.0907  0.0903 0.1133 88.2 -0.0281 -0.0742  0.0194
HaidaGwai      0.0353 0.0778  0.0394 0.1066 44.0 -0.1147 -0.1817 -0.0405
OVERALL  0.0655 0.0867  0.0688 0.1100   NA -0.0513      NA      NA"

### pop stats output #####

write.table(forb2.pop.stats.hwe, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/hierfstat_popstats_region_601INDIV_577SNPS_hwe44pops.txt",sep = "\t",row.names = TRUE,col.names = TRUE)


########################################
############ Fis per region CI ################
#########################################
str(forb2.pop.stats.hwe)

forb2.pop.stats.hwe$Fis

Fis_95CI_region<-cbind(forb2.pop.stats.hwe$Fis,forb2.pop.stats.hwe$Fis_ll,forb2.pop.stats.hwe$Fis_hl)
dfFis_95CI_region<-as.data.frame(Fis_95CI_region)

#delete last row
dfFis_95CI_region<-dfFis_95CI_region[-c(4), ] 


## add region column
dfFis_95CI_region$Region<-c("Vancouver Island","Lower Mainland", "Haida Gwaii")

#rename column
colnames(dfFis_95CI_region) <- c("Fis","Fis_ll","Fis_hl","Region")    # Applying colnames

## reorder
dfFis_95CI_region$Region<- factor(dfFis_95CI_region$Region, levels = c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
dfFis_95CI_region<- droplevels(dfFis_95CI_region)

## plot - no legend
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Inbreeding (Fis ± 95% CIs)")+theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p

### save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_95CI_per_region_601INDIV_577SNPS_hwe44pops.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Inbreeding (Fis ± 95% CIs)")+theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p
dev.off()



#################################
###        FIS filter            ###################
#########################


#### genind dartR from genligght ####
genidtoad.nosibs.hwe<-dartR::gl2gi(gl.toad.nosibs.hwe, probar = FALSE, verbose = NULL)
genidtoad.nosibs.hwe
""

pop(genidtoad.nosibs.hwe)<-pop.data$pop
pop(genidtoad.nosibs.hwe)

### check ###


#### heirfstat conversion from genind ########
hierftoad.nosibs.hwe<-genind2hierfstat(genidtoad.nosibs.hwe)

#hierftoad.nosibs.hwe<-genind2hierfstat(genidtoad.nosibs.hwe,pop=TRUE)
#hierftoad.nosibs.hwe$pop<-pop.data$pop

hierftoad.nosibs.hwebasic<-basic.stats(hierftoad.nosibs.hwe)

hierftoad.nosibs.hwebasic$Fis

head(hierftoad.nosibs.hwebasic$Fis)

FISdataperpond<-hierftoad.nosibs.hwebasic$Fis

FISdataperpondavg <- data.frame("SNP"=rownames(FISdataperpond), "Average"=rowMeans(FISdataperpond, na.rm=T))
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))

#png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/FISperpond_plot_419INDIV_1399SNPS_hwe0.01_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-2,2,l=100))
dev.off()

#png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/FISperpond_plot_7419INDIV_1399SNPS_hwe0.01_50bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-2,2,l=50))
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/FISperpond_plot_419INDIV_1399SNPS_hwe0.01_50bins_formatted_BEFOREfilter.png", width = 8, height = 6.5, units = 'in', res = 300)
par(mar = c(5, 5,2, 2)+ 0.3, xaxs = "i", yaxs = "i")
hist(FISdataperpondavg$Average,breaks=seq(-2,2,l=50),bty = "l",
     las = 1,cex.lab=1.6,main="",ylab="Frequency",xlab=expression(paste(italic(F)[plain(IS)],paste(" per locus"))))
dev.off()




######### filter

# snps I want to keep
nrow(FISdataperpondavg[which(FISdataperpondavg$Average<0.4 & FISdataperpondavg$Average>-0.35),]) 
# 1367

# snps I want to delete
nrow(FISdataperpondavg[which(FISdataperpondavg$Average>0.4 & FISdataperpondavg$Average>(-0.35)),]) 
#

str(FISdataperpondavg)

SNPkeepaveFIS<-FISdataperpondavg[which(FISdataperpondavg$Average<0.4 & FISdataperpondavg$Average>-0.35),]

SNPkeep<-SNPkeepaveFIS$SNP
str(SNPkeep) #chr [1:1367]

SNPkeep2<- sub( '\\.', ':', SNPkeep)



gl.toad.nosibs.hwe.FIS <- gl.toad.nosibs.hwe[ , SNPkeep2]
gl.toad.nosibs.hwe.FIS


" /// GENLIGHT OBJECT /////////

 // 419 genotypes,  1,367 binary SNPs, size: 1.2 Mb
 87224 (15.23 %) missing data

 // Basic content
   @gen: list of 419 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  419 individual labels
   @loc.names:  1367 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 32-160)
   @other: a list containing: elements without names  "



#### genind dartR from genligght ####
genidtoad.nosibs.hwe.FIS<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS, probar = FALSE, verbose = NULL)
genidtoad.nosibs.hwe.FIS
""

pop(genidtoad.nosibs.hwe.FIS)<-pop.data$pop
pop(genidtoad.nosibs.hwe.FIS)

### check ###


#### heirfstat conversion from genind ########
hierftoad.nosibs.hwe.FIS<-genind2hierfstat(genidtoad.nosibs.hwe.FIS)

#hierftoad.nosibs.hwe.FIS<-genind2hierfstat(genidtoad.nosibs.hwe.FIS,pop=TRUE)
#hierftoad.nosibs.hwe.FIS$pop<-pop.data$pop

hierftoad.nosibs.hwe.FISbasic<-basic.stats(hierftoad.nosibs.hwe.FIS)

hierftoad.nosibs.hwe.FISbasic$Fis

head(hierftoad.nosibs.hwe.FISbasic$Fis)

FISdataperpond<-hierftoad.nosibs.hwe.FISbasic$Fis

FISdataperpondavg <- data.frame("SNP"=rownames(FISdataperpond), "Average"=rowMeans(FISdataperpond, na.rm=T))
hist(FISdataperpondavg$Average,breaks=seq(-0.5,0.5,l=100))

#png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/FISperpond_plot_419NDIV_1367SNPS_hwe0.01_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-0.5,0.5,l=100))
dev.off()

#png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/FISperpond_plot_419NDIV_1367SNPS_hwe0.01_50bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-0.5,0.5,l=50))
dev.off()


## formatted after filter
#png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/FISperpond_plot_419NDIV_1367SNPS_hwe0.01_50bins_formatted_AFTEREfilter.png", width = 8, height = 6.5, units = 'in', res = 300)
par(mar = c(5, 5,2, 2)+ 0.3, xaxs = "i", yaxs = "i")
hist(FISdataperpondavg$Average,breaks=seq(-0.5,0.5,l=100)
,bty = "l",
     las = 1,cex.lab=1.6,main="",ylab="Frequency",xlab=expression(paste(italic(F)[plain(IS)],paste(" per locus"))))
dev.off()


###########################################
#### convert to structure ####
###########################################

########### genind2structureL ######
#https://github.com/lvclark/R_genetics_conv/blob/master/genind2structure.R

source("F:/GBS_data_03_02_21/genind2structureL_function.R")


#### genind dartR from genligght ####
genid.toad.nosibs.hwe.FIS<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS, probar = FALSE, verbose = NULL)


pop(genid.toad.nosibs.hwe.FIS)<-pop.data$fourclusters

#has pop
genid.toad.nosibs.hwe.FIS@pop


genind2structureL(genid.toad.nosibs.hwe.FIS, file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/genind_structure_419INDIV_1367SNPS_regionID.str", pops=TRUE)

#############################
######### PCA FIS ##########
###############################

##per region
pop(gl.toad.nosibs.hwe.FIS) <- pop.data$fourclusters


my_pca <- glPca(gl.toad.nosibs.hwe.FIS, nf=6, parallel = require("parallel"))


my_scores <- as.data.frame(my_pca$scores)
my_scores$pop <- pop(gl.toad.nosibs.hwe.FIS)

#storing eigenvalues on a file


#Checking how many eigenvalues to keep
#plotting PCA
toad.pca.scores <- as.data.frame(my_pca$scores)
toad.pca.scores$pop <- pop(gl.toad.nosibs.hwe.FIS)

#gl.toad.nosibs.hwe.FIS.subset
"
"

library(ggplot2)

# get % pca for 1st two axes ????? it created 402 PCA axes, why????



my_pca$eig[c(1,2)]

#not a percent!!!!!!!!!
sum(my_pca$eig)

#### percentage variance explained for PC axes
100*my_pca$eig/sum(my_pca$eig)


100*my_pca$eig[c(1,2)]/sum(my_pca$eig)
#  31.14959  2.71185

##per region
pop(gl.toad.nosibs.hwe.FIS) <- pop.data$fourclusters
toad.pca.scores$pop <- pop(gl.toad.nosibs.hwe.FIS)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_4clusters_4axes_419NDIV_1367SNPS_hweperpond_FISperpond.png", width = 8, height = 6.5, units = 'in', res = 300)
set.seed(9)
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                            name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='Region') 
p<-p + ylab("PC2 (2.71% explained variance)") + xlab("PC1 (31.1% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="top",legend.text=element_text(size=16))
p<-p + guides(colour=guide_legend(nrow=2))
p
dev.off()



##per pop
pop(gl.toad.nosibs.hwe.FIS) <- pop.data$pop
toad.pca.scores$pop <- pop(gl.toad.nosibs.hwe.FIS)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_perpop_4axes_419NDIV_1367SNPS_hweperpond_FISperpond.png", width = 8, height = 6.5, units = 'in', res = 300)
set.seed(9)
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
#p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
#                            name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='pop') 
p<-p + ylab("PC2 (2.71% explained variance)") + xlab("PC1 (31.1% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + guides(colour=guide_legend(nrow=5))
p
dev.off()



## NEW COLS 
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_MOEpop_4axes_419NDIV_1367SNPS_hweperpond_FISperpond_orangegrad_legright2.pdf",  bg = "transparent",width =300, height = 150, units = c("mm"))
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=0.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
#p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#ab4e03", "#fca45d","#009E73"),
                            name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='moe region') 
p<-p + ylab("PC2 (2.71% explained variance)") + xlab("PC1 (31.1% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="right",legend.text=element_text(size=18),legend.key.height=unit(2,"line"),panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank())
p
dev.off()





##############################
####### PC3 and 4 #############
##per region
pop(gl.toad.nosibs.hwe.FIS) <- pop.data$fourclusters


my_pca <- glPca(gl.toad.nosibs.hwe.FIS, nf=6, parallel = require("parallel"))


my_scores <- as.data.frame(my_pca$scores)
my_scores$pop <- pop(gl.toad.nosibs.hwe.FIS)

#storing eigenvalues on a file


#Checking how many eigenvalues to keep
#plotting PCA
toad.pca.scores <- as.data.frame(my_pca$scores)
toad.pca.scores$pop <- pop(gl.toad.nosibs.hwe.FIS)

#gl.toad.nosibs.hwe.FIS.subset
"
"

library(ggplot2)

# get % pca for 1st two axes ????? it created 402 PCA axes, why????



my_pca$eig[c(3,4)]

#not a percent!!!!!!!!!
sum(my_pca$eig)

#### percentage variance explained for PC axes
100*my_pca$eig/sum(my_pca$eig)


100*my_pca$eig[c(3,4)]/sum(my_pca$eig)
# 1.323785 1.263350

##per region
pop(gl.toad.nosibs.hwe.FIS) <- pop.data$fourclusters
toad.pca.scores$pop <- pop(gl.toad.nosibs.hwe.FIS)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_MOEpop_4axes_419NDIV_1367SNPS_hweperpond_FISperpond_PC3_PC4.png", width = 8, height = 6.5, units = 'in', res = 300)
set.seed(9)
p <- ggplot(toad.pca.scores, aes(x=PC3, y=PC4, colour=pop))
p <- p + geom_point(size=4, alpha=.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73"),
                            name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='moe region') 
p<-p + ylab("PC4 (1.3% explained variance)") + xlab("PC3 (1.3% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="top",legend.text=element_text(size=16))
p
dev.off()

## NEW COL
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_MOEpop_419NDIV_1367SNPS_hweperpond_FISperpond_orangegrad_legright2_PC3_PC4.pdf",  bg = "transparent",width =300, height = 150, units = c("mm"))
p <- ggplot(toad.pca.scores, aes(x=PC3, y=PC4, colour=pop))
p <- p + geom_point(size=4, alpha=0.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
#p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#ab4e03", "#fca45d","#009E73"),
                            name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='moe region') 
p<-p + ylab("PC4 (1.3% explained variance)") + xlab("PC3 (1.3% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="right",legend.text=element_text(size=18),legend.key.height=unit(2,"line"),panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank())
p
dev.off()


######################################
########## FST FIS ##############
########################################

genid.toad.nosibs.hwe.FIS<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS, probar = FALSE, verbose = NULL)
genid.toad.nosibs.hwe.FIS

"/// GENIND OBJECT /////////

 // 419 individuals; 1,367 loci; 2,734 alleles; size: 5.2 Mb

 // Basic content
   @tab:  419 x 2734 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 2734 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: df2genind(X = xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
    pop = gl@pop, NA.char = "-", ploidy = 2)

 // Optional content
   @pop: population of each individual (group size range: 1-22)"

#### heirfstat conversion from genind ########
hierf.toad.nosibs.hwe.FIS<-genind2hierfstat(genid.toad.nosibs.hwe.FIS,pop=NULL)





#has pop - as moe regions
hierf.toad.nosibs.hwe.FIS$pop

hierf.toad.nosibs.hwe.FIS$pop <- pop.data$fourclusters
hierf.toad.nosibs.hwe.FIS$pop <- factor(pop.data$fourclusters)

trial<-pairwise.WCfst(hierf.toad.nosibs.hwe.FIS,diploid=T)
trial

"      HaidaGwai  LowerMain Northwest  VanIsland
HaidaGwai        NA 0.50939527 0.6277858 0.55075195
LowerMain 0.5093953         NA 0.1045395 0.05055287
Northwest 0.6277858 0.10453954        NA 0.11580875
VanIsland 0.5507519 0.05055287 0.1158087         NA         "

write.table(trial, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_4clusters_419NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1.txt",sep = "\t",row.names = TRUE,col.names = TRUE )






#has pop - as three areas
hierf.toad.nosibs.hwe.FIS$pop

hierf.toad.nosibs.hwe.FIS$pop <- pop.data$threeclusters
hierf.toad.nosibs.hwe.FIS$pop <- factor(pop.data$threeclusters)

trial2<-pairwise.WCfst(hierf.toad.nosibs.hwe.FIS,diploid=T)
trial2

"             HaidaGwai  Northwest       swBC
HaidaGwai        NA 0.62778584 0.46341032
Northwest 0.6277858         NA 0.09491858
swBC      0.4634103 0.09491858         NA        "

write.table(trial2, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_3clusters_419NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1.txt",sep = "\t",row.names = TRUE,col.names = TRUE )


# as pop
hierf.toad.nosibs.hwe.FIS$pop

hierf.toad.nosibs.hwe.FIS$pop <- pop.data$pop
hierf.toad.nosibs.hwe.FIS$pop <- factor(pop.data$pop)

trial3<-pairwise.WCfst(hierf.toad.nosibs.hwe.FIS,diploid=T)
trial3

"               "

write.csv(trial3, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_POP_419NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1_take2.csv")




##### p value pop
pop(gl.toad.nosibs.hwe.FIS) <- pop.data$pop

(Fstp<-stamppFst(gl.toad.nosibs.hwe.FIS, nboots = 100, percent = 95, nclusters = 1))


write.csv(Fstp$Pvalues, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_pvalue_POP_419NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1.csv" )



##### p value two areas
pop(gl.toad.nosibs.hwe.FIS) <- pop.data$twoclusters

(Fstp<-stamppFst(gl.toad.nosibs.hwe.FIS, nboots = 100, percent = 95, nclusters = 1))


write.csv(Fstp$Pvalues, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_pvalue_2clusters_419NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1.csv" )

##### p value four areas
pop(gl.toad.nosibs.hwe.FIS) <- pop.data$fourclusters

(Fstp<-stamppFst(gl.toad.nosibs.hwe.FIS, nboots = 100, percent = 95, nclusters = 1))


write.csv(Fstp$Pvalues, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_pvalue_3regions_419NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1.csv" )





### FST per loc overall pops ########

hierf.toad.nosibs.hwe.FIS$pop <- pop.data$twoclusters
hierf.toad.nosibs.hwe.FIS$pop <- factor(pop.data$twoclusters)

basic.nosibs.hwe.FIS<-basic.stats(hierf.toad.nosibs.hwe.FIS)

perloc.basic.nosibs.hwe.FIS<-basic.nosibs.hwe.FIS$perloc

par(mar = c(3,3,3,3))

summary(perloc.basic.nosibs.hwe.FIS$Fst)

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Hist_Fst_perlocus_2regions_419NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1_2.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(perloc.basic.nosibs.hwe.FIS$Fst,breaks=seq(-0.003100,1,l=50),
     main = " ",
     xlab="Pairwise Fst per locus between Haida Gwaii and southwest BC",
     ylab="Frequency")
dev.off()


ggsave("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Hist_Fst_perlocus_2regions_419NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1.svg", width = 8, height = 6.5, units = 'in')
hist(perloc.basic.nosibs.hwe.FIS$Fst,breaks=seq(-0.003100,1,l=50),
     xlab="Pairwise Fst between Haida Gwaii and southwest BC",
     ylab="Frequency")
dev.off()

######## wc

(FSTperloc<-hierfstat::wc(hierf.toad.nosibs.hwe.FIS,diploid=TRUE))
perlocobject<-FSTperloc$per.loc

perlocobject$FST

summary(perlocobject$FST)

hist(perlocobject$FST,breaks=seq(-0.006251 ,1,l=50))

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Hist_Fst_perlocus_2regions_419NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1_wcfunc.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(perlocobject$FST,breaks=seq(-0.006251 ,1,l=50),
main = " ",
     xlab="Pairwise Fst per locus between Haida Gwaii and southwest BC",
     ylab="Frequency")
dev.off()

    
     
     ####### FST per loc per region #####

#### genind dartR from genligght ####
genidtoad.nosibs.hwe<-dartR::gl2gi(gl.toad.nosibs.hwe, probar = FALSE, verbose = NULL)
genidtoad.nosibs.hwe
""

pop(genidtoad.nosibs.hwe)<-pop.data$pop
pop(genidtoad.nosibs.hwe)

### check ###


#### heirfstat conversion from genind ########
hierftoad.nosibs.hwe<-genind2hierfstat(genidtoad.nosibs.hwe)









#################
######### fst per locus ######
#########################

#library(devtools)
#install_github('j-a-thia/genomalicious')

## run on vcftools

fstperloc_swbc_hg_vcftools<-read.csv("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/FST_VCFTOOLS_766_1503_swbc_HG_perloc.csv")

hist(fstperloc_swbc_hg_vcftools$WEIR_AND_COCKERHAM_FST)

summary(fstperloc_swbc_hg_vcftools$WEIR_AND_COCKERHAM_FST)


################################
##### FST within region  -  FST overall per pop FIS ########
#################################

########### three areas #############


pop(gl.toad.nosibs.hwe.FIS) <- pop.data$fourclusters

gl.toad.nosibs.hwe.FIS.bypop<-seppop(gl.toad.nosibs.hwe.FIS)


gl.toad.nosibs.hwe.FIS.VanIsland<-gl.toad.nosibs.hwe.FIS.bypop$VanIsland
gl.toad.nosibs.hwe.FIS.LowerMain<-gl.toad.nosibs.hwe.FIS.bypop$LowerMain
gl.toad.nosibs.hwe.FIS.HaidaGwai<-gl.toad.nosibs.hwe.FIS.bypop$HaidaGwai



### convert genlight to genind objects 
gind.toad.nosibs.hwe.FIS.VanIsland<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS.VanIsland, probar = FALSE, verbose = NULL)


""


gind.toad.nosibs.hwe.FIS.LowerMain<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS.LowerMain, probar = FALSE, verbose = NULL)

""


gind.toad.nosibs.hwe.FIS.HaidaGwai<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS.HaidaGwai, probar = FALSE, verbose = NULL)

''


### make new pop files for each dataset
pop.data_VanIsland <- filter(pop.data, 4clusters == "VanIsland")
pop.data_LowerMain <- filter(pop.data, 4clusters == "LowerMain")
pop.data_HaidaGwai <- filter(pop.data, 4clusters == "HaidaGwai")




##### VanIsland #############

# make pop breeding pond
pop(gind.toad.nosibs.hwe.FIS.VanIsland)<-pop.data_VanIsland$pop

#### heirfstat conversion from genind ########
hierf.toad.nosibs.hwe.FIS.VanIsland<-genind2hierfstat(gind.toad.nosibs.hwe.FIS.VanIsland,pop=NULL)

#has pop - as breeding pond
hierf.toad.nosibs.hwe.FIS.VanIsland$pop


(VanIslandFST<-hierfstat::wc(hierf.toad.nosibs.hwe.FIS.VanIsland,diploid=TRUE))




##### LowerMain #############


# make pop breeding pond
pop(gind.toad.nosibs.hwe.FIS.LowerMain)<-pop.data_LowerMain$pop

#### heirfstat conversion from genind ########
hierf.toad.nosibs.hwe.FIS.LowerMain<-genind2hierfstat(gind.toad.nosibs.hwe.FIS.LowerMain,pop=NULL)

#has pop - as breeding pond
hierf.toad.nosibs.hwe.FIS.LowerMain$pop


(LowerMainFST<-hierfstat::wc(hierf.toad.nosibs.hwe.FIS.LowerMain,diploid=TRUE))



##### HaidaGwai #############


# make pop breeding pond
pop(gind.toad.nosibs.hwe.FIS.HaidaGwai)<-pop.data_HaidaGwai$pop

#### heirfstat conversion from genind ########
hierf.toad.nosibs.hwe.FIS.HaidaGwai<-genind2hierfstat(gind.toad.nosibs.hwe.FIS.HaidaGwai,pop=NULL)

#has pop - as breeding pond
hierf.toad.nosibs.hwe.FIS.HaidaGwai$pop


(HaidaGwaiFST<-hierfstat::wc(hierf.toad.nosibs.hwe.FIS.HaidaGwai,diploid=TRUE))


######## save values

VanIslandFST1<-VanIslandFST$FST
LowerMainFST1<-LowerMainFST$FST
HaidaGwaiFST1<-HaidaGwaiFST$FST


FSTperpop<-cbind(VanIslandFST1,LowerMainFST1,HaidaGwaiFST1)


write.table(FSTperpop[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/FST_within_4clusters_766_1503.txt")



###################################
########### two clusters #############
########################################

pop(gl.toad.nosibs.hwe.FIS) <- pop.data$twoclusters

gl.toad.nosibs.hwe.FIS.bypop<-seppop(gl.toad.nosibs.hwe.FIS)


gl.toad.nosibs.hwe.FIS.HG<-gl.toad.nosibs.hwe.FIS.bypop$HG
gl.toad.nosibs.hwe.FIS.swBC<-gl.toad.nosibs.hwe.FIS.bypop$swBC



### convert genlight to genind objects 
gind.toad.nosibs.hwe.FIS.HG<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS.HG, probar = FALSE, verbose = NULL)


""


gind.toad.nosibs.hwe.FIS.swBC<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS.swBC, probar = FALSE, verbose = NULL)

""





### make new pop files for each dataset
pop.data_HG <- filter(pop.data, twoclusters == "HG")
pop.data_swBC <- filter(pop.data, twoclusters == "swBC")




##### HG #############

# make pop breeding pond
pop(gind.toad.nosibs.hwe.FIS.HG)<-pop.data_HG$pop

#### heirfstat conversion from genind ########
hierf.toad.nosibs.hwe.FIS.HG<-genind2hierfstat(gind.toad.nosibs.hwe.FIS.HG,pop=NULL)

#has pop - as breeding pond
hierf.toad.nosibs.hwe.FIS.HG$pop


(HGFST<-hierfstat::wc(hierf.toad.nosibs.hwe.FIS.HG,diploid=TRUE))




##### swBC #############


# make pop breeding pond
pop(gind.toad.nosibs.hwe.FIS.swBC)<-pop.data_swBC$pop

#### heirfstat conversion from genind ########
hierf.toad.nosibs.hwe.FIS.swBC<-genind2hierfstat(gind.toad.nosibs.hwe.FIS.swBC,pop=NULL)

#has pop - as breeding pond
hierf.toad.nosibs.hwe.FIS.swBC$pop


(swBCFST<-hierfstat::wc(hierf.toad.nosibs.hwe.FIS.swBC,diploid=TRUE))


######## save values

HGFST1<-HGFST$FST
swBCFST1<-swBCFST$FST


FSTperpop<-cbind(HGFST1,swBCFST1)


write.table(FSTperpop[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/FST_within_twoclusters_766_1503.txt")















########################## BASIC STATS FIS ###################

#### genind dartR from genligght ####
genid.toad.nosibs.hwe.FIS<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS, probar = FALSE, verbose = NULL)

#set pop as region
pop(genid.toad.nosibs.hwe.FIS)<-pop.data$fourclusters


#has pop
genid.toad.nosibs.hwe.FIS@pop

#### heirfstat conversion from genind ########
hierf.toad.nosibs.hwe.FIS<-genind2hierfstat(genid.toad.nosibs.hwe.FIS,pop=NULL)

#has pop - as region
hierf.toad.nosibs.hwe.FIS$pop




dataIN<-hierf.toad.nosibs.hwe.FIS
pop.stats <- function(dataIN = dataIN, boots = 10000){
  
  # dataIN = population genomic data in hierfstat format
  
  #- calculate stats and fis bootstrapping
  tmp.stats <- hierfstat::basic.stats( dataIN )
  tmp.boot <- hierfstat::boot.ppfis( dataIN, nboot = boots )
  
  #- create output df
  
  tmp.summary <- matrix( nrow = 8, ncol = length( levels(dataIN$pop)),
                         dimnames = list( c( "Hs_mean", "Hs_SD", "Ho_mean", "Ho_SD", "P", "Fis", "Fis_ll", "Fis_hl" ),
                                          levels(dataIN$pop) )
  )
  
  #- calculate mean and sd for Hs and Ho from individual loci estimates
  
  tmp.summary[ "Hs_mean", ] <- apply( tmp.stats$Hs, 2, function(x) mean( x, na.rm = T ) )
  tmp.summary[ "Hs_SD", ] <- apply( tmp.stats$Hs, 2, function(x) sd( x, na.rm = T) )
  
  tmp.summary[ "Ho_mean", ] <- apply( tmp.stats$Ho, 2, function(x) mean( x, na.rm = T ) )
  tmp.summary[ "Ho_SD", ] <- apply( tmp.stats$Ho, 2, function(x) sd( x, na.rm = T) )
  
  #- calculate percentage polymorphic loci from 'pop.freq' data
  
  tmp.freq <- do.call( rbind, lapply( tmp.stats$pop.freq, function(x) x[1,] ))
  tmp.P <- apply( tmp.freq, 2, function(x) 1 - ( length( x[ x == 0 | x == 1 ] ) / length(x) ) )
  tmp.summary[ "P", ] <- round( tmp.P*100, 1)
  
  #- calculate per population Fis as per hierfstat's formula: Fis = 1 - Ho/Hs.
  #- get 95% CI for Fist from bootstrap calculations
  
  tmp.summary[ "Fis", ] <- 1 - ( tmp.summary[ 'Ho_mean', ] / tmp.summary[ 'Hs_mean', ] )
  tmp.summary[ "Fis_ll", ] <- tmp.boot$fis.ci$ll
  tmp.summary[ "Fis_hl", ] <- tmp.boot$fis.ci$hl
  
  #- add overall stats
  
  tmp.overall <- c( tmp.stats$overall[[ 'Hs' ]],
                    sd( tmp.stats$Hs, na.rm = T ),
                    tmp.stats$overall[[ 'Ho' ]],
                    sd( tmp.stats$Ho, na.rm = T),
                    NA, # polymorphic loci
                    tmp.stats$overall[[ 'Fis' ]],
                    NA, NA) # overall Fis 95% CI
  
  tmp.summary <- cbind( tmp.summary, "OVERALL" = tmp.overall )
  
  return( as.data.frame( t( round( tmp.summary, 4 ) )))
  
}

forb2.pop.stats.hwe.FIS <- pop.stats(hierf.toad.nosibs.hwe.FIS)
forb2.pop.stats.hwe.FIS

### filter max het!! no negative FIS

"          Hs_mean  Hs_SD Ho_mean  Ho_SD    P     Fis  Fis_ll  Fis_hl
VanIsland  0.1211 0.0954  0.1203 0.1002 92.3  0.0060 -0.0038  0.0160
LowerMain  0.1387 0.0919  0.1368 0.0958 97.9  0.0143  0.0065  0.0220
HaidaGwai  0.0132 0.0577  0.0156 0.0754 24.1 -0.1816 -0.2431 -0.1184
Northwest  0.1104 0.1484  0.1174 0.1645 55.2 -0.0637 -0.0819 -0.0457
OVERALL    0.0958 0.1144  0.0975 0.1236   NA -0.0179      NA      NA "

### pop stats output #####

write.table(forb2.pop.stats.hwe.FIS, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/hierfstat_popstats_region_419NDIV_1367SNPS_hwe44pops_FISperpond.txt",sep = "\t",row.names = TRUE,col.names = TRUE)


########################################
############ Fis per region CI ################
#########################################
str(forb2.pop.stats.hwe.FIS)

forb2.pop.stats.hwe.FIS$Fis

Fis_95CI_region<-cbind(forb2.pop.stats.hwe.FIS$Fis,forb2.pop.stats.hwe.FIS$Fis_ll,forb2.pop.stats.hwe$Fis_hl)
dfFis_95CI_region<-as.data.frame(Fis_95CI_region)

#delete last row
dfFis_95CI_region<-dfFis_95CI_region[-c(4), ] 


## add region column
dfFis_95CI_region$Region<-c("Vancouver Island","Lower Mainland", "Haida Gwaii")

#rename column
colnames(dfFis_95CI_region) <- c("Fis","Fis_ll","Fis_hl","Region")    # Applying colnames

## reorder
dfFis_95CI_region$Region<- factor(dfFis_95CI_region$Region, levels = c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
dfFis_95CI_region<- droplevels(dfFis_95CI_region)

## plot - no legend
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Inbreeding (Fis ± 95% CIs)")+theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p

### save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_95CI_per_region_419NDIV_1367SNPS_hwe44pops_FISperpond0.1.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Inbreeding (Fis ± 95% CIs)")+theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p
dev.off()





##  OLD with FIS all italic
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  labs(x="Region",y=expression(bold(paste("Inbreeding  ( ",italic(F[i][s]),paste("  ±  95% CIs)"))))) +
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p



## formatted correctly
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_95CI_per_region_419NDIV_1367SNPS_hwe44pops_FISperpond0.1_symbolsforamtted.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  labs(x="Region",y=expression(bold(paste("Inbreeding  ( ",italic(F)[is],paste(" )  ±  95% CIs"))))) +
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p
dev.off()


## new COL
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_95CI_per_region_419NDIV_1367SNPS_hwe44pops_FISperpond0.1_neworange.png",  bg = "transparent",width =200, height = 180, units = c("mm"))
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 8)+
  theme_classic() +
  scale_color_manual(values=c( "#ab4e03", "#fca45d","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  labs(x="Region",y=expression(bold(paste("Inbreeding  (",italic(F)[plain(IS)],paste(")")))))+theme(axis.title = element_text(face = "bold",size=20), axis.text = element_text(size = 18))
p<-p+theme(legend.position = "none",axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
p
dev.off()



####################
#### MAF VS Ho ######
#######################

genid.toad<-dartR::gl2gi(gl.toad, probar = FALSE, verbose = NULL)

pop(genid.toad)<-pop.data$onecluster

#### heirfstat conversion from genind ########
hierf.toad<-genind2hierfstat(genid.toad,pop=NULL)

#has pop
hierf.toad$pop


hierf.toad.basic<-basic.stats(hierf.toad)
hierf.toad.basic.perloc<-as.data.frame(hierf.toad.basic$perloc)

hierf.toad.basic.perloc$loci <- rownames(hierf.toad.basic.perloc)
head(hierf.toad.basic.perloc)
hierf.toad.basic.perloc$SNP<-hierf.toad.basic.perloc$loci 
#rownames(hierf.toad.basic.perloc) <- NULL

#edit loci names to match how formatted in maf
#hierf.toad.basic.perloc$loci <- sub( 'X', '', sub( '\\.', '\\/', sub( '\\.', '-', sub( '\\.', ':', hierf.toad.basic.perloc$loci))))
#head(hierf.toad.basic.perloc)


#minor allele freq
mafvcf<-maf(filtered.VCF,element=2)



mafvcf.df<-as.data.frame(mafvcf)

head(mafvcf.df)
dim()



# edit maf loci names to match row names of hierfstat

mafvcf.df$SNP <- rownames(mafvcf.df)
head(mafvcf.df)


# remove + and - from name
mafvcf.df$SNP<-gsub("\\+", "", mafvcf.df$SNP)
head(mafvcf.df)

mafvcf.df$SNP<-gsub("\\-", "", mafvcf.df$SNP)
head(mafvcf.df)

mafvcf.df$SNP<-paste("X", mafvcf.df$SNP, sep="")
head(mafvcf.df)


mafvcf.df$SNP<-gsub(":", ".", mafvcf.df$SNP)
head(mafvcf.df)


mafvcf.df$SNP<-paste(mafvcf.df$SNP, sep="",".")
head(mafvcf.df)



head(hierf.toad.basic.perloc)


## cind

maf_Ho<-cbind(hierf.toad.basic.perloc,mafvcf.df)
head(maf_Ho)

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAF_vs_Ho_plot_766INDIV_3496SNPS.png", width = 8, height = 6.5, units = 'in', res = 300)
plot(maf_Ho$Frequency,maf_Ho$Ho, xlab="Minor allele frequency",
     ylab="Observed heterozygosity (Ho)")
dev.off()







####################################
######### FIS PER LOCUS ############
#########################

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fis_perlocus_plot_766INDIV_3496SNPS_title.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(maf_Ho$Fis,breaks=seq(-1,1,l=50), xlab="Fis",
     ylab="Frequency",main ="Fis per locus no max obs het filter")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fis_perlocus_plot_766INDIV_3496SNPS.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(maf_Ho$Fis,breaks=seq(-1,1,l=50), xlab="Fis",
     ylab="Frequency",main ="")
dev.off()

#count nymber loci that have FIS under 0
sum(maf_Ho$Fis<0) #177

177/602
[1] 0.8428571

# 84% of loci have Fis below zero




###### what if I split into regions before calc across loci? 




#############################
##### obs het vs missingness #########
#######################################
hierf.toadobs<-basic.stats(hierf.toad)


hierf.toadobs$Ho









###############################
####### He vs Ho FIS ##########
#############################

##################################
###########  Ext Het m- max obs het as one ####################
############################


#### split matrix into 3 regions ###########

pop(gl.toad.nosibs.hwe.FIS) <- pop.data$fourclusters

gl.toad.nosibs.hwe.FISbypop<-seppop(gl.toad.nosibs.hwe.FIS)


gl.toad.nosibs.hwe.FISVanIsland<-gl.toad.nosibs.hwe.FISbypop$VanIsland
gl.toad.nosibs.hwe.FISLowerMain<-gl.toad.nosibs.hwe.FISbypop$LowerMain
gl.toad.nosibs.hwe.FISHaidaGwai<-gl.toad.nosibs.hwe.FISbypop$HaidaGwai
gl.toad.nosibs.hwe.FISNorthwest<-gl.toad.nosibs.hwe.FISbypop$Northwest




### convert genlight to genind objects 
genind.toad.nosibs.hwe.FISVanIsland<-dartR::gl2gi(gl.toad.nosibs.hwe.FISVanIsland, probar = FALSE, verbose = NULL)


""


genind.toad.nosibs.hwe.FISLowerMain<-dartR::gl2gi(gl.toad.nosibs.hwe.FISLowerMain, probar = FALSE, verbose = NULL)

""


genind.toad.nosibs.hwe.FISHaidaGwai<-dartR::gl2gi(gl.toad.nosibs.hwe.FISHaidaGwai, probar = FALSE, verbose = NULL)

''
genind.toad.nosibs.hwe.FISNorthwest<-dartR::gl2gi(gl.toad.nosibs.hwe.FISNorthwest, probar = FALSE, verbose = NULL)


############### summary on each one
genind.toad.nosibs.hwe.FISVanIsland_sum<-adegenet::summary(genind.toad.nosibs.hwe.FISVanIsland)

genind.toad.nosibs.hwe.FISLowerMain_sum<-adegenet::summary(genind.toad.nosibs.hwe.FISLowerMain)

genind.toad.nosibs.hwe.FISHaidaGwai_sum<-adegenet::summary(genind.toad.nosibs.hwe.FISHaidaGwai)

genind.toad.nosibs.hwe.FISNorthwest_sum<-adegenet::summary(genind.toad.nosibs.hwe.FISNorthwest)


###### Hexp per region ############

## rename
divVanIsland<-genind.toad.nosibs.hwe.FISVanIsland_sum
divLowerMain<-genind.toad.nosibs.hwe.FISLowerMain_sum
divHaidaGwai<-genind.toad.nosibs.hwe.FISHaidaGwai_sum
divNorthwest<-genind.toad.nosibs.hwe.FISNorthwest_sum


divLowerMain$LowerMain_Hexp<-divLowerMain$Hexp
divVanIsland$VanIsland_Hexp<-divVanIsland$Hexp
divHaidaGwai$HaidaGwai_Hexp<-divHaidaGwai$Hexp
divNorthwest$Northwest_Hexp<-divNorthwest$Hexp


divallfour<-cbind(divLowerMain$LowerMain_Hexp,divVanIsland$VanIsland_Hexp,divHaidaGwai$HaidaGwai_Hexp,divNorthwest$Northwest_Hexp)
head(divallfour)

## rename column names
names(divallfour)
colnames(divallfour) <- c("LowerMain_Hexp", "VanIsland_Hexp","HaidaGwai_Hexp","Northwest_Hexp")    # Applying colnames
head(divallfour)
str(divallfour)

# turn to dataframe
dfdivallfour<-as.data.frame(divallfour)
head(dfdivallfour)


### Hexp means + stdev ######
VanIsland_mean<-mean(dfdivallfour$VanIsland_Hexp) #
VanIsland_sd<-sd(dfdivallfour$VanIsland_Hexp) #

LowerMain_mean<-mean(dfdivallfour$LowerMain_Hexp) #
LowerMain_sd<-sd(dfdivallfour$LowerMain_Hexp) #

HaidaGwai_mean<-mean(dfdivallfour$HaidaGwai_Hexp) #
HaidaGwai_sd<-sd(dfdivallfour$HaidaGwai_Hexp) 


Northwest_mean<-mean(dfdivallfour$Northwest_Hexp) #
Northwest_sd<-sd(dfdivallfour$Northwest_Hexp) 

### save output
Hetsum<-cbind(VanIsland_mean, VanIsland_sd, LowerMain_mean,LowerMain_sd,HaidaGwai_mean,HaidaGwai_sd,Northwest_mean,Northwest_sd)

#write.table(Hetsum[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exp_Het_summary_3regions_601INDIV_173SNPS_maxhetobsasone_nohwe.txt")



dfdivallfour$locus<-row.names(dfdivallfour)
head(dfdivallfour)


colnames(dfdivallfour) <- c("LowerMain", "VanIsland","HaidaGwai","Northwest","locus")  
head(dfdivallfour)

### turn wide to long
df_long <- tidyr::gather(dfdivallfour,
                         key = Region,
                         value = Expect_Het,
                         LowerMain ,  VanIsland ,  HaidaGwai,Northwest)

head(df_long)
dim(df_long)

df_long_Hexp<-df_long
head(df_long_Hexp)


########## H obs #####

#### genind dartR from genligght ####
genind.toad.nosibs.hwe.FIS<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS, probar = FALSE, verbose = NULL)
genind.toad.nosibs.hwe.FIS
""

pop(genind.toad.nosibs.hwe.FIS)<-pop.data$fourclusters
pop(genind.toad.nosibs.hwe.FIS)

### check ###


#### heirfstat conversion from genind ########
hierf.toad.nosibs.hwe.FIS<-genind2hierfstat(genind.toad.nosibs.hwe.FIS)

#hierf.toad.nosibs.hwe.FIS<-genind2hierfstat(genind.toad.nosibs.hwe.FIS,pop=TRUE)
#hierf.toad.nosibs.hwe.FIS$pop<-pop.data$pop

hierf.toad.nosibs.hwe.FISbasic<-basic.stats(hierf.toad.nosibs.hwe.FIS)

hierf.toad.nosibs.hwe.FISbasic$Ho

head(hierf.toad.nosibs.hwe.FISbasic$Ho)

Hodf<-as.data.frame(hierf.toad.nosibs.hwe.FISbasic$Ho)

Hodf$locus<-row.names(Hodf)

head(Hodf)

Hodf$locus<- sub( 'X', '', sub( '\\.', '\\/', sub( '\\.', '-', sub( '\\.', ':', Hodf$locus))))

head(Hodf)

## rename 
colnames(Hodf) <- c("VanIsland_Hobs", "LowerMain_Hobs","HaidaGwai_Hobs","Northwest_Hobs","locus")  
head(Hodf)


colnames(Hodf) <- c("VanIsland", "LowerMain","HaidaGwai","Northwest","locus")  
head(Hodf)

### turn wide to long
df_long_Hobs <- tidyr::gather(Hodf,
                              key = Region,
                              value = Obs_Het,
                              VanIsland ,  LowerMain ,  HaidaGwai,Northwest)

head(df_long_Hobs)
dim(df_long_Hobs)


######### join Hobs and Hexp ####

Ho_He_df<-merge(x = df_long_Hobs, y = df_long_Hexp, by = c("locus","Region"), all = TRUE)

dim(Ho_He_df) # 360

head(Ho_He_df)
#View(Ho_He_df)

plot(Ho_He_df$Expect_Het,Ho_He_df$Obs_Het)

## rename regions

labels <- c(VanIsland = "Vancouver Island", LowerMain = "Lower mainland",HaidaGwai="Haida Gwaii",Northwest= "Northwest BC")

p<-ggplot(Ho_He_df, aes(x=Expect_Het, y=Obs_Het,fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_point()+
  facet_grid(Region ~ .,labeller=labeller(Region = labels))+
  theme_bw()+
  scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                     name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
xlab("Expected heterozygosity") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
# scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
# ylim(0,1)
p<-p+theme(legend.position = "none")
p



## save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/He_vs_Ho_per_region_419INDIV_1367SNPS_maxhetobs0.6asone.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(Ho_He_df, aes(x=Expect_Het, y=Obs_Het,fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_point()+
  facet_grid(Region ~ .,labeller=labeller(Region = labels))+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Expected heterozygosity") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
# scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
# ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()



## just Ho
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Ho_per_region_419INDIV_1367SNPS_maxhetobs0.6asone.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(Ho_He_df, aes(x=Region, y=Obs_Het,fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()

# just He
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/He_per_region_419INDIV_1367SNPS_maxhetobs0.6asone.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(Ho_He_df, aes(x=Region, y=Expect_Het,fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Expected heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()







##################################
###########  Ext Het & Obs het FIS ####################
############################


#### split matrix into 3 regions ###########

pop(gl.toad.nosibs.hwe.FIS) <- pop.data$fourclusters

gl.toad.nosibs.hwe.FIS.bypop<-seppop(gl.toad.nosibs.hwe.FIS)


gl.toad.nosibs.hwe.FIS.VanIsland<-gl.toad.nosibs.hwe.FIS.bypop$VanIsland


" /// GENLIGHT OBJECT /////////

 // 33 genotypes,  1503 binary SNPs, size: 286.8 Kb
 5026 (20.81 %) missing data

 // Basic content
   @gen: list of 33 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  33 individual labels
   @loc.names:  1503 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 33-33)
   @other: a list containing: elements without names "



gl.toad.nosibs.hwe.FIS.LowerMain<-gl.toad.nosibs.hwe.FIS.bypop$LowerMain

"/// GENLIGHT OBJECT /////////

 // 71 genotypes,  1503 binary SNPs, size: 374.8 Kb
 10519 (20.24 %) missing data

 // Basic content
   @gen: list of 71 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  71 individual labels
   @loc.names:  1503 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 71-71)
   @other: a list containing: elements without names 
"

gl.toad.nosibs.hwe.FIS.HaidaGwai<-gl.toad.nosibs.hwe.FIS.bypop$HaidaGwai

" /// GENLIGHT OBJECT /////////

 // 267 genotypes,  1503 binary SNPs, size: 860.7 Kb
 46991 (24.04 %) missing data

 // Basic content
   @gen: list of 267 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  267 individual labels
   @loc.names:  1503 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 267-267)
   @other: a list containing: elements without names "

### convert genlight to genind objects 
gl.toad.nosibs.hwe.FIS.VanIsland<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS.VanIsland, probar = FALSE, verbose = NULL)


""


gl.toad.nosibs.hwe.FIS.LowerMain<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS.LowerMain, probar = FALSE, verbose = NULL)

""


gl.toad.nosibs.hwe.FIS.HaidaGwai<-dartR::gl2gi(gl.toad.nosibs.hwe.FIS.HaidaGwai, probar = FALSE, verbose = NULL)

''


############### summary on each one
gl.toad.nosibs.hwe.FIS.VanIsland_sum<-adegenet::summary(gl.toad.nosibs.hwe.FIS.VanIsland)

gl.toad.nosibs.hwe.FIS.LowerMain_sum<-adegenet::summary(gl.toad.nosibs.hwe.FIS.LowerMain)

gl.toad.nosibs.hwe.FIS.HaidaGwai_sum<-adegenet::summary(gl.toad.nosibs.hwe.FIS.HaidaGwai)


## rename
divVanIsland<-gl.toad.nosibs.hwe.FIS.VanIsland_sum
divLowerMain<-gl.toad.nosibs.hwe.FIS.LowerMain_sum
divHaidaGwai<-gl.toad.nosibs.hwe.FIS.HaidaGwai_sum

###############################################
############  Obs Het per region ############
##########################################



hist(divVanIsland$Hobs)

VanIslandHobs<-divVanIsland$Hobs
LowerMainHobs<-divLowerMain$Hobs
HaidaGwaiHobs<-divHaidaGwai$Hobs

str(HaidaGwaiHobs)

lsNames=c("VanIslandHobs","LowerMainHobs","HaidaGwaiHobs")


## creates list for each locus
do.call(mapply, c(FUN=c, sapply(lsNames, as.symbol), SIMPLIFY=FALSE))

""


### just join as columns - cbind
HaidaGwaiHobs<-as.data.frame(HaidaGwaiHobs)
LowerMainHobs<-as.data.frame(LowerMainHobs)
VanIslandHobs<-as.data.frame(VanIslandHobs)

Hobsdf <- data.frame(cbind(LowerMainHobs,VanIslandHobs,HaidaGwaiHobs))

head(Hobsdf)

"  "


### turn wide to long
df_longHobs <- tidyr::gather(Hobsdf,
                             key = Region,
                             value = Obs_Het,
                             LowerMainHobs ,  VanIslandHobs ,  HaidaGwaiHobs)

head(df_longHobs)




### plot ###

ggplot(df_longHobs, aes(x=Region, y=Obs_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)

## save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Boxplot_obs_Het_per_region_419NDIV_1367SNPS_hwe44pops_FISperpond0.1.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(df_longHobs, aes(x=Region, y=Obs_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()

## save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Boxplot_obs_Het_per_region_419NDIV_1367SNPS_hwe44pops_FISperpond0.1_noylim.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(df_longHobs, aes(x=Region, y=Obs_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
  #ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()


##################################
###########  Ext Het FIS filter ####################
############################



###### Hexp per region ############

divLowerMain$LowerMain_Hexp<-divLowerMain$Hexp
divVanIsland$VanIsland_Hexp<-divVanIsland$Hexp
divHaidaGwai$HaidaGwai_Hexp<-divHaidaGwai$Hexp

divallthree<-cbind(divLowerMain$LowerMain_Hexp,divVanIsland$VanIsland_Hexp,divHaidaGwai$HaidaGwai_Hexp)
head(divallthree)

## rename column names
names(divallthree)
colnames(divallthree) <- c("LowerMain_Hexp", "VanIsland_Hexp","HaidaGwai_Hexp")    # Applying colnames
head(divallthree)
str(divallthree)

# turn to dataframe
dfdivallthree<-as.data.frame(divallthree)
head(dfdivallthree)


### Hexp means + stdev ######
VanIsland_mean<-mean(dfdivallthree$VanIsland_Hexp) #
VanIsland_sd<-sd(dfdivallthree$VanIsland_Hexp) #

LowerMain_mean<-mean(dfdivallthree$LowerMain_Hexp) #
LowerMain_sd<-sd(dfdivallthree$LowerMain_Hexp) #

HaidaGwai_mean<-mean(dfdivallthree$HaidaGwai_Hexp) #
HaidaGwai_sd<-sd(dfdivallthree$HaidaGwai_Hexp) 


### save output
Hetsum<-cbind(VanIsland_mean, VanIsland_sd, LowerMain_mean,LowerMain_sd,HaidaGwai_mean,HaidaGwai_sd)

write.table(Hetsum[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exp_Het_summary_3regions_419NDIV_1367SNPS_hwe44pops_FISperpond0.1.txt")





### turn wide to long
df_long <- tidyr::gather(dfdivallthree,
                         key = Region,
                         value = Expect_Het,
                         LowerMain_Hexp ,  VanIsland_Hexp ,  HaidaGwai_Hexp)

head(df_long)







### plot ###

p<-ggplot(df_long, aes(x=Region, y=Expect_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Expected heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p


## save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Boxplot_expect_Het_per_region_419NDIV_1367SNPS_hwe44pops_FISperpond0.1.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(df_long, aes(x=Region, y=Expect_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Expected heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()

## save no ylim
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Boxplot_expect_Het_per_region_419NDIV_1367SNPS_hwe44pops_FISperpond0.1_noylim.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(df_long, aes(x=Region, y=Expect_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Expected heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
#ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()





# with SE and 95% CI ###########


exphet_SE<-group_by(df_long, Region) %>%
  summarise_each(funs(mean=mean(Expect_Het),n=n(),sd=sd(Expect_Het),se=sd(.)/sqrt(n())),Expect_Het)


exphet_SEdf<-as.data.frame(exphet_SE)

exphet_SEdf$Exphet_highCI<-exphet_SEdf$mean+1.96*(exphet_SEdf$sd/sqrt(exphet_SEdf$n))
exphet_SEdf$Exphet_lowCI<-exphet_SEdf$mean-1.96*(exphet_SEdf$sd/sqrt(exphet_SEdf$n))

exphet_SEdf$Exphet_highSE<-exphet_SEdf$mean+exphet_SEdf$se

exphet_SEdf$Exphet_lowSE<-exphet_SEdf$mean-exphet_SEdf$se

exphet_SEdf



## reorder
exphet_SEdf$Region<- factor(exphet_SEdf$Region, levels = c("VanIsland_Hexp", "LowerMain_Hexp", "HaidaGwai_Hexp"))
exphet_SEdf<- droplevels(exphet_SEdf)

write.csv(exphet_SEdf,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exphet_mean_sd_SE_CI_MOE_419NDIV_1367SNPS_hwe44pops_FISperpond0.1.csv")



### NEW COL
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exphet_95pctCI_MOE_419NDIV_1367SNPS_hwe44pops_FISperpond0.1_orange_correctCI.pdf",  bg = "transparent",width =200, height = 180, units = c("mm"))
p<-ggplot(exphet_SEdf, aes(x = Region, color = Region, y=mean)) +
  geom_point(size = 8)+
  theme_classic() +
  scale_color_manual(values=c("#ab4e03", "#fca45d","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  
  geom_errorbar(aes(ymax = Exphet_highCI, ymin = Exphet_lowCI),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Expected Heterozygosity")+theme(axis.title = element_text(face = "bold",size=20), axis.text = element_text(size = 18))
p<-p+theme(legend.position = "none",axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
p
dev.off()


#############################
##### missingness #########
#######################################


#### genind dartR from genligght ####
genid.toadsub.hwe.nolas<-dartR::gl2gi(gl.toad.subset.nolas.hwe, probar = FALSE, verbose = NULL)

#set pop as region
pop(genid.toadsub.hwe.nolas)<-pop.data$fourclusters


#has pop
genid.toadsub.hwe.nolas@pop

#### heirfstat conversion from genind ########
hierf.toadsub.hwe.nolas<-genind2hierfstat(genid.toadsub.hwe.nolas,pop=NULL)



basicstat4clusters<-basic.stats(hierf.toadsub.hwe.nolas,diploid=TRUE)
hierf.toad.ind<-basicstat4clusters$n.ind.samp


### count number indiv per region

numVanIsland<-sum(pop.data$fourclusters=="VanIsland") #54
numLowerMain<-sum(pop.data$fourclusters=="LowerMain") #104
numHaidaGwai<-sum(pop.data$fourclusters=="HaidaGwai") #243


hierf.toad.ind




dfhierf.toad.ind<-as.data.frame(hierf.toad.ind)

dfhierf.toad.ind$VanIsland_pct<-dfhierf.toad.ind$VanIsland/numVanIsland
dfhierf.toad.ind$LowerMain_pct<-dfhierf.toad.ind$LowerMain/numLowerMain
dfhierf.toad.ind$HaidaGwai_pct<-dfhierf.toad.ind$HaidaGwai/numHaidaGwai



dfhierf.toad.ind$VanIsland_pctmissing<-(1-(dfhierf.toad.ind$VanIsland/numVanIsland))*100
dfhierf.toad.ind$LowerMain_pctmissing<-(1-(dfhierf.toad.ind$LowerMain/numLowerMain))*100
dfhierf.toad.ind$HaidaGwai_pctmissing<-(1-(dfhierf.toad.ind$HaidaGwai/numHaidaGwai))*100


head(dfhierf.toad.ind)



hist(dfhierf.toad.ind$HaidaGwai_pctmissing)

hist(dfhierf.toad.ind$VanIsland_pctmissing)

hist(dfhierf.toad.ind$LowerMain_pctmissing)

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/VanIsland_LowerMain_HaidaGwai_percentmissingloci_xlim_ylim_601INDIV_601SNPS.png", width = 10, height = 6.5, units = 'in', res = 600)
par(mfrow=c(1,3))
hist(dfhierf.toad.ind$HaidaGwai_pctmissing,xlim = c(0,100),ylim = c(0,1600))
hist(dfhierf.toad.ind$VanIsland_pctmissing,xlim = c(0,100),ylim = c(0,1600))
hist(dfhierf.toad.ind$LowerMain_pctmissing,xlim = c(0,100),ylim = c(0,1600))
dev.off()



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/VanIsland_LowerMain_HaidaGwai_percentmissingloci_xlim_ylim_40pct_601INDIV_601SNPS.png", width = 10, height = 6.5, units = 'in', res = 600)
par(mfrow=c(1,3))
hist(dfhierf.toad.ind$HaidaGwai_pctmissing,xlim = c(0,100),ylim = c(0,1600))
abline(v=40,col="blue")
hist(dfhierf.toad.ind$VanIsland_pctmissing,xlim = c(0,100),ylim = c(0,1600))
abline(v=40,col="blue")
hist(dfhierf.toad.ind$LowerMain_pctmissing,xlim = c(0,100),ylim = c(0,1600))
abline(v=40,col="blue")
dev.off()










###############################
####### He vs Ho ##########
#############################

##################################
###########  Ext Het m- max obs het as one ####################
############################


#### split matrix into 3 regions ###########

pop(gl.toad.subset.nolas) <- pop.data$fourclusters

gl.toad.subset.nolasbypop<-seppop(gl.toad.subset.nolas)


gl.toad.subset.nolasVanIsland<-gl.toad.subset.nolasbypop$VanIsland
gl.toad.subset.nolasLowerMain<-gl.toad.subset.nolasbypop$LowerMain
gl.toad.subset.nolasHaidaGwai<-gl.toad.subset.nolasbypop$HaidaGwai



### convert genlight to genind objects 
gl.toad.subset.nolasVanIsland<-dartR::gl2gi(gl.toad.subset.nolasVanIsland, probar = FALSE, verbose = NULL)


""


gl.toad.subset.nolasLowerMain<-dartR::gl2gi(gl.toad.subset.nolasLowerMain, probar = FALSE, verbose = NULL)

""


gl.toad.subset.nolasHaidaGwai<-dartR::gl2gi(gl.toad.subset.nolasHaidaGwai, probar = FALSE, verbose = NULL)

''


############### summary on each one
gl.toad.subset.nolasVanIsland_sum<-adegenet::summary(gl.toad.subset.nolasVanIsland)

gl.toad.subset.nolasLowerMain_sum<-adegenet::summary(gl.toad.subset.nolasLowerMain)

gl.toad.subset.nolasHaidaGwai_sum<-adegenet::summary(gl.toad.subset.nolasHaidaGwai)



###### Hexp per region ############

## rename
divVanIsland<-gl.toad.subset.nolasVanIsland_sum
divLowerMain<-gl.toad.subset.nolasLowerMain_sum
divHaidaGwai<-gl.toad.subset.nolasHaidaGwai_sum

divLowerMain$LowerMain_Hexp<-divLowerMain$Hexp
divVanIsland$VanIsland_Hexp<-divVanIsland$Hexp
divHaidaGwai$HaidaGwai_Hexp<-divHaidaGwai$Hexp

divallthree<-cbind(divLowerMain$LowerMain_Hexp,divVanIsland$VanIsland_Hexp,divHaidaGwai$HaidaGwai_Hexp)
head(divallthree)

## rename column names
names(divallthree)
colnames(divallthree) <- c("LowerMain_Hexp", "VanIsland_Hexp","HaidaGwai_Hexp")    # Applying colnames
head(divallthree)
str(divallthree)

# turn to dataframe
dfdivallthree<-as.data.frame(divallthree)
head(dfdivallthree)


### Hexp means + stdev ######
VanIsland_mean<-mean(dfdivallthree$VanIsland_Hexp) #
VanIsland_sd<-sd(dfdivallthree$VanIsland_Hexp) #

LowerMain_mean<-mean(dfdivallthree$LowerMain_Hexp) #
LowerMain_sd<-sd(dfdivallthree$LowerMain_Hexp) #

HaidaGwai_mean<-mean(dfdivallthree$HaidaGwai_Hexp) #
HaidaGwai_sd<-sd(dfdivallthree$HaidaGwai_Hexp) 


### save output
Hetsum<-cbind(VanIsland_mean, VanIsland_sd, LowerMain_mean,LowerMain_sd,HaidaGwai_mean,HaidaGwai_sd)

#write.table(Hetsum[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exp_Het_summary_3regions_601INDIV_173SNPS_maxhetobsasone_nohwe.txt")



dfdivallthree$locus<-row.names(dfdivallthree)
head(dfdivallthree)


colnames(dfdivallthree) <- c("LowerMain", "VanIsland","HaidaGwai","locus")  
head(dfdivallthree)

### turn wide to long
df_long <- tidyr::gather(dfdivallthree,
                         key = Region,
                         value = Expect_Het,
                         LowerMain ,  VanIsland ,  HaidaGwai)

head(df_long)
dim(df_long)

df_long_Hexp<-df_long
head(df_long_Hexp)


########## H obs #####

#### genind dartR from genligght ####
genid.toadsub.nolas<-dartR::gl2gi(gl.toad.subset.nolas, probar = FALSE, verbose = NULL)
genid.toadsub.nolas
""

pop(genid.toadsub.nolas)<-pop.data$fourclusters
pop(genid.toadsub.nolas)

### check ###


#### heirfstat conversion from genind ########
hierf.toad.nolas<-genind2hierfstat(genid.toadsub.nolas)

#hierf.toad.nolas<-genind2hierfstat(genid.toadsub.nolas,pop=TRUE)
#hierf.toad.nolas$pop<-pop.data$pop

hierf.toad.nolasbasic<-basic.stats(hierf.toad.nolas)

hierf.toad.nolasbasic$Ho

head(hierf.toad.nolasbasic$Ho)

Hodf<-as.data.frame(hierf.toad.nolasbasic$Ho)

Hodf$locus<-row.names(Hodf)

head(Hodf)

Hodf$locus<- sub( 'X', '', sub( '\\.', '\\/', sub( '\\.', '-', sub( '\\.', ':', Hodf$locus))))

head(Hodf)

## rename 
colnames(Hodf) <- c("VanIsland_Hobs", "LowerMain_Hobs","HaidaGwai_Hobs","locus")  
head(Hodf)


colnames(Hodf) <- c("VanIsland", "LowerMain","HaidaGwai","locus")  
head(Hodf)

### turn wide to long
df_long_Hobs <- tidyr::gather(Hodf,
                              key = Region,
                              value = Obs_Het,
                              VanIsland ,  LowerMain ,  HaidaGwai)

head(df_long_Hobs)
dim(df_long_Hobs)


######### join Hobs and Hexp ####

Ho_He_df<-merge(x = df_long_Hobs, y = df_long_Hexp, by = c("locus","Region"), all = TRUE)

dim(Ho_He_df) # 360

head(Ho_He_df)
#View(Ho_He_df)

plot(Ho_He_df$Expect_Het,Ho_He_df$Obs_Het)

## rename regions

labels <- c(VanIsland = "Vancouver Island", LowerMain = "Lower mainland",HaidaGwai="Haida Gwaii")

p<-ggplot(Ho_He_df, aes(x=Expect_Het, y=Obs_Het,fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_point()+
  facet_grid(Region ~ .,labeller=labeller(Region = labels))+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Expected heterozygosity") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
# scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
# ylim(0,1)
p<-p+theme(legend.position = "none")
p



## save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/He_vs_Ho_per_region_601INDIV_173SNPS_maxhetobs0.6asone.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(Ho_He_df, aes(x=Expect_Het, y=Obs_Het,fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_point()+
  facet_grid(Region ~ .,labeller=labeller(Region = labels))+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Expected heterozygosity") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
# scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
# ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()









































############# no additional filters #######


##################################
###########  Ext Het no additional filters ####################
############################


#### split matrix into 3 regions ###########

pop(gl.toad) <- pop.data$fourclusters

gl.toad.bypop<-seppop(gl.toad)


gl.toad.VanIsland<-gl.toad.bypop$VanIsland
gl.toad.LowerMain<-gl.toad.bypop$LowerMain
gl.toad.HaidaGwai<-gl.toad.bypop$HaidaGwai



### convert genlight to genind objects 
gl.toad.VanIsland<-dartR::gl2gi(gl.toad.VanIsland, probar = FALSE, verbose = NULL)


""


gl.toad.LowerMain<-dartR::gl2gi(gl.toad.LowerMain, probar = FALSE, verbose = NULL)

""


gl.toad.HaidaGwai<-dartR::gl2gi(gl.toad.HaidaGwai, probar = FALSE, verbose = NULL)

''


############### summary on each one
gl.toad.VanIsland_sum<-adegenet::summary(gl.toad.VanIsland)

gl.toad.LowerMain_sum<-adegenet::summary(gl.toad.LowerMain)

gl.toad.HaidaGwai_sum<-adegenet::summary(gl.toad.HaidaGwai)



###### Hexp per region ############

## rename
divVanIsland<-gl.toad.VanIsland_sum
divLowerMain<-gl.toad.LowerMain_sum
divHaidaGwai<-gl.toad.HaidaGwai_sum

divLowerMain$LowerMain_Hexp<-divLowerMain$Hexp
divVanIsland$VanIsland_Hexp<-divVanIsland$Hexp
divHaidaGwai$HaidaGwai_Hexp<-divHaidaGwai$Hexp

divallthree<-cbind(divLowerMain$LowerMain_Hexp,divVanIsland$VanIsland_Hexp,divHaidaGwai$HaidaGwai_Hexp)
head(divallthree)

## rename column names
names(divallthree)
colnames(divallthree) <- c("LowerMain_Hexp", "VanIsland_Hexp","HaidaGwai_Hexp")    # Applying colnames
head(divallthree)
str(divallthree)

# turn to dataframe
dfdivallthree<-as.data.frame(divallthree)
head(dfdivallthree)


### Hexp means + stdev ######
VanIsland_mean<-mean(dfdivallthree$VanIsland_Hexp) #
VanIsland_sd<-sd(dfdivallthree$VanIsland_Hexp) #

LowerMain_mean<-mean(dfdivallthree$LowerMain_Hexp) #
LowerMain_sd<-sd(dfdivallthree$LowerMain_Hexp) #

HaidaGwai_mean<-mean(dfdivallthree$HaidaGwai_Hexp) #
HaidaGwai_sd<-sd(dfdivallthree$HaidaGwai_Hexp) 


### save output
Hetsum<-cbind(VanIsland_mean, VanIsland_sd, LowerMain_mean,LowerMain_sd,HaidaGwai_mean,HaidaGwai_sd)

write.table(Hetsum[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exp_Het_summary_3regions_766INDIV_3496SNPS.txt")


### turn wide to long
df_long <- tidyr::gather(dfdivallthree,
                         key = Region,
                         value = Expect_Het,
                         LowerMain_Hexp ,  VanIsland_Hexp ,  HaidaGwai_Hexp)

head(df_long)


### plot ###

p<-ggplot(df_long, aes(x=Region, y=Expect_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Expected heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p


## save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Boxplot_expect_Het_per_region_766INDIV_3496SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(df_long, aes(x=Region, y=Expect_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Expected heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()

## save no ylim
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Boxplot_expect_Het_per_region_766INDIV_3496SNPS_noylim.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(df_long, aes(x=Region, y=Expect_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Expected heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
#ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()






############### basic stats to get FIS per locus averaged across pops ##################

#### genind dartR from genligght ####
genid.toad<-dartR::gl2gi(gl.toad, probar = FALSE, verbose = NULL)
genid.toad
""

pop(genid.toad)<-pop.data$pop
pop(genid.toad)

### check ###


#### heirfstat conversion from genind ########
hierf.toad<-genind2hierfstat(genid.toad)

#hierf.toad<-genind2hierfstat(genid.toad,pop=TRUE)
#hierf.toad$pop<-pop.data$pop

hierf.toadbasic<-basic.stats(hierf.toad)

hierf.toadbasic$Fis

head(hierf.toadbasic$Fis)

FISdataperpond<-hierf.toadbasic$Fis

FISdataperpondavg <- data.frame("SNP"=rownames(FISdataperpond), "Average"=rowMeans(FISdataperpond, na.rm=T))
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/FISperpond_plot_766INDIV_3496SNPS_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))
dev.off()

# snps I want to keep
nrow(FISdataperpondavg[which(FISdataperpondavg$Average<0.1 & FISdataperpondavg$Average>-0.1),]) 
#120

# snps I want to delete
nrow(FISdataperpondavg[which(FISdataperpondavg$Average>0.1 & FISdataperpondavg$Average>-0.1),]) 
#15

nrow(FISdataperpondavg) #602





##### per region ############

FISdataperponddf<-as.data.frame(FISdataperpond)

## HG ###
FISdataperpondHG<-FISdataperponddf[ , grepl( "HaidaGwai" , names( FISdataperponddf ) ) ]

FISdataperpondHGavg <- data.frame("SNP"=rownames(FISdataperpondHG), "Average"=rowMeans(FISdataperpondHG, na.rm=T))
hist(FISdataperpondHGavg$Average,breaks=seq(-1,1,l=100))

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/FISperpond_HG_plot_766INDIV_3496SNPS_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondHGavg$Average,breaks=seq(-1,1,l=100))
dev.off()

## VI ###
FISdataperpondVI<-FISdataperponddf[ , grepl( "VanIsland" , names( FISdataperponddf ) ) ]

FISdataperpondVIavg <- data.frame("SNP"=rownames(FISdataperpondVI), "Average"=rowMeans(FISdataperpondVI, na.rm=T))
hist(FISdataperpondVIavg$Average,breaks=seq(-1,1,l=100))

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/FISperpond_VI_plot_766INDIV_3496SNPS_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondVIavg$Average,breaks=seq(-1,1,l=100))
dev.off()

## LM ###
FISdataperpondLM<-FISdataperponddf[ , grepl( "LowerMain" , names( FISdataperponddf ) ) ]

FISdataperpondLMavg <- data.frame("SNP"=rownames(FISdataperpondLM), "Average"=rowMeans(FISdataperpondLM, na.rm=T))
hist(FISdataperpondLMavg$Average,breaks=seq(-1,1,l=100))

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/FISperpond_LM_plot_766INDIV_3496SNPS_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondLMavg$Average,breaks=seq(-1,1,l=100))
dev.off()


### per individual - doesn't work ######

hierf.toadindiv<-hierf.toad
hierf.toadindiv$pop<-pop.data$sample.id

hierf.toadindivbasic<-basic.stats(hierf.toadindiv)

hierf.toadindivbasic$Fis

head(hierf.toadbasic$Fis)

FISdataperpond<-hierf.toadbasic$Fis

FISdataperpondavg <- data.frame("SNP"=rownames(FISdataperpond), "Average"=rowMeans(FISdataperpond, na.rm=T))
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))



## compare av across all individuals - metapop

hierf.toad2<-genind2hierfstat(genid.toad,pop=TRUE)
#hierf.toad$pop<-pop.data$pop

hierf.toadbasic2<-basic.stats(hierf.toad2)

#one pop
head(hierf.toadbasic2$Fis)



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/FISacrossalllocusmetapop_plot_766INDIV_3496SNPS_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(hierf.toadbasic2$Fis,breaks=seq(-1,1,l=100))
dev.off()

# snps I want to keep
nrow(hierf.toadbasic2$Fis[which(hierf.toadbasic2$Fis<0.1 & hierf.toadbasic2$Fis>-0.1),]) 
#132

nrow(hierf.toadbasic2$Fis) #602






######################################
########## NJ TREE ##################
###############################################

library(ape)
tre <- nj(dist(as.matrix(gl.toad.nosibs.hwe.FIS)))
tre
plot(tre, typ="fan", cex=0.7)
#title("NJ tree of the 18 test samples filter 22")

#par(mar=c(5.1,4.1,4.1,2.1))



#colour 4clusters
pop(gl.toad.nosibs.hwe.FIS) <- pop.data$fourclusters

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/NJTREE/NJTREE_col4clusters_419NDIV_1367SNPS.png", width = 6.5, height = 6.5, units = 'in', res = 300)
plot(tre, typ="fan", show.tip=TRUE, cex=0.5)
tiplabels(pch=20, col=brewer.pal(10, "Paired")[gl.toad.nosibs.hwe.FIS$pop], cex=4)
#title("NJ tree of Western toads for the 20 test samples/n colour from subregion/n denovo test 20 samples rm sibs gapped 0.9 M3 rm sibs fitlering may 11th")
dev.off()






###################################
######## DAPC ###############
#########################




# https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf


# 1. find correct no. clusters with find.clusters()
# 2. run DAPC for each cluster that looks reasonable - using 200 PCs - using dapc()
# 3. find the optimal number of PCs using optim.a.score()
# 4. membership prob plots

#################################
######### find correct number clusters ##########
####################################

############################
####### 200 PCs, 3 groups

grp <- find.clusters(gl.toad.nosibs.hwe.FIS, max.n.clust=26)


# save variance explained PC plot
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_varianceexplainedbyPCA_findclusters_maxclusters26_419NDIV_1367SNPS.png")



# https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
# pick largest
#Choose the number PCs to retain (>=1): 
200

# save BIC plot

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_BIC_findclusters_maxclusters26_200PCs_419NDIV_1367SNPS.png")



#Choose the number of clusters (>=2): 
3

grp$grp

# saving cluster allocation
write.csv(grp$grp, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters26_200PCs_3clusters_419NDIV_1367SNPS.csv")



############################
####### 200 PCs, 2 groups #######

grp <- find.clusters(gl.toad.nosibs.hwe.FIS, max.n.clust=26)


#Choose the number PCs to retain (>=1): 
200


#Choose the number of clusters (>=2): 
2

grp$grp

# saving cluster allocation
write.csv(grp$grp, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters26_200PCs_2clusters_419NDIV_1367SNPS.csv")



############################
####### 200 PCs, 4 groups #######

grp <- find.clusters(gl.toad.nosibs.hwe.FIS, max.n.clust=26)


#Choose the number PCs to retain (>=1): 
200


#Choose the number of clusters (>=2): 
4

grp$grp

# saving cluster allocation
write.csv(grp$grp, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters26_200PCs_4clusters_419NDIV_1367SNPS.csv")




##############################
######## DAPC CLUSTER 4 ###########
###########################



clustass4clusters<-read.csv( file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters26_200PCs_4clusters_419NDIV_1367SNPS.csv")



clustass4clusters<-as.data.frame(clustass4clusters)

str(clustass4clusters)

clustass4clusters$sample.id<-clustass4clusters$X

head(clustass4clusters)

clustass4clusters$fourclusters200pca<-clustass4clusters$x

#JUST SAVE 2 COLS
clustass4clusterssub<-clustass4clusters[c(3:4)]

head(clustass4clusterssub)




### join w pop data 
pop.data4clusters<-left_join(pop.data,clustass4clusterssub)

head(pop.data4clusters)


# work out which pop corresponds to which cluster number
#View(pop.data4clusters)


pop.data4clusters$fourclusters200pcaNAMES<-pop.data4clusters$fourclusters200pca

pop.data4clusters$fourclusters200pcaNAMES[which(pop.data4clusters$fourclusters200pca=="1")] <- "SwBCmix1"

pop.data4clusters$fourclusters200pcaNAMES[which(pop.data4clusters$fourclusters200pca=="4")] <- "HaidaGwaii"

pop.data4clusters$fourclusters200pcaNAMES[which(pop.data4clusters$fourclusters200pca=="2")] <- "SwBCmix2"

pop.data4clusters$fourclusters200pcaNAMES[which(pop.data4clusters$fourclusters200pca=="3")] <- "Chilliwack"

head(pop.data4clusters)



pop(gl.toad.nosibs.hwe.FIS) <- pop.data4clusters$fourclusters200pcaNAMES




## dapc with max pcs and das
pnw.dapc <- dapc(gl.toad.nosibs.hwe.FIS, n.pca = 200, n.da = 10, parallel = require("parallel"))




#### check optimal number of PCs = 23
temp <- optim.a.score(pnw.dapc)


# run dapc with optimal no. pcs
pnw.dapc2 <- dapc(gl.toad.nosibs.hwe.FIS, n.pca = 23, n.da = 10, parallel = require("parallel"))


scatter(pnw.dapc2, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_419NDIV_1367SNPS_optimal_23PCaxes_defaultcols.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc2, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)
dev.off()



# reduce da axes from 10 to 4 to test - makes no difference
pnw.dapc3 <- dapc(gl.toad.nosibs.hwe.FIS, n.pca = 16, n.da = 4, parallel = require("parallel"))

scatter(pnw.dapc3, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)






labs <- c("Chikundal", "Gudal Lake","Gwaii Haanas","Northern Haida Gwaii")
cols=c("#58d19f","pink","#0072B2", "#D55E00" )


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_419NDIV_1367SNPS_optimal_23PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc2, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_419NDIV_1367SNPS_200PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()



### DAPC with 80% variation explained####

highdapc <- dapc(gl.toad.nosibs.hwe.FIS, parallel = require("parallel"), n.da = 10, pca.select = "percVar", perc.pca=80)



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_419NDIV_1367SNPS_80pctPCs_defaultcols.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(highdapc, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)
dev.off()


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_419NDIV_1367SNPS_80pctPCs.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(highdapc, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()





#Loci contributing to observed differences, threshold set arbitrarily?
set.seed(4)
contrib <- loadingplot(pnw.dapc2$var.contr, axis = 2, thres = 0.02, lab.jitter = 1)

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_loadings_col_4_clusters_419NDIV_1367SNPS_optimal_23PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
contrib <- loadingplot(pnw.dapc2$var.contr, axis = 2, thres = 0.02, lab.jitter = 1)
dev.off()




#### MEMBERSHIP PROP ######
#probability of population belonging to the pop it's assigned

set.seed(999)
pramx <- xvalDapc(tab(gl.toad.nosibs.hwe.FIS), pop(gl.toad.nosibs.hwe.FIS), parallel = "snow")
###--->40



compoplot(pnw.dapc2,col = brewer.pal(4, "Paired"), posi = 'top')


dapc.results <- as.data.frame(pnw.dapc2$posterior)
dapc.results$oldpop <- pop(gl.toad.nosibs.hwe.FIS)
dapc.results$indNames <- rownames(dapc.results)


dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample",
                            "Assigned_Pop","Posterior_membership_probability")


###### no sep areas 
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_4_clusters_419NDIV_1367SNPS_optimalPC_23PCaxes_nosepareas.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_blank(),axis.title.y = element_text( size = 10),axis.ticks.x = element_blank(),panel.background = element_blank())
t <- t + scale_fill_manual(values=c( "lightblue","#009E73","#7949a6","#9999CC"),
                           name="Assigned Pop K=4 cluster", breaks = c("Chilliwack","HaidaGwaii","SwBCmix1","SwBCmix2"), labels=c("Chilliwack","HaidaGwaii","SwBCmix1","SwBCmix2"))
t
dev.off()


# with sample id
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_4_clusters_419NDIV_1367SNPS_optimalPC_23PCaxes_nosepareas_sampleid.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_text(angle = 90, size = 2,vjust=0.2),axis.title.y = element_text( size = 8))
t <- t + scale_fill_manual(values=c( "lightblue","#009E73","#7949a6","#9999CC"),
                           name="Assigned Pop K=4 cluster", breaks = c("Chilliwack","HaidaGwaii","SwBCmix1","SwBCmix2"), labels=c("Chilliwack","HaidaGwaii","SwBCmix1","SwBCmix2"))
t
dev.off()




t <- t + scale_fill_manual(values=c("#009E73", "#9999CC",  "#7949a6"),
                           name="Assigned Pop K=3 cluster", breaks = c("Haida_Gwaii","SwBC1", "SwBC_Sea2Sky"), labels=c("Haida_Gwaii", "SwBC1", "SwBC_Sea2Sky"))
t
dev.off()



######## as four original pops 

pop(gl.toad.nosibs.hwe.FIS)<-pop.data4clusters$fourclusters200pcaNAMES




#make labs
supp.labs <- c( Gwaii_Haanas=" Gwaii Haanas (original pop)",Chikundal="Chikundal (original pops)",Northern_Haida_Gwaii="Northern Haida Gwaii (original pop)",Gudal_Lake ="Gudal Lake (original pop)")


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_4_clusters_419NDIV_1367SNPS_optimalPC_23PCaxes_4originalareas.png", width = 12, height = 6.5, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + facet_grid(~Original_Pop, scales = "free",  labeller = labeller(Original_Pop = supp.labs))
t <- t + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
t <- t + scale_fill_manual(values=c( "#0072B2","#58d19f","#D55E00","pink"),
                           name="Assigned Pop K4 cluster", breaks = c("Gwaii_Haanas","Chikundal","Northern_Haida_Gwaii","Gudal_Lake"), labels=c("Gwaii Haanas","Chikundal","Northern Haida Gwaii","Gudal Lake"))
t
dev.off()



## as two original pops....
pop(gl.toad.nosibs.hwe.FIS)<-pop.data4clusters$two_areas






#make labs
supp.labs <- c( Gwaii_Haanas=" Gwaii Haanas (original pop)",Northern_HaidaGwaii="Northern Haida Gwaii (original pop)")


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_4_clusters_419NDIV_1367SNPS_optimalPC_23PCaxes_2originalareas.png", width = 12, height = 6.5, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + facet_grid(~Original_Pop, scales = "free",  labeller = labeller(Original_Pop = supp.labs))
t <- t + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
t <- t + scale_fill_manual(values=c( "#0072B2","#58d19f","#D55E00","pink"),
                           name="Assigned Pop K4 cluster", breaks = c("Gwaii_Haanas","Chikundal","Northern_Haida_Gwaii","Gudal_Lake"), labels=c("Gwaii Haanas","Chikundal","Northern Haida Gwaii","Gudal Lake"))
t
dev.off()






##############################
######## DAPC CLUSTER 3 ###########
###########################



clustass3clusters<-read.csv( file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters26_200PCs_3clusters_419NDIV_1367SNPS.csv")
dim(clustass3clusters)
dim(pop.data)

clustass3clusters<-as.data.frame(clustass3clusters)

str(clustass3clusters)

clustass3clusters$sample.id<-clustass3clusters$X

head(clustass3clusters)

clustass3clusters$threeclusters200pca<-clustass3clusters$x

#JUST SAVE 2 COLS
clustass3clusterssub<-clustass3clusters[c(3:4)]

head(clustass3clusterssub)




### join w pop data 
pop.data3clusters<-left_join(pop.data,clustass3clusterssub)

head(pop.data3clusters)


# work out which pop corresponds to which cluster number
#View(pop.data3clusters)


pop.data3clusters$threeclusters200pcaNAMES<-pop.data3clusters$threeclusters200pca

pop.data3clusters$threeclusters200pcaNAMES[which(pop.data3clusters$threeclusters200pca=="2")] <- "Haida_Gwaii"

pop.data3clusters$threeclusters200pcaNAMES[which(pop.data3clusters$threeclusters200pca=="1")] <- "SwBC1"

pop.data3clusters$threeclusters200pcaNAMES[which(pop.data3clusters$threeclusters200pca=="3")] <- "SwBC_Sea2Sky"


head(pop.data3clusters)



pop(gl.toad.nosibs.hwe.FIS) <- pop.data3clusters$threeclusters200pcaNAMES




## dapc with max pcs and das
pnw.dapc <- dapc(gl.toad.nosibs.hwe.FIS, n.pca = 200, n.da = 10, parallel = require("parallel"))




#### check optimal number of PCs = 24
temp <- optim.a.score(pnw.dapc)


# run dapc with optimal no. pcs
pnw.dapc2 <- dapc(gl.toad.nosibs.hwe.FIS, n.pca = 24, n.da = 10, parallel = require("parallel"))


scatter(pnw.dapc2, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_3_clusters_419NDIV_1367SNPS_optimal_24PCaxes_defaultcols.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc2, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)
dev.off()




# reduce da axes from 10 to 4 to test - makes no difference
pnw.dapc3 <- dapc(gl.toad.nosibs.hwe.FIS, n.pca = 16, n.da = 4, parallel = require("parallel"))

scatter(pnw.dapc3, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)






labs <- c( "Gudal Lake","Gwaii Haanas","Northern Haida Gwaii")
cols=c("pink","#0072B2", "#D55E00" )



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_3_clusters_419NDIV_1367SNPS_optimal_24PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc2, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_3_clusters_419NDIV_1367SNPS_200PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc,col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()



## #DAPC with 80% variation explained####

highdapc <- dapc(gl.toad.nosibs.hwe.FIS, parallel = require("parallel"), n.da = 10, pca.select = "percVar", perc.pca=80)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_3_clusters_419NDIV_1367SNPS_80pctPCs.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(highdapc,col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()





#Loci contributing to observed differences, threshold set arbitrarily?
set.seed(4)
contrib <- loadingplot(pnw.dapc2$var.contr, axis = 2, thres = 0.02, lab.jitter = 1)

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_loadings_col_3_clusters_419NDIV_1367SNPS_optimal_24PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
contrib <- loadingplot(pnw.dapc2$var.contr, axis = 2, thres = 0.02, lab.jitter = 1)
dev.off()



#### MEMBERSHIP PROP ######
#probability of population belonging to the pop it's assigned


pop(gl.toad.nosibs.hwe.FIS)<-pop.data3clusters$threeclusters200pcaNAMES

set.seed(999)
pramx <- xvalDapc(tab(gl.toad.nosibs.hwe.FIS), pop(gl.toad.nosibs.hwe.FIS), parallel = "snow")
###--->40



compoplot(pnw.dapc2,col = brewer.pal(4, "Paired"), posi = 'top')



dapc.results <- as.data.frame(pnw.dapc2$posterior)
dapc.results$oldpop <- pop(gl.toad.nosibs.hwe.FIS)
dapc.results$indNames <- rownames(dapc.results)


dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample",
                            "Assigned_Pop","Posterior_membership_probability")

###### no sep areas 
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_3_clusters_419NDIV_1367SNPS_optimalPC_24PCaxes_nosepareas.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_blank(),axis.title.y = element_text( size = 10),axis.ticks.x = element_blank(),panel.background = element_blank())
t <- t + scale_fill_manual(values=c("#009E73", "#9999CC",  "#7949a6"),
                           name="Assigned Pop K=3 cluster", breaks = c("Haida_Gwaii","SwBC1", "SwBC_Sea2Sky"), labels=c("Haida_Gwaii", "SwBC1", "SwBC_Sea2Sky"))
t
dev.off()


# with sample id
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_3_clusters_419NDIV_1367SNPS_optimalPC_24PCaxes_nosepareas_sampleid.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_text(angle = 90, size = 3),axis.title.y = element_text( size = 8))
t <- t + scale_fill_manual(values=c("#009E73", "#9999CC",  "#7949a6"),
                           name="Assigned Pop K=3 cluster", breaks = c("Haida_Gwaii","SwBC1", "SwBC_Sea2Sky"), labels=c("Haida_Gwaii", "SwBC1", "SwBC_Sea2Sky"))
t
dev.off()


######## as four original pops 

pop(gl.toad.nosibs.hwe.FIS)<-pop.data4clusters$fourclusters200pcaNAMES



#make labs
supp.labs <- c( Gwaii_Haanas=" Gwaii Haanas (original pop)",Chikundal="Chikundal (original pops)",Northern_Haida_Gwaii="Northern Haida Gwaii (original pop)",Gudal_Lake ="Gudal Lake (original pop)")


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_4_clusters_419NDIV_1367SNPS_optimalPC_24PCaxes_4originalareas.png", width = 12, height = 6.5, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + facet_grid(~Original_Pop, scales = "free",  labeller = labeller(Original_Pop = supp.labs))
t <- t + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
t <- t + scale_fill_manual(values=c( "#0072B2","#58d19f","#D55E00","pink"),
                           name="Assigned Pop K4 cluster", breaks = c("Gwaii_Haanas","Chikundal","Northern_Haida_Gwaii","Gudal_Lake"), labels=c("Gwaii Haanas","Chikundal","Northern Haida Gwaii","Gudal Lake"))
t
dev.off()



## as two original pops....
pop(gl.toad.nosibs.hwe.FIS)<-pop.data3clusters$two_areas





#make labs
supp.labs <- c( Gwaii_Haanas=" Gwaii Haanas (original pop)",Northern_HaidaGwaii="Northern Haida Gwaii (original pop)")


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_3_clusters_419NDIV_1367SNPS_optimalPC_24PCaxes_2originalareas.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + facet_grid(~Original_Pop, scales = "free",  labeller = labeller(Original_Pop = supp.labs))
t <- t + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),axis.title.y = element_text( size = 8))
t <- t + scale_fill_manual(values=c( "#0072B2","#D55E00","pink"),
                           name="Assigned Pop K=3 cluster", breaks = c("Gwaii_Haanas","Northern_Haida_Gwaii","Gudal_Lake"), labels=c("Gwaii Haanas","Northern Haida Gwaii","Gudal Lake"))
t
dev.off()









##############################
######## DAPC CLUSTER 2 ###########
###########################



clustass2clusters<-read.csv( file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters26_200PCs_2clusters_419NDIV_1367SNPS.csv")



clustass2clusters<-as.data.frame(clustass2clusters)

str(clustass2clusters)

clustass2clusters$sample.id<-clustass2clusters$X

head(clustass2clusters)

clustass2clusters$twoclusters200pca<-clustass2clusters$x

#JUST SAVE 2 COLS
clustass2clusterssub<-clustass2clusters[c(3:4)]

head(clustass2clusterssub)




### join w pop data 
pop.data2clusters<-left_join(pop.data,clustass2clusterssub)

head(pop.data2clusters)


# work out which pop corresponds to which cluster number
#View(pop.data2clusters)


pop.data2clusters$twoclusters200pcaNAMES<-pop.data2clusters$twoclusters200pca


pop.data2clusters$twoclusters200pcaNAMES[which(pop.data2clusters$twoclusters200pca=="1")] <- "Southwest_BC"

pop.data2clusters$twoclusters200pcaNAMES[which(pop.data2clusters$twoclusters200pca=="2")] <- "Haida_Gwaii"


head(pop.data2clusters)



pop(gl.toad.nosibs.hwe.FIS) <- pop.data2clusters$twoclusters200pcaNAMES




## dapc with max pcs and das
pnw.dapc <- dapc(gl.toad.nosibs.hwe.FIS, n.pca = 200, n.da = 10, parallel = require("parallel"))


scatter(pnw.dapc, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)

#### check optimal number of PCs = 1
temp <- optim.a.score(pnw.dapc)


# run dapc with optimal no. pcs
pnw.dapc2 <- dapc(gl.toad.nosibs.hwe.FIS, n.pca = 1, n.da = 10, parallel = require("parallel"))

scatter(pnw.dapc2, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)

# reduce da axes from 10 to 4 to test - makes no difference
pnw.dapc3 <- dapc(gl.toad.nosibs.hwe.FIS, n.pca = 1, n.da = 4, parallel = require("parallel"))


labs <- c("Haida Gwaii","SouthWest BC")
cols=c("#009E73","#E69F00" )



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_2_clusters_419NDIV_1367SNPSSNPS_optimal_1PCaxis.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc2, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_2_clusters_419NDIV_1367SNPSSNPS_200PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc,col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()



## #DAPC with 80% variation explained####

highdapc <- dapc(gl.toad.nosibs.hwe.FIS, parallel = require("parallel"), n.da = 10, pca.select = "percVar", perc.pca=80)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_2_clusters_419NDIV_1367SNPSSNPS_80pctPCs.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(highdapc,col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()



########## membership prob

dapc.results <- as.data.frame(pnw.dapc2$posterior)
dapc.results$oldpop <- pop(gl.toad.nosibs.hwe.FIS)
dapc.results$indNames <- rownames(dapc.results)


dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample",
                            "Assigned_Pop","Posterior_membership_probability")



###### no sep areas 
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_2_clusters_419NDIV_1367SNPS_optimalPC_1PCaxis_nosepareas.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_blank(),axis.title.y = element_text( size = 10),axis.ticks.x = element_blank(),panel.background = element_blank())
t <- t + scale_fill_manual(values=c("#009E73","#E69F00" ),
                           name="Assigned Pop K=2 cluster", breaks = c("Haida_Gwaii","Southwest_BC"), labels=c("Haida_Gwaii","Southwest_BC"))
t
dev.off()


# with sample id
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_2_clusters_419NDIV_1367SNPS_optimalPC_1PCaxis_nosepareas_sampleid.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_text(angle = 90, size = 3),axis.title.y = element_text( size = 8))
t <- t + scale_fill_manual(values=c("#009E73","#E69F00" ),
                           name="Assigned Pop K=2 cluster", breaks = c("Haida_Gwaii","Southwest_BC"), labels=c("Haida_Gwaii","Southwest_BC"))
t
dev.off()



