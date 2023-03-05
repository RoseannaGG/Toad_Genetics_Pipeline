


################ no max het filter ##########

######## 0.125 2ND DEGREE OF MORE RELATED ################

ibddf <- read.table("F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/SambaR/run_noadditionalfilters/SambaR_output/Kinship/Close_kin.threshold0.2.txt",header = TRUE)

dim(ibddf)

summary(ibddf)

### subset to just to full sibs 
ibd<-ibddf[ibddf$kingrobust >= 0.25,]
dim(ibd) # 188763     22


#r1
ibd<-ibd[ibd$r1 >= 0.45,]
dim(ibd) # 188763     22


#ibd<-ibd[ibd$r1 <= 0.1,]
#dim(ibd) # 188763     22



### just within ponds
subset = ibd[substr(ibd$name1,1,9)==substr(ibd$name2,1,9),]
subset 
dim(subset) #7246   22

head(subset)

## save
write.table(subset, "F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/kingrobust_0.25_r1_0.45_relatedpairs_justwithinponds.txt", sep="\t",row.names = FALSE, col.names = TRUE,quote=FALSE )




### just isolate one row of names
namesib7<-ibd$name1
namesib7

dim(namesib7)


#deletes the repeats --- gives list of individuals to remove that I can feed back into vcf tools
namesib9<-unique(namesib7)
dim(namesib9) #260

#
write.table(namesib9, "F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/listindiv_kingrobust_0.25_r1_0.45.txt", row.names = FALSE, col.names = FALSE,quote=FALSE )


str(namesib9)

" chr [1:260] "R02-CW-KA-05" "

## get both pairs list
namesib8<-subset[c(3,4)]
namesib8

#write.table(namesib8, "F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/populationssnps_ibd_BOTHPAIRS_individualstoremovesibs2nddgree_hasduplicates_DP3_minGQ20_RMSIBSAFTER_05_07_21.txt", row.names = FALSE, col.names = FALSE,quote=FALSE )

#deletes the repeats -- same as above
namesib10<-unique(namesib8)
namesib10
write.table(namesib10, "F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/populationssnps_ibd_BOTHPAIRS_individualstoremovesibs2nddgree_ormore.txt", row.names = FALSE, col.names = FALSE,quote=FALSE )



################ king robust 0.25, and R1 = 1.9 ##########


ibddf <- read.table("F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/SambaR/run_noadditionalfilters/SambaR_output/Kinship/Close_kin.threshold0.2.txt",header = TRUE)

dim(ibddf)

summary(ibddf)

### subset to just to full sibs 
ibd<-ibddf[ibddf$kingrobust >= 0.25,]
dim(ibd) # 188763     22


#r1
ibd<-ibd[ibd$r1 >= 1.9,]
dim(ibd) # 188763     22


#ibd<-ibd[ibd$r1 <= 0.1,]
#dim(ibd) # 188763     22



### just within ponds
subset = ibd[substr(ibd$name1,1,9)==substr(ibd$name2,1,9),]
subset 
dim(subset) #7246   22

head(subset)

## save
write.table(subset, "F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/kingrobust_0.25_r1_1.9_relatedpairs_justwithinponds.txt", sep="\t",row.names = FALSE, col.names = TRUE,quote=FALSE )




### just isolate one row of names
namesib7<-ibd$name1
namesib7

dim(namesib7)


#deletes the repeats --- gives list of individuals to remove that I can feed back into vcf tools
namesib9<-unique(namesib7)
dim(namesib9) #260


write.table(namesib9, "F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/listindiv_kingrobust_0.25_r1_1.9.txt", row.names = FALSE, col.names = FALSE,quote=FALSE )


str(namesib9)

" chr [1:260] "R02-CW-KA-05" "

## get both pairs list
namesib8<-subset[c(3,4)]
namesib8

#write.table(namesib8, "F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/populationssnps_ibd_BOTHPAIRS_individualstoremovesibs2nddgree_hasduplicates_DP3_minGQ20_RMSIBSAFTER_05_07_21.txt", row.names = FALSE, col.names = FALSE,quote=FALSE )

#deletes the repeats -- same as above
namesib10<-unique(namesib8)
namesib10
write.table(namesib10, "F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/populationssnps_ibd_BOTHPAIRS_individualstoremove_kingrobust_0.25_r1_1.9.txt", row.names = FALSE, col.names = FALSE,quote=FALSE )




####### PLOT raw ######


ibddf <- read.table("F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/SambaR/run_noadditionalfilters/SambaR_output/Kinship/Close_kin.threshold0.txt",header = TRUE)

dim(ibddf)

summary(ibddf)

ibd<-ibddf

#ibd<-ibd[ibd$r1 <= 0.1,]
#dim(ibd) # 188763     22



### just within ponds
subset = ibd[substr(ibd$name1,1,9)==substr(ibd$name2,1,9),]
subset 
dim(subset) #7246   22

head(subset)

## save
#write.table(subset, "F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/kingrobust_0.25_r1_0.45_relatedpairs_justwithinponds.txt", sep="\t",row.names = FALSE, col.names = TRUE,quote=FALSE )

plot(subset$r1,subset$kincoef)

plot(subset$r1,subset$kingrobust)


plot(subset$r0,subset$kingrobust)


plot(subset$r1,subset$r0)



plot(subset$k1,subset$k2)


png("F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/r1vskingrobust_601INDIV_602SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
plot(subset$r1,subset$kingrobust)
dev.off()


png("F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/r0vskingrobust_601INDIV_602SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
plot(subset$r0,subset$kingrobust)
dev.off()

png("F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/r1vsr0_601INDIV_602SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
plot(subset$r1,subset$r0)
dev.off()

png("F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/k1vsk2_601INDIV_602SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
plot(subset$k1,subset$k2)
dev.off()



png("F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/r1vskincoef_601INDIV_602SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
plot(subset$r1,subset$kincoef)
dev.off()


