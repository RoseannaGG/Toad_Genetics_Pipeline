
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