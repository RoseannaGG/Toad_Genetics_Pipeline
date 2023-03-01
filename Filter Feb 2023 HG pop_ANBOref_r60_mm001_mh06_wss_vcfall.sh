


######## BELUGA #########




# http://www.ddocent.com/filtering/
salloc --time=03:00:00 --cpus-per-task=16 --mem=1000M --ntasks=1 --account=def-saitken



module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.14
src=/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/
cd /home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_mm001_mh06_wss_vcfall




#############################################################
#############################  feb 24th 2023#########################

########## remove R02-SS-EF samples #############


############ EXPORT unfiltered data #########################


#rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_mm001_mh06_wss_vcfall/gapped0.9.recode.vcf /drives/f/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/pop_ANBOref_r60_mm001_mh06_wss_vcfall/ --progress




### per site snp quality ###

vcftools --vcf $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/populations.snps.vcf --site-quality --out $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/sitequality


#################################################################################
## 1. prelim filtering before removing sibs
module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.14
src=/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/

cd /home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_mm001_mh06_wss_vcfall


 #compute imissing
vcftools --vcf $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/populations.snps.vcf --missing-indv --out $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/


"
"

######### prelim fitlering
vcftools --vcf $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/populations.snps.vcf --remove-indels --minDP 3 --min-alleles 2 --max-alleles 2 --minGQ 20 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic

"




After filtering, kept 385 out of 385 Individuals
Outputting VCF file...
After filtering, kept 150716 out of a possible 150716 Sites
Run Time = 85.00 seconds





"



##############################
####### INDIV > 30% #########
################################


 #compute imissing
vcftools --vcf $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic.recode.vcf --missing-indv --out $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS

"
"

mawk '!/IN/' pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS.imiss | cut -f5 > pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS_totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

                                       Histogram of % missing data per individual                      
                                                                                                       
  16 +-+-------------+---------------+---------------+--------------+---------------+---------------+-------------+-+
     +               +               +               +              +           **  +               +               +
     |                                                                          **                          ******* |
  14 +-+                                                                      ****  ***                           +-+
     |                                                                        * **  * *                       **    |
     |                                                                        * **  * *                       **    |
  12 +-+                                                                      * **  * *                   *** **  +-+
     |                                                                        * ***** *  *********    ** ** ****    |
     |                                                                        * ** ** *  * ** ** *    ** ** * **    |
  10 +-+                                                                 **** * ** ** *  * ** ** *    ** ** * **  +-+
     |                                                   *** ****     **** ** * ** ** *  * ** ** *    ** ** * **    |
   8 +-+                                        **       * * ** *   *** ** ** * ** ** **** ** ** ********** * **  +-+
     |                                          **       * * ** *   * * ** ** * ** ** ** * ** ** * ** ** ** * **    |
     |                                          **       * * ** **  * * ** ** * ** ** ** * ** ** * ** ** ** * **    |
   6 +-+                                 ***    ****    ** * ** **  * * ** **** ** ** ** * ** ** * ** ** ** * ** **-+
     |                                   * *    ** *    ** * ** **  * * ** ** * ** ** ** * ** ** * ** ** ** * ** ** |
     |                                   * **** ** *  **** **** ** ** * ** ** * ** ** ** * ** ** * ** ** ** * ** ** |
   4 +-+                                 * * ** ** *  * ** * ** ** ** * ** ** * ** ** ** * ** ** * ** ** ** * *****-+
     |                                   * * ** ** *  * ** * ** ** ** * ** ** * ** ** ** * ** ** * ** ** ** * ** ** |
     |                    **     ********************************************************************************** |
   2 +-+             **** **  ****  *  * * * ** ** **** ** * ** ** ** * ** ** * ** ** ** * ** ** * ** ** ** * ** **-+
     |               *  * **  *  *  *  * * * ** ** * ** ** * ** ** ** * ** ** * ** ** ** * ** ** * ** ** ** * ** ** |
     +   *************  *******  *  **** * * ** ** * ** ** * ** ***** * ** ** * ** ** ** * ** ** * ** ** ** * ** ** +
   0 +-+-**********************************************************************************************************-+
    0.3             0.4             0.5             0.6            0.7             0.8             0.9              1
                                                    % of missing data                                  
                                                                                                       



# 80%
# remove snps with over 30% missing data across all individuals
#exclude snps that are not present in at least 20% of the samples
vcftools --vcf $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic.recode.vcf --max-missing 0.8 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP0.8



"After filtering, kept 385 out of 385 Individuals
Outputting VCF file...
After filtering, kept 81 out of a possible 150716 Sites
Run Time = 11.00 seconds
"


# 70%
# remove snps with over 30% missing data across all individuals
#exclude snps that are not present in at least 30% of the samples
vcftools --vcf $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic.recode.vcf --max-missing 0.7 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP70%

"After filtering, kept 385 out of 385 Individuals
Outputting VCF file...
After filtering, kept 449 out of a possible 150716 Sites
Run Time = 15.00 seconds
"



# 60%
# remove snps with over 40% missing data across all individuals
#exclude snps that are not present in at least 40% of the samples
vcftools --vcf $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic.recode.vcf --max-missing 0.6 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%


"
After filtering, kept 385 out of 385 Individuals
Outputting VCF file...
After filtering, kept 3323 out of a possible 150716 Sites
Run Time = 10.00 seconds
"



 #compute imissing
vcftools --vcf $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%.recode.vcf --missing-indv --out $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS

"
"

mawk '!/IN/' pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS.imiss | cut -f5 > pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS_totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF




                                       Histogram of % missing data per individual                      
                                                                                                       
  25 +-+-------+---------+---------+---------+---------+----------+---------+---------+---------+---------+-------+-+
     +         +         +         +         +         +          +         +         +         +         +         +
     |                                                                                                      ******* |
     |                                                                                                              |
     |    **                                                                                                        |
  20 +-+  **                                                                                                      +-+
     |    **                                                                                                        |
     |    **                                                                                                        |
     |    **                                                                                                        |
  15 +-+  **                                                                                                      +-+
     |    **                                                                                                        |
     |    **      **          **                                                                                    |
     |    **      **          **                                                                                    |
     |    **** ** **          **                                                                                    |
  10 +-+  **********          ****                                                                                +-+
     |   ***********     **   ******                                                                                |
     |   ***********     **   ******                                                                                |
     |   ************ ** ***  **********  **  **                                                                    |
   5 +-+********************* **********  **  **                 **                                               +-+
     |  ********************* **************  ********        ** ** **             **       **                      |
     |  ********************* ************************  ***  ************ ****     ******** ***                     |
     |  ********************************************************************** **************** ******              |
     +  ***************************************************************************************************         +
   0 +-+***************************************************************************************************-------+-+
     0        0.1       0.2       0.3       0.4       0.5        0.6       0.7       0.8       0.9        1        1.1
                                                    % of missing data                                  
                     








## 30%  ## remove 369
mawk '$5 > 0.3' pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS.imiss | cut -f1 > pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_30%missINDV.indv

cat pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_30%missINDV.indv



### ACTUALLY REMOVE INDIVI
vcftools --vcf pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%.recode.vcf --remove  pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_30%missINDV.indv --recode --recode-INFO-all --out pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%

"



Excluding individuals in 'exclude' list
After filtering, kept 204 out of 385 Individuals
Outputting VCF file...
After filtering, kept 3323 out of a possible 3323 Sites
Run Time = 1.00 seconds






"

### MAF
vcftools --vcf $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%.recode.vcf --maf 0.01 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%_maf001


"

After filtering, kept 204 out of 204 Individuals
Outputting VCF file...
After filtering, kept 2702 out of a possible 3323 Sites
Run Time = 1.00 seconds




"

vcftools --vcf $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%.recode.vcf --maf 0.05 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%_maf005


"After filtering, kept 204 out of 204 Individuals
Outputting VCF file...
After filtering, kept 715 out of a possible 3323 Sites
Run Time = 1.00 seconds
"

########### EXPORT 602 INDIV 602 SNPS ###########

### pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%_maf001 ###


# maf 0.01 second time
rsync -zv roseanna@beluga.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%_maf001.recode.vcf roseanna@beluga.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_mm001_mh06_wss_vcfall/pop_ANBOref_r60_mm001_mh06_wss_vcfall_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_30%missINDV.indv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_mm001_mh06_wss_vcfall/204INDIV_2702SNPS/ --progress



rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_mm001_mh06_wss_vcfall/populations.log roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_mm001_mh06_wss_vcfall/populations.sumstats_summary.tsv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_mm001_mh06_wss_vcfall/populations.log.distribs roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_mm001_mh06_wss_vcfall/populations.snps.vcf /drives/f/GBS_data_03_02_21/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_mm001_mh06_wss_vcfall/ --progress


rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/catalog.calls roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/catalog.fa.gz /drives/f/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/ --progress





