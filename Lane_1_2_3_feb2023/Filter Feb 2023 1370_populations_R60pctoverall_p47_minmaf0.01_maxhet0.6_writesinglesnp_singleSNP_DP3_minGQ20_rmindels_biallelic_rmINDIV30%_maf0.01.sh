# http://www.ddocent.com/filtering/
salloc --time=01:00:00 --cpus-per-task=2 --mem=1000M --ntasks=1 --account=def-saitken

module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.14
src=/scratch/roseanna/Denovo_May2021/stacks_gappedmin0.9.M3/
cd /scratch/roseanna/Denovo_May2021/stacks_gappedmin0.9.M3/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/



#############################################################
#############################  aug 8th 2021 #########################

########## remove R02-SS-EF samples #############


############ EXPORT unfiltered data #########################


#rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Denovo_May2021/stacks_gappedmin0.9.M3/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/gapped0.9.recode.vcf /drives/f/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/ --progress




### per site snp quality ###

vcftools --vcf $src/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/gapped0.9.recode.vcf --site-quality --out $src/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/gapped0.9_sitequality


#################################################################################
## 1. prelim filtering before removing sibs
module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.14
src=/scratch/roseanna/Denovo_May2021/stacks_gappedmin0.9.M3/

cd /scratch/roseanna/Denovo_May2021/stacks_gappedmin0.9.M3/1370_pops_R60pctoverall_p46_mmf001_mh06_wsp/


 #compute imissing
vcftools --vcf $src/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/populations.snps.vcf --missing-indv --out $src/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/


"
"

######### prelim fitlering
vcftools --vcf $src/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/populations.snps.vcf --remove-indels --minDP 3 --min-alleles 2 --max-alleles 2 --minGQ 20 --remove-filtered-all --recode --recode-INFO-all --out $src/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic

"


After filtering, kept 805 out of 805 Individuals
Outputting VCF file...
After filtering, kept 4697 out of a possible 4697 Sites
Run Time = 8.00 seconds




"



##############################
####### INDIV > 30% #########
################################


 #compute imissing
vcftools --vcf $src/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic.recode.vcf --missing-indv --out $src/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS

"







"

## 30%  ## remove 369
mawk '$5 > 0.3' 805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS.imiss | cut -f1 > 805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS30%missing.indv

cat 805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS30%missing.indv



### ACTUALLY REMOVE INDIVI
vcftools --vcf 805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic.recode.vcf --remove  805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS30%missing.indv --recode --recode-INFO-all --out 805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic_rmINDIV30%

"


Excluding individuals in 'exclude' list
After filtering, kept 371 out of 805 Individuals
Outputting VCF file...
After filtering, kept 4697 out of a possible 4697 Sites
Run Time = 4.00 seconds





"

### MAF
vcftools --vcf $src/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic_rmINDIV30%.recode.vcf --maf 0.01 --remove-filtered-all --recode --recode-INFO-all --out $src/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic_rmINDIV30%_maf0.01


"
After filtering, kept 371 out of 371 Individuals
Outputting VCF file...
After filtering, kept 2392 out of a possible 4697 Sites
Run Time = 2.00 seconds



"




########### EXPORT 602 INDIV 602 SNPS ###########

### 805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic_rmINDIV30%_maf0.01 ###


rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Denovo_May2021/stacks_gappedmin0.9.M3/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic_rmINDIV30%_maf0.01.recode.vcf roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Denovo_May2021/stacks_gappedmin0.9.M3/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS30%missing.indv /drives/f/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/371INDIV_2392SNPS/ --progress



rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Denovo_May2021/stacks_gappedmin0.9.M3/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/populations.log roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Denovo_May2021/stacks_gappedmin0.9.M3/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/populations.sumstats_summary.tsv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Denovo_May2021/stacks_gappedmin0.9.M3/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/populations.log.distribs roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Denovo_May2021/stacks_gappedmin0.9.M3/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/populations.snps.vcf /drives/f/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R60pctoverall_p27_minmaf0.01_maxhet0.6_writesinglesnp/ --progress


rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Denovo_May2021/stacks_gappedmin0.9.M3/catalog.calls roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Denovo_May2021/stacks_gappedmin0.9.M3/catalog.fa.gz /drives/f/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/ --progress





