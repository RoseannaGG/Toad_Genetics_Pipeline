#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=168:00:00
#SBATCH --mem=100000
#SBATCH --job-name=jobsub_lane1_lane2_lane3_gstacks_populations
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane1_lane2_lane3_gstacks_populations_%j.out
#
# Run gstacks: build a paired-end contig from the metapopulation data (if paired-reads provided),
# align reads per sample, call variant sites in the population, genotypes in each individual.


gstacks -I $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignmentsmapq20/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20/ --min-mapq 20 -t $cpu



echo “Starting populations_ANBOref_maxhet0.6_writesinglesnp at: `date`”

# mkdir populations_ANBOref_maxhet0.6_writesinglesnp populations_ANBOref_writesinglesnp populations_ANBOref_r10_R70pctoverall_p44_minmaf0.01_maxhet0.6_writesinglesnp populations_ANBOref_r10_R70pctoverall_p42_minmaf0.01_maxhet0.6_writesinglesnp

#
# Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics
# export several output files.


populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20/populations_ANBOref_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --plink --structure --genepop --treemix --verbose -t $cpu



echo “Starting populations_ANBOref_writesinglesnp no filter except write single snp at: `date`”

#
# Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics
# export several output files.


populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20/populations_ANBOref_writesinglesnp/ --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --plink --structure --genepop --treemix --verbose -t $cpu



echo “Job finished with exit code $? at: `date`”



echo “Starting populations_ANBOref_r10_R70pctoverall_p44_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20/populations_ANBOref_r10_R70pctoverall_p44_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-overall 0.7 --min-samples-per-pop 0.1 --min-populations 44 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --plink --structure --genepop --treemix --verbose -t $cpu


echo “Job finished with exit code $? at: `date`”

echo “Starting populations_ANBOref_r10_R70pctoverall_p42_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20/populations_ANBOref_r10_R70pctoverall_p42_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-overall 0.7 --min-samples-per-pop 0.1 --min-populations 42 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --plink --structure --genepop --treemix --verbose -t $cpu


echo “Job finished with exit code $? at: `date`”


## download populations outputs to: F:\GBS_data_03_02_21\Lane_1_2_3_feb2023\HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW
