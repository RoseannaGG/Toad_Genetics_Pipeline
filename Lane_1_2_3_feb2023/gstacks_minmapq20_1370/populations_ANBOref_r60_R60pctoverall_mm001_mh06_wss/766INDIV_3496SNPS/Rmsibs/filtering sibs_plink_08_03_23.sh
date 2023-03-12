
rsync -zv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3356SNPS_editlocinamesinR_obshet06Perregion_R.vcf --progress roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/plink/

cd /home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/plink/


module load StdEnv/2020 plink/1.9b_6.21-x86_64

declare -a ponds=("R01-CR-CE
R01-CR-RA
R01-NI-CX
R01-NI-OC
R01-SI-CA
R01-SI-DB
R01-SI-FR
R01-SI-GL
R01-SI-HR
R01-SI-KE
R01-SI-LC
R01-SI-MO
R01-SI-SJ
R02-CW-AG
R02-CW-AL
R02-CW-BE
R02-CW-DE
R02-CW-JO
R02-CW-KA
R02-CW-LA
R02-CW-MI
R02-CW-RY
R02-CW-SI
R02-SC-IN
R02-SS-CR
R02-SS-FA
R02-SS-LC
R02-SS-LT
R06-GH-DL
R06-GH-GW
R06-GH-LL
R06-GH-PQ
R06-GH-PT
R06-GI-CK
R06-GI-CN
R06-GI-EV
R06-GI-GL
R06-GI-LV
R06-GI-MY
R06-GI-RR
R06-MI-MM
R06-SK-FE
R06-SK-ML
R06-SK-ON
")

for j in "${ponds[@]}"
do
  plink --vcf 766INDIV_3356SNPS_editlocinamesinR_obshet06Perregion_R.vcf --make-rel square --allow-extra-chr --out 766INDIV_3356SNPS_editlocinamesinR_obshet06Perregion_R_relatedness
done





rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/plink/766INDIV_3356SNPS_editlocinamesinR_obshet06Perregion_R_relatedness.* /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/plink/766INDIV_3356SNPS/ --progress








conda activate VCFTOOLS

for c in "${ponds [@]}"
do
  vcftools --vcf $DATDIR/"$c".vcf --missing-ind $OUTDIR/"$c"
done


#AFTER R script
ROUT=R OUTPUT FOLDER
FINALOUT=FINAL OUTPUT DIRECTORY

for y in "${ponds[@]}"
do
  vcftools --remove-indv $ROUT/"$y"_NoSibs.txt --vcf $DATDIR/"$y".vcf --out $FINALOUT/"$y"_NoSibs.vcf
