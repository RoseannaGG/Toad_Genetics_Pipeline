module load plink

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
  plink --file 766INDIV_3496SNPS --make-rel square --allow-extra-chr --out 766INDIV_3496SNPS_relatedness
done


"
[roseanna@cedar5 populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss]$ for j in "${ponds[@]}"
> do
>   plink --file 766INDIV_3496SNPS --make-rel square --allow-extra-chr --out 766INDIV_3496SNPS_relatedness
> done
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to 766INDIV_3496SNPS_relatedness.log.
Options in effect:
  --allow-extra-chr
  --file 766INDIV_3496SNPS
  --make-rel square
  --out 766INDIV_3496SNPS_relatedness

257860 MB RAM detected; reserving 128930 MB for main workspace.
Allocated 7259 MB successfully, after larger attempt(s) failed.
.ped scan complete (for binary autoconversion).
Performing single-pass .bed write (3496 variants, 766 people).
--file: 766INDIV_3496SNPS_relatedness-temporary.bed +
766INDIV_3496SNPS_relatedness-temporary.bim +
766INDIV_3496SNPS_relatedness-temporary.fam written.
3496 variants loaded from .bim file.
766 people (0 males, 0 females, 766 ambiguous) loaded from .fam.
Ambiguous sex IDs written to 766INDIV_3496SNPS_relatedness.nosex .
Using up to 63 threads (change this with --threads).
Before main variant filters, 766 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.869599.
3496 variants and 766 people pass filters and QC.
Note: No phenotypes present.
Relationship matrix calculation complete.
Relationship matrix written to 766INDIV_3496SNPS_relatedness.rel , and IDs
written to 766INDIV_3496SNPS_relatedness.rel.id .
"






rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS_relatedness.* /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/plink/ --progress






"
[roseanna@cedar5 populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss]$ plink --file 766INDIV_3496SNPS --make-rel square --allow-extra-chr --out 766INDIV_3496SNPS_relatedness
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to 766INDIV_3496SNPS_relatedness.log.
Options in effect:
  --allow-extra-chr
  --file 766INDIV_3496SNPS
  --make-rel square
  --out 766INDIV_3496SNPS_relatedness

257860 MB RAM detected; reserving 128930 MB for main workspace.
Allocated 7259 MB successfully, after larger attempt(s) failed.
.ped scan complete (for binary autoconversion).
Performing single-pass .bed write (3496 variants, 766 people).
--file: 766INDIV_3496SNPS_relatedness-temporary.bed +
766INDIV_3496SNPS_relatedness-temporary.bim +
766INDIV_3496SNPS_relatedness-temporary.fam written.
3496 variants loaded from .bim file.
766 people (0 males, 0 females, 766 ambiguous) loaded from .fam.
Ambiguous sex IDs written to 766INDIV_3496SNPS_relatedness.nosex .
Using up to 63 threads (change this with --threads).
Before main variant filters, 766 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.869599.
3496 variants and 766 people pass filters and QC.
Note: No phenotypes present.
Relationship matrix calculation complete.
Relationship matrix written to 766INDIV_3496SNPS_relatedness.rel , and IDs
written to 766INDIV_3496SNPS_relatedness.rel.id .
"




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
