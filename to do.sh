## to do

# tranfer raw gbs files from hard drive to compute canada server using globus



# generate md5sum for all files:
a) on server

md5sum NS.2053.003.B716---D502.Hamelin__20230109-Plate-2_R2.fastq.gz

b) on hard drive (using windows commander)

CertUtil -hashfile "I:\Lane_3_GBS_data_13_02_23\NS.2053.003.B716---D502.Hamelin__20230109-Plate-2_R2.fastq.gz" MD5

CertUtil -hashfile "F:\GBS_data_03_02_21\Lane_3_GBS_data_13_02_23\NS.2053.003.B716---D502.Hamelin__20230109-Plate-2_R2.fastq.gz" MD5


# check md5sum of files on server

# get md5sum from Daryl

# run md5sum on anaxy ref genome

# upload to server

# run md5sum on genome on server



####### issues to resolve

# check lane 1 working with Pstl  (YES)

# check lane 1 working with anaxyrus boreas ref (YES)

# check lane 1 working with PstI and anaxyrus boreas -(12 sample test D701 - CURRENTLY RUNNING)



# sort lane 2 process_radtags error (YES)

# trim and demultiplex one plate lane 2 (B503 - CURRENTLY RUNNING

# copy demultiplexed forward reads from lane 2 plate to Demultiplexing_stacks/lane1_lane2_lane3 -DONE USING GLOBUS SYNC SETTING AND VERIFY

# align that one plate from lane 2 (B503) - SEND TO MAIN FOLDER SCRATCH alignments - ANBO_refassembly_HGthesis_lane1_lane2_lane3
# had 50000 mem, it ran out of memory - so give 100000 next time aligning a plate

# test gstacks and pops with stacks 2.62 using some samples from B503 (8) and D701 (12)  (yes) 

 rerun B503 and D701 with pstI and stacks 2.62 (CHECK if used) 



 # trim and demultiplex and align ALL lane 1 and 2 plates (# when done)

  rsync -zv /drives/f/GBS_data_03_02_21/Lane2_GBSdata_2022/Trim_demultiplex_Oct2022/*.txt /drives/f/GBS_data_03_02_21/Lane2_GBSdata_2022/Trim_demultiplex_Oct2022/*.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Demultiplexing_stacks/  --progress

 rsync -zv /drives/f/GBS_data_03_02_21/Lane1_GBSdata_2021/Trim_demultiplex_Oct2021/*.txt /drives/f/GBS_data_03_02_21/Lane1_GBSdata_2021/Trim_demultiplex_Oct2021/*.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Demultiplexing_stacks/  --progress


  rsync -zv /home/roseanna/scratch/tests.ref.ANBO/jobsub_D701_trim_demultiplex_norenz2_trimR2_PstI_test12samples* /home/roseanna/scratch/Demultiplexing_stacks/  --progress

  rsync -zv /home/roseanna/scratch/tests.ref.ANBO/process_radtags_standoutputerror_D701_norenz2_trimR2_PstI.oe /home/roseanna/scratch/Demultiplexing_stacks/  --progress



 
## currenly running all these in blue below (9/2/23)
 # D701 # run but need to check if psti and 2.62 (yes psti but not 2.62)
 # D702
 # D703
 # D704
 # D705
 # D706
 # D707
 # D708
 # D709

 #B501
 #B502
 #B503 # run but need to check if psti and 2.62 (yes psti but not 2.62)
 #B504
 #D502
 #D503
 #D504 

# copy all FORWARD DEMULIPLEXED reads into one folder Demultiplexing_stacks/lane1_lane2_lane3 (ACTUALLY, THIS ISN'T NECCESSARY NOW)



## ALIGNMENT ###

# align whole D701 and B503 and send to specific folders (B503 DONE)

 rsync -zv /drives/f/GBS_data_03_02_21/Lane1_GBSdata_2021/jobsub_lane1_D701_align.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/  --progress

 re align B503 with updated samtools and more memory, to see if errors are fixed

transfer both D701 and B503 align to one folder in tests.ref.ANBO
  rsync -zv /home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignmentsmapq20_B503_moremem_updatedsamtools/*.bam /home/roseanna/scratch/tests.ref.ANBO/alignmentsmapq20_B503_D701_moremem_updatedsamtools/ --progress

 rsync -zv /home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignmentsmapq20_D701/*.bam /home/roseanna/scratch/tests.ref.ANBO/alignmentsmapq20_B503_D701_moremem_updatedsamtools/ --progress


# run gstacks jobsub_lane1_D701_12_lane2_B503_8_gstacks_pop_fullfolders.sh to see if can call only some files form align folder with all samples (YES)


 # align lane 1 and 2 to toad genome per plate (once errors sorted)


  rsync -zv /drives/f/GBS_data_03_02_21/Lane2_GBSdata_2022/*.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/  --progress

  rsync -zv /drives/f/GBS_data_03_02_21/Lane1_GBSdata_2021/*.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/  --progress

 - SEND TO MAIN FOLDER SCRATCH - ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_minmq20

 #D701 
 #D702
 #D703
 #D704
 #D705
 #D706
 #D707
 #D708
 #D709

 #B501
 #B502
# B503
 #B504
 #D502
 #D503
 #D504 



# check alignment

# b) check alignment using samtools

module load samtools/1.16.1 StdEnv/2020

x) samtools view -q 30 -c R02-SC-IN-01.bam # counts reads mapped with over 30 quality
y) samtools flagstat R02-SC-IN-01.bam -O tsv # gives summary of total numner of reads, % mapped etc. 

Divide (x) by (y) to get % mapped reads over 30 quality



# create barcode files for lane 3 plates (6 files) - (YES)

# create pop map file for all samples (lane 1-3) (YES)


 when lane 3 samples arrive

 module load fastqc

fastqc NS.1760.001.B711---D503.Hamelin_202110_plate2_R1.fastq.gz

 fastqc plates using interactive node salloc

 rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/rawreads_lane3_2023/*.html /drives/f/GBS_data_03_02_21/Lane3_GBSdata_2023/FastQC/ --progress

 # create 6 trim + demuliplex scripts for lane 3 plates

 demultiplex lane 3 - one plate test - copy demultiplexed forward reads to Demultiplexing_stacks/lane1_lane2_lane3


 align lane 3 - one plate test - SEND TO MAIN FOLDER SCRATCH alignments - ANBO_refassembly_HGthesis_lane1_lane2_lane3


run gstacks and populations on lane 3 test

trim demultiplex all plates lane 3

 run gstacks and poulations on lane 1-3 combined using stacks v 2.62 (exports vcf file with all sites)

check missing data with different settings

 download vcf files and do analysis









