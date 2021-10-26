# Toad Genetics Pipeline Oct 2021

## Data info
We have GBS double-digested (Sbfl and Mspl), paired-end (2x101bp), Illumina Novaseq data with an individual level barcode on the forward read only. The data provided from the sequencing platform came back to us as demultiplexed plates with the Illumina adapter sequences and plate barcodes trimmed off (just leaving the individual barcode and cutsite on the forward read, and just the cutsite on the reverse read). We had nine plates - equating to 835 samples. DNA extractions done at UBC Hamelin lab, library prep done at Laval, sequencing done at Genome QC. Research funded by National Geographic, BC Ministry of Forests and the UBC Public Scholar Intiative. RGGs stipend provided by Govt. Canada Vanier Scholarship and BC Ministry of Forests. 

## 1. Create directories on server

## 2. Upload files to server  
   - raw reads forwad and reverse plates .gz files - do md5sum check before downloading from nanuq and also upon download to laptop and after uploaded to server 
   - stacks barcode files 
   - pop.map .tsv  
   - job scripts 

## 3. FastQC plates

Fine except low quality score for reverse read cut site

## 4. Trim reverse reads

- Take out the restriction enzyme cut site on reverse read

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 4 -trimlog /scratch/roseanna/Trimmed_reverseplates/D709_P9_R2_trimmed_log.txt /project/def-saitken/roseanna/rawreads_md5checked_feb272021/NS.1470.002.D709.Hamelin_202010_plate__5_R2.fastq.gz /scratch/roseanna/Trimmed_reverseplates/NS.1470.002.D709.Hamelin_202010_plate__5_R2.trimmed.fastq.gz HEADCROP:3



## 5. Demultiplex/clean plates - stacks/process_radtags

- discards reads with low quality (below 20) in sliding window of 25% of read length

process_radtags -1 /project/def-saitken/roseanna/rawreads_md5checked_feb272021/NS.1470.002.D707.Hamelin_202010_plate__3_R1.fastq.gz -2 /scratch/roseanna/Trimmed_reverseplates/NS.1470.002.D707.Hamelin_202010_plate__3_R2.trimmed.fastq.gz -b /scratch/roseanna/Demultiplexing_stacks/stacks_barcode_D707.txt -o /scratch/roseanna/Demultiplexing_stacks/D707_norenz2_trimR2 -w 0.25 -s 20 -y gzfastq --inline_null --renz_1 sbfI --quality --rescue --barcode_dist_1 1 -D &> process_radtags_standoutputerror_D707_norenz2_trimR2.oe


## 6. FastQC samples

## 7. Align the reads - stacks/ustacks/cstacks

## 8. Build loci - stacks/gstacks

## 9. Call snps - stacks/populations
- take first snp in locus (locus =  DNA between cut sites)
- - run with various dif filtering criteria

NB - run steps 7-9 on small test dataset first and play with values of M and n to see which to choose. As per Rochette and Catchan paper

### 9b.  nucleotide diverstiy 
  - run populations with some good filtering to get that has had some filtering



## 10. Filtering - O'Leary paper
  ### Filter biallelic, mac=3, min depth =5, mean min depth =5, GQ =20, rm indels
  
  VCF tools
  
  ### Fitler missing SNPs and missing individuals iteratively until got dataset with SNPS missing in less than 10% of indiv, and indivs missing less than 25% of SNPs.

  VCF tools
  
  ### Fitler samples over 0.6 max obs het - at regional level

  R - hierfstat

  ### Fitler out siblings
  
  R SNPrelate OR COLONY
  
  ### Fitler SNPs out of HWE

## 11. Plot PCA of filtered snps


## 12. FST comparison 

- between Haida Gwaii and Mainland & Vancouver Island
- within Haida Gwaii

## 13. Isolation-by-distance

## 14. Structure

## 15. FIS

## 16. Expected het





