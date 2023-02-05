#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=1000
#SBATCH --job-name=jobsub_lane3_plate1_trim_demultiplex_norenz2_trimR2
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane3_plate1_trim_demultiplex_norenz2_trimR2_%j.out

# rsync -zv /drives/f/GBS_data_03_02_21/Trim_demultiplex_Oct2022/jobsub_B501_trim_demultiplex_norenz2_trimR2.sh /drives/f/GBS_data_03_02_21/Trim_demultiplex_Oct2022/stacks_barcode_B501.txt roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Demultiplexing_stacks/ --progress


module load StdEnv/2020
module load trimmomatic/0.39
module load stacks/2.60


echo “Starting run at: `date`”

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 4 -trimlog /scratch/roseanna/Trimmed_reverseplates/lane3_P1_R2_trimmed_log.txt /home/roseanna/scratch/rawreads_lane3_2023/xxxxxxxxx_R2.fastq.gz /scratch/roseanna/Trimmed_reverseplates/xxxxxxxxxxx.fastq.gz HEADCROP:3


echo “Starting demultiplex,clean at: `date`”

process_radtags -1 /home/roseanna/scratch/rawreads_lane3_2023/xxxxxx_R1.fastq.gz -2 /scratch/roseanna/Trimmed_reverseplates/xxxxxxx_R2.trimmed.fastq.gz -b /scratch/roseanna/Demultiplexing_stacks/stacks_barcode_lane3_plate1.txt -o /scratch/roseanna/Demultiplexing_stacks/lane3_plate1 -w 0.25 -s 20 -y gzfastq --inline_null --renz_1 pstI --quality --rescue --barcode_dist_1 1 -D &> process_radtags_standoutputerror_lane3_plate1.oe 

echo “Job finished with exit code $? at: `date`”