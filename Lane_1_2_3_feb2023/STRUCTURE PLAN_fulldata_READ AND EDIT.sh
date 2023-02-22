1. create structure file in R ... set pops as regions

2. open structure file in excel and ****delete pop column***** completely and delete the name above the indiv column ('ind').... IF OPEN IN TXT EDITOR SHOUDL ONLY HAVE ONE TAB AT THE START OF THE HEADER ROW

3. open edited struc fiel in txt editor and check if has ONE empty tab at the start of the header row and only one tab between indiv namne and first locus value... 

4. edit numbered pop id file to only have samples in struc file (using pop edit file) 

5. edit python script to have pop id file the file that i just edited and the struc fiel the file i just edited

6. put both pop id and str file in same folder with python script.

7. double click python file to run it...

NB file names CANNOT hvae dots in them

8. check if created new str file which has numbered pop id column and see if has any extra rows at the end of file - if so, delte rows, and change to unix line endings... ADD ANOTHER TAB AT THE START OF THE HEADER SO THAT THERE ARE TWO TABS BEFORE THE LOCUS NAMES START (ONE FOR SAMPLE ID AND ONE FOR POP ID)


9. edit other files 

	a) main params file, 
	b) extra params file, 
	c) file with all the structure commands for k1 to k12, 
	d) and then the job script to run it on the cluster - remmeber toi change number of nodes and parallel to number of lines in file c)

10. check if runs locally 

cmd



cd C:\Users\Roseanna\Downloads\structure_aug18th2021\164INDIV_1644SNPS

C:\Users\Roseanna\Downloads\structure_windows_console\console\structure.exe -K 2 -m mainparams_164_1644_POPID -e extraparams -i genind_structure_164INDIV_1644SNPS_POPIDlabels.str -o genind_structure_164INDIV_1644SNPS_POPIDlabels_k2r1 -D 37847122




### SENDS OUTPUT TO LOG FILE
C:\Users\Roseanna\Downloads\structure_windows_console\console\structure.exe -K 2 -m mainparams_164_1644_POPID -e extraparams -i genind_structure_164INDIV_1644SNPS_POPIDlabels.str -o genind_structure_164INDIV_1644SNPS_POPID_k2r1 -D 37847122 > genind_structure_164INDIV_1644SNPS_POPID_screenoutput_k2r1.log




11. upload files a-d and also the structure file to a folder on cluster

######## upload file _nodirforstrfile - i.e. has dir (file path) for all files except the structure file
#python and job script
rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_denovo_June2021/stacks_gapped0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_3124SNPS/structure/job_structure_python_nov18th_2021_164_1644_k1to13_5reps_nodirforstrfile /drives/f/GBS_data_03_02_21/HaidaGwaii_denovo_June2021/stacks_gapped0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_3124SNPS/structure/jobSUB_structure_python_nov18_2021parallel_164_1644_k1to13_5reps_nodirforstrfile.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/ --progress

#main param, extra param, str file
rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_denovo_June2021/stacks_gapped0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_3124SNPS/structure/mainparams_164_1644_POPID /drives/f/GBS_data_03_02_21/HaidaGwaii_denovo_June2021/stacks_gapped0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_3124SNPS/structure/genind_structure_164INDIV_1644SNPS_POPIDlabels.str /drives/f/GBS_data_03_02_21/HaidaGwaii_denovo_June2021/stacks_gapped0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_3124SNPS/structure/extraparams roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/ --progress





12. test run job

## make dir 
164INDIV_1644SNPS/structure/

# mkdir 
outputs_k1to13_5rep 
joblogs_k1to13

cd /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/

module load nixpkgs/16.09 intel/2018.3 structure/2.3.4






########## testing whether the directory is the problem - it is!!! #######


### works
structure -K 1 -m mainparams_164_1644_POPID -e extraparams -i genind_structure_164INDIV_1644SNPS_POPIDlabels.str -o genind_structure_164INDIV_1644SNPS_POPIDlabels_k2r1 -D 37847122


# works
structure -K 1 -m mainparams_164_1644_POPID -e extraparams -i genind_structure_164INDIV_1644SNPS_POPIDlabels.str -o /outputs_k1to13_5rep/genind_structure_164INDIV_1644SNPS_POPIDlabels_k2r1 -D 37847122


## works
structure -K 1 -m mainparams_164_1644_POPID -e extraparams -i genind_structure_164INDIV_1644SNPS_POPIDlabels.str -o /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/outputs_k1to13_5rep/genind_structure_164INDIV_1644SNPS_POPIDlabels_k2r1 -D 37847122



# works
structure -K 1 -m /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/mainparams_164_1644_POPID -e extraparams -i genind_structure_164INDIV_1644SNPS_POPIDlabels.str -o /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/outputs_k1to13_5rep/genind_structure_164INDIV_1644SNPS_POPIDlabels_k2r1 -D 37847122


# works
structure -K 1 -m /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/mainparams_164_1644_POPID -e /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/extraparams -i genind_structure_164INDIV_1644SNPS_POPIDlabels.str -o /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/outputs_k1to13_5rep/genind_structure_164INDIV_1644SNPS_POPIDlabels_k2r1 -D 37847122


### DOESN'T WORK
structure -K 1 -m /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/mainparams_164_1644_POPID -e /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/extraparams -i /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/genind_structure_164INDIV_1644SNPS_POPIDlabels.str -o /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/outputs_k1to13_5rep/genind_structure_164INDIV_1644SNPS_POPIDlabels_k2r1 -D 37847122



#end job
CRTL + C

cd joblogs_k1to13

# yes there is a file there and there wasn't an error! I can submit the job!



## download
rsync -zv roseanna@cedar.computecanada.ca:/scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/outputs_k1to13_5rep/* /drives/f/GBS_data_03_02_21/HaidaGwaii_denovo_June2021/stacks_gapped0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_3124SNPS/structure/outputs_k1to13_5rep/  --progress

## download 3 weeks 1 node verion
rsync -zv roseanna@cedar.computecanada.ca:/scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/outputs_k1to13_5rep_3/* /drives/f/GBS_data_03_02_21/HaidaGwaii_denovo_June2021/stacks_gapped0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_3124SNPS/structure/outputs_k1to13_5rep_3/  --progress



#### now compress outputs_k1to13_5rep folder

### upload compressed folder outputs_k1to13_5rep to structure harvester
http://taylor0.biology.ucla.edu/structureHarvester/




## upload compressed folder outputs_k1to13_5rep to clumpak
http://clumpak.tau.ac.il/index.html

## save pdf job pipeline summary and zip file output



########## 3 weeks version ###
module load nixpkgs/16.09 intel/2018.3 structure/2.3.4


rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_denovo_June2021/stacks_gapped0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_3124SNPS/structure/take2_3weeks/job_structure_python_nov18th_2021_164_1644_k1to13_5reps_nodirforstrfile_3weeks /drives/f/GBS_data_03_02_21/HaidaGwaii_denovo_June2021/stacks_gapped0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_3124SNPS/structure/take2_3weeks/jobSUB_structure_python_nov18_2021parallel_164_1644_k1to13_5reps_nodirforstrfile_3weeks.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/ --progress


structure -K 1 -m /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/mainparams_164_1644_POPID -e /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/extraparams -i genind_structure_164INDIV_1644SNPS_POPIDlabels.str -o /scratch/roseanna/HaidaGwaii_denovo_June2021/stacks_gappedmin0.9.M3/HG_populations_R60pctoverall_p13_minmaf0.01_maxhet0.6_writesinglesnp/164INDIV_1644SNPS/structure/outputs_k1to13_5rep_2/genind_structure_164INDIV_1644SNPS_POPIDlabels_output_k1r1 -D 94109968 

















