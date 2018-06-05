#!/bin/bash

#specify paths to your working directory and your reads directory
WORKING_DIR='/home/lgardiner/Jemima/Intro_to_NGS_exercises/Exercise_1'
READ_DIR=$WORKING_DIR/reads

cd $WORKING_DIR
#create an overall output directory which will be unique with date and time
OUTPUT_DIR=$(date +out_%Y_%m_%d_%H_%M)

#make separate output directories within the overall output directories; -p indicates to create the whole path
mkdir -p ./$OUTPUT_DIR/run_logs
mkdir -p ./$OUTPUT_DIR/slurm_scripts
mkdir -p ./$OUTPUT_DIR/trimming_output

READS=($(ls $READ_DIR))

#the reads are in pairs in folders, iterate over each pair
for READ in ${READS[@]}
	do
	READS_PATH=$READ_DIR/$READ
	#specify the input reads and output names
	R1=$READ'_1.fastq.gz'
	R2=$READ'_2.fastq.gz'
	R1_trimmed=$READ'_1_trimmed.fastq.gz'
	R1_unpaired=$READ'_1_unpaired.fastq.gz'
	R2_trimmed=$READ'_2_trimmed.fastq.gz'
	R2_unpaired=$READ'_2_unpaired.fastq.gz'

#generate a unique slurm script for this pair	
echo \
"#!/bin/bash -e
#SBATCH -p batch
#SBATCH -n 4
#SBATCH -o ./$OUTPUT_DIR/run_logs/trimmomatic_$READ.%N.%j.out
#SBATCH -e ./$OUTPUT_DIR/run_logs/trimmomatic_$READ.%N.%j.err

#automatically loads 0.38
module load trimmomatic

#PE for paired end
trimmomatic PE $READS_PATH/$R1 $READS_PATH/$R2 \
./$OUTPUT_DIR/trimming_output/$R1_trimmed ./$OUTPUT_DIR/trimming_output/$R1_unpaired \
./$OUTPUT_DIR/trimming_output/$R2_trimmed ./$OUTPUT_DIR/trimming_output/$R2_unpaired  \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
" > ./$OUTPUT_DIR/slurm_scripts/$READ'_trimmomatic_slurm.sh'

sbatch ./$OUTPUT_DIR/slurm_scripts/$READ'_trimmomatic_slurm.sh'
done



