#!/bin/bash

#specify paths to your working directory and your reads directory
WORKING_DIR='/home/lgardiner/Jemima/Intro_to_NGS_exercises/Exercise_1'
#this is where we previously saved the trimmed reads
READ_DIR=$WORKING_DIR/out_2018_06_04_18_15/trimming_output/
REFERENCE=$WORKING_DIR/reference/Zea_mays.AGPv4.dna.toplevel.fa.gz

cd $WORKING_DIR
#create an overall output directory which will be unique with date and time
OUTPUT_DIR=$(date +out_%Y_%m_%d_%H_%M)

#make separate output directories within the overall output directories; -p indicates to create the whole path
mkdir -p ./$OUTPUT_DIR/run_logs
mkdir -p ./$OUTPUT_DIR/slurm_scripts
mkdir -p ./$OUTPUT_DIR/alignment_output

#iterate over just the trimmed R1 reads (so we don't dupicate alignments by iterating over each file in the three pairs)
READS=($(ls $READ_DIR | grep "_1_trimmed.fastq.gz"))

#the reads are in pairs in directories with the sample name, iterate over each pair (u)
for READ in ${READS[@]}
	do
	READ_NAME=${READ//_1_trimmed.fastq.gz}
	#specify the input reads and output names
	R1=$READ
	#replace the ending for read 1 with the ending for read 2
	R2=$READ_NAME'_2_trimmed.fastq.gz'

	
#generate a unique slurm script for this pair	
echo \
"#!/bin/bash -e
#SBATCH -p batch
#SBATCH -n 4
#SBATCH -o ./$OUTPUT_DIR/run_logs/$READ_NAME.bwa-mem.%N.%j.out
#SBATCH -e ./$OUTPUT_DIR/run_logs/$READ_NAME.bwa-mem.%N.%j.err


module load bwa

bwa mem -t 4 \
-R '@RG\tID:$READ_NAME\tPL:Illumina' \
-M $REFERENCE \
$READ_DIR/$R1 $READ_DIR/$R2 > $OUTPUT_DIR/alignment_output/$READ_NAME'_trimmed.sam'
" > ./$OUTPUT_DIR/slurm_scripts/$READ'_bwa-mem.sh'

chmod u+x ./$OUTPUT_DIR/slurm_scripts/*
chmod u+x ./reference/*

sbatch ./$OUTPUT_DIR/slurm_scripts/$READ'_bwa-mem.sh'
done



