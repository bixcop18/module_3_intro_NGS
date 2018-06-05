#!/bin/bash -e
#SBATCH -p batch
#SBATCH -n 4
#SBATCH -o ./trimmomatic_SRR3134441.%N.%j.out
#SBATCH -e ./trimmomatic_SRR3134441.%N.%j.err

cd /home/lgardiner/Jemima/Intro_to_NGS_exercises/Exercise_1/
#automatically loads 0.38
module load trimmomatic

#make directory for trimmed reads
#careful! If this directory already exists then this will error and kill the job...if this is the case you can remove this mkdir command
mkdir ./reads/SRR3134441/trimmed/

#PE for paired end
trimmomatic PE \
./reads/SRR3134441/SRR3134441_1.fastq.gz \
./reads/SRR3134441/SRR3134441_2.fastq.gz \
./reads/SRR3134441/trimmed/SRR3134441_1_trimmed.fastq.gz \
./reads/SRR3134441/trimmed/SRR3134441_1_unpaired.fastq.gz \
./reads/SRR3134441/trimmed/SRR3134441_2_trimmed.fastq.gz \
./reads/SRR3134441/trimmed/SRR3134441_2_unpaired.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

