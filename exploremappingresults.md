

After mapping the fastq file to the reference genome you will end up with a SAM or BAM alignment file. SAM and BAM files contain a lot of different information which can be difficult to understand at first, but knowing how to extract the data you want from these files is essential for a bioinformatician. 

* **SAM format:** 
SAM stands for Sequence Alignment/Map format.
A single SAM file can store mapped, unmapped, and even QC-failed reads from a sequencing run, and indexed to allow rapid access: this means that the raw sequencing data can be fully recapitulated from the SAM/BAM file.

* **BAM format:** SAM is rarely helpful and really takes up too much space which is why we use only the BAM in principle. A BAM file (.bam) is the binary version of a SAM file (saving storage and faster manipulation)

* SAM tools to explore SAM and BAM files

You can use **samtools**: a free software package for manipulating SAM/BAM files to manipulate your SAM/BAM files and extract different kind of information from them.
Samtools provide utilities for: <br/>

* Viewing and formatting
* Extracting statistics
* Indexing 
* Manipulating SAM/BAM files
* Editing 

```
Usage:  samtools <command> [options]
```
samtools offers many options, you can check them [here](http://www.htslib.org/doc/samtools.html)

One of the most used tools  since BAM files are often the input files needed for many different analysis programs.
```
samtools view
```

**from SAM to BAM**
```
samtools view -b file.sam > file.bam
```
or 
```
samtools view -bT file.sam > file.bam
#if the header is absent from the SAM file
```
or 
```
samtools view -bS file.sam > file.bam
#if the header is header information is available
```

**from BAM to SAM**
```
samtools view file.bam > file.sam
```
Use options –h and –H to deal with the header, this indicates that the header should be put into the SAM file (contains information about the number of chromosomes etc.)



* **samtools sort** 
```
#sorting a bam file
samtools sort file.bam –o file_sorted.bam
```
```
#converting SAM directly to a sorted BAM file
samtools view file.sam |samtools sort –o file_sorted.bam
```
SAM/BAM files can be sorted in multiple ways, e.g. by location of alignment on the chromosome, by read name, etc. Note that different alignment tools will output differently sorted SAM/BAM, and you might need differently sorted alignment files as input for different downstream analysis tools.<br/>
Look at the sizes of the unsorted and sorted BAM files, notice the sorted bam is around ¼ smaller than the unsorted: don’t worry! This is normal, sorting the BAM puts the reads in an order that allows for more efficient compression leading to a smaller output file. 

* **Mapping statistics** 

```
samtools idxstats file.bam
```
The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads.<br/>

Easy ways to extract 

```
#number of reads
samtools idxstats file.bam | awk '{s+=$3+$4} END {print s}'
```
```
#number of mapped reads
samtools idxstats file.bam | awk '{s+=$3} END {print s}'
```
```
samtools flagstat file.bam
```
Does a full pass through the input file to calculate and print statistics such as: %reads mapped, % unmapped reads, % reads properly paired and Other information. <br/>


```
#Count the total number of alignments
samtools view file.bam | wc -l
```



Many tools require a BAM Index file to more efficiently access reads in a BAM file.  
To create a BAM index, you  must first sort the BAM file to create a sorted.bam and then run samtools index with the sorted.bam as input
This will create a file named sorted.bam.bai which contains the index.

```
samtools view  file.sam >file.bam
samtools sort file.bam -o file_sorted.bam
samtools index file_sorted.bam file_sorted.bai
```

* **Filtering out unmapped reads from BAM files**

```
samtools view -h -F 4  file.bam > file_only_mapped.sam
# output back to BAM
samtools view -h -F 4 –b file.bam > file_only_mapped.bam
```
This leaves you with a BAM file that contains ONLY sequences reads that DID map to the reference genome.


* **Extracting unmapped reads from BAM files**

```
samtools view -h -f 4 -b  file.bam > file_only_unmapped.sam
# output back to BAM
samtools view -h -f 4 file.bam > file_only_unmapped.bam
```
This leaves you with a BAM file that contains ONLY sequences reads that DID NOT map to the reference genome. <br/>

Note: using –F tells SAMtools view to “skip” sequence alignments with a specific flag (in this case it skips all of the unmapped sequences and only outputs the mapped sequences). Using –f tells SAMtools to “skip” sequence alignments that DON’T match a specific flag 

* **Run flagstat again on these files containing only unmapped and only mapped reads and see what changed**

**Extracting SAM entries mapping to a specific region**

```
#index the bam file first
samtools index file.bam 
```
```
samtools view file.bam Chr1:200000-500000 >file_chr1_200k-500k.sam
#all reads mapping on chr1 as another bam 
samtools view –b file.bam chr1 > file_chr1.bam
```
You can visualise the content of file_chr1_200k-500k.sam using less, more, etc. or any text editor. This is very useful when you have a small region of interest allowing you to interrogate your region with much smaller files (reducing analysis time!). This analysis requires a sorted BAM file and an index file for the BAM file.

* **Subsample alignments from a BAM file**

 ```
 samtools view -s 0.5 -b ./file_sorted.bam  > random_half_of_file.bam
 ```

This selects 50% of the alignments from the BAM file at random. This is useful when you want to know how your alignment might have looked if you did 50% less sequencing, very useful if you are planning future experiments.


* **Computing the depth**

Samtools allows computing the depth at each position using the depth option

```
samtools depth options file.bam
```

```
# –a allows to output all positions (including those with zero depth) 
samtools depth –a file.bam
```

```
#–q INT only count reads with base quality greater than INT
samtools depth –q int file.bam
```

* **Computing the coverage per region**

Samtools allows computing the read depth per genomic region, as specified in the supplied BED file using Samtools bedcov

```
samtools bedcov options region.bed file.bam
```

* **SAMtools fillmd – visualising alignment matches/mismatches**
```
samtools view -b file_sorted.bam | samtools fillmd -e - ./yourreference.fasta > file_fillmd.sam
```

This changes the sequence section of the SAM file to a string of “=” when the alignment was a match while leaving the base code e.g. “A” if there was not a match at that position. Use the “less” command to explore the file_fillmd.sam file. 

* **Extract fastq sequence reads from a BAM file**

```
samtools fastq file.bam > file_extract.fastq
```

You can also extract the sequences into R1 and R2 files… check the options for this tool by entering “samtools fastq” into the command line

* **Extract fasta sequence reads from a BAM file**

```
samtools fasta file_sorted.bam > file_extract_test.fasta
```

