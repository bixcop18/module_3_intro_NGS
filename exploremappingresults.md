

After mapping the fastq file to the reference genome you will end up with a SAM or BAM alignment file. 

* **SAM format:** 
SAM stands for Sequence Alignment/Map format.
A single SAM file can store mapped, unmapped, and even QC-failed reads from a sequencing run, and indexed to allow rapid access: this means that the raw sequencing data can be fully recapitulated from the SAM/BAM file.

* **BAM format:** SAM is rarely helpful and really takes up too much space which is why we use only the BAM in principle. A BAM file (.bam) is the binary version of a SAM file (saving storage and faster manipulation)

* SAM tools to explore SAM and BAM files

You can use **samtools**: a free software package for manipulating SAM/BAM files to manipulate your SAM/BAM files and extract different kind of information from them:
Samtools provide utilities for:

|Viewing and formatting|
|----------------------|
|Extracting statistics |
|Indexing              |
|Manipulating SAM/BAM files|
|Editing               |

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
samtools view -b test.sam > test.bam
```
or 
```
samtools view -bT test.sam > test.bam
#if the header is absent from the SAM file
```
or 
```
samtools view -bS test.sam > test.bam
#if the header is header information is available
```

**from BAM to SAM**
```
samtools view test.bam > test.sam
```
Use options –h and –H to deal with the header



* **samtools sort** 
```
#sorting a bam file
samtools sort test.bam –o test.bam
```
```
#converting SAM directly to a sorted BAM file
samtools view test.sam |samtools sort –o test.bam
```
SAM/BAM files can be sorted in multiple ways, e.g. by location of alignment on the chromosome, by read name, etc. Note that different alignment tools will output differently sorted SAM/BAM, and you might need differently sorted alignment files as input for different downstream analysis tools.
* **mapping statistics** 

```
samtools flagstat file.bam
```
does a full pass through the input file to calculate and print statistics such as:

|% reads mapped|
|--------------|
|% unmapped reads|
|% reads properly paired|
|Other information      |


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

**Extracting SAM entries mapping to a specific region**

```
#index the bam file first
samtools index file.bam 
```
```
samtools view file.bam chr1:200000-500000
#all reads mapping on chr1 as another bam 
samtools view –b file.bam chr1 > file_chr1.bam
```

* **Computing the depth**

Samtools allows computing the depth at each position using the depth option

```
samtools depth options file.bam
```

```
# –a allows to output all positions (including those with zero depth) 
samtools depth –a test.bam
```

```
#–q INT only count reads with base quality greater than INT
samtools depth –q int test.bam
```

* **Computing the coverage per region**

Samtools allows computing the read depth per genomic region, as specified in the supplied BED file using Samtools bedcov

```
samtools bedcov options region.bed file.bam
```
