Before getting started with the downstream analysis, it's always good to do some quality checks on your raw sequences using [FastQC tool](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to assess the quality of raw sequence data, the fastq files.
We will be using the  Seattle arabidopsis samples for the session. The raw sequences (fastq) files are in the Seattle_reads folder (under arabidopsis). 

* Create a new folder (arabidopsis) and then a subfolder Seattle_reads. Download fastq files there (links available [here](https://github.com/bixcop18/module_3_intro_NGS/blob/master/arabidopsis/Files.md))
* If you would like to run FastQC as a GUI just like in Windows type:

```
fastqc 
```
and you will be able to upload files using the GUI.

* When you have many fastq files to check, it's better to running FastQC from command line so you can loop over your files and process the reports automatically. Running it from the command line also allows you to combine and run a set of commands to perform complete analysis.
* FASTQC can process different types of files fastq, SAM and BAM files, you can use the -f option to tell fastqc upfront which format to expect (fastq if you are using fastq format). By default FastQC will try to guess the file format from the name of the input file. Anything 
ending in .sam or .bam will be opened as a SAM/BAM file (more about these formats on Day 4) and everything else will be treated as fastq format.

* Navigate to Seattle_reads and run the command (all fastq files in the directory)

```
fastqc ./*.fastq
```
or (you can specify as many files as you want)

```
fastqc seattle0_r1_filter.fastq seattle0_r2_filter.fastq
```
An html report page will be automatically created for each fastq file (Seattle_reads) without launching a user interface. You can load up these html pages in your browser to assess your data through graphs and summary tables.

* You can also use loops

```
for file in ./*.fastq
do
fastqc -f fastq ${file}
done
```

* As for most of the commands, the -h option allows you to browse the help file
```
fastqc -h
```

* For more information about the different QC report sections check this [FastQC mannual](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf)

* After running the fastqc command here, 2 .html files contain the FASTQC reports and graphs and can be opened in a browser. A zip file containing individual graphs and additional data files will be also created. By default, these files will be generated in the same folder, you might want to create a dedicated Fastqcresults folder. In this case, you have to use the option -o and specify the results directory.

```
mkdir Fastqcresults
fastqc -f fastq -o Fastqcresults seattle0_r1_filter.fastq seattle0_r2_filter.fastq
```
**Trimming and filtering**

After checking your fastq files quality, you have an idea about the quality of your reads, some of your reads might not be of a very good quality or the quality might drop at some positions (near the begining or end of reads) across all reads and this requires to clean up your library to minimize biaises in your analysis by  filtering poor quality reads and/or trim poor quality bases from our samples. There are many trimming tools that will allow to do so. 

* Trimmomatic (a program written in Java) performs a variety of trimming tasks for Illumina paired end and single end data. <br/>
A basic command to run trimmomatic, PE or SE specifies if your reads are Paired end or Single End:

```
java -jar trimmomatic-XX.jar PE 
bash java -jar Trimmomatic-XX.jar inputfile outputfile OPTION:VALUE...
```

You can use different options and parameters. The selection of trimming steps and their associated parameters are supplied on the command line.
Examples of parameters:

  **ILLUMINACLIP**: Cut adapter and other illumina-specific sequences from the read.<br/>
  **SLIDINGWINDOW**: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.<br/>
  **LEADING**: Cut bases off the start of a read, if below a threshold quality<br/>
  **TRAILING**: Cut bases off the end of a read, if below a threshold quality<br/>
  **CROP**: Cut the read to a specified length<br/>
  **HEADCROP**: Cut the specified number of bases from the start of the read<br/>
  **MINLEN**: Drop the read if it is below a specified length <br/>
  
Command example
```
java -jar trimmomatic-XX.jar PE seattle0_r1_filter.fastq seattle0_r2_filter.fastq seattle0_r1_filter_trimmed.fastq seattle0_r1_filter_trimmed_unpaired.fastq seattle0_r2_filter_trimmed.fastq seattle0_r2_filter_trimmed_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
or 
```
module load trimmomatic
trimmomatic PE seattle0_r1_filter.fastq seattle0_r2_filter.fastq seattle0_r1_filter_trimmed.fastq seattle0_r1_filter_trimmed_unpaired.fastq seattle0_r2_filter_trimmed.fastq seattle0_r2_filter_trimmed_unpaired.fastq  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```




You can specify the order in which you would like to perform the different cleaning steps. It is recommended in most cases that adapter clipping, if required, is done as early as possible, since correctly identifying adapters using partial matches is more difficult.


The output (trimmed reads) are in gzipped fastq format.

More information about Trimmomatic [here](http://www.usadellab.org/cms/?page=trimmomatic)


* Cutadapt allows to search and remove these adapter sequences from your reads

NGS data sets contain a high number of contaminating adapter sequences. You can also use cutadapt to search and remove these adapters from your reads. All reads that were present in the input file will also be present in the output file, some of them trimmed, some of them not. <br/>
Command to run cutadapt:
```
cutadapt -a ADAPTER -o output.fastq input.fastq
```
Many different options are available with cutadapt, you can find more information [here](https://cutadapt.readthedocs.io/en/stable/)

* After performing filtering and trimming, your samples should perform better on the quality tests run by FastQC. You can re-run FastQC on your trimmed fastq files to check the base sequence quality, it should be higher after trimming.
