analyzing your samples, it's always good to do some quality checks on your raw sequences using [FastQC tool](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to assess the quality of raw sequence data, the fastq files.
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

* Go under Seattle_reads and run the command (all fastq files in the directory)

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
* After running the fastqc command here, 2 .html files contain the FASTQC reports and graphs and can be opened in a browser. A zip file containing individual graphs and additional data files will be also created. By default, these files will be generated in the same folder, you might want to create a dedicated Fastqcresults folder. In this case, you have to use the option -o and specify the results directory.

```
mkdir Fastqcresults
fastqc -f fastq -o Fastqcresults seattle0_r1_filter.fastq seattle0_r2_filter.fastq
```


