**Exercise 3: Using other alignment methods**

BLAST is a very useful tool for sequence alignment, but other programs may be better suited to certain applications.

For example, say we have the assembled genome sequences of multiple varieties of a certain species but only one variety has been annotated with gene models. 
However, we would like to know where the genes are located in the other genomes. How should we go about finding the positions of these gene models in the genome assemblies of the other varieties?

In this exercise, we will take the genome assemblies from multiple different accessions of Arabidopsis (seattle, chi2, C24 and Bur) and try using BLAST and gmap to identify to locations of the Col-0 gene models in the other assemblies.

1. Download/copy the genome assemblies for the different Arabidopsis accessions
	- you can copy the files from here /home/lgardiner/Jemima/Intro_to_NGS_exercises/Exercise_3/arabidopsis_accessions/
  
  *Solution: e.g.*
  
 ```bash
 cp /home/lgardiner/Jemima/Intro_to_NGS_exercises/Exercise_3/arabidopsis_accessions/Bur.fa /my/home/dir/
 ```
 2. Produce some summary statistics about the genome sequences. 
		e.g. 	How many chromosomes?
				Are the genomes the same size?
				Are all the chromosomes the same length in the different accessions?
        
  *Solution: you can use samtools to generate a fasta index that will have the lengths of all the chromosomes*
  
  ```bash
  samtools faidx ref_genome.fa
  ```

3. Download the Col-0 (Reference) Arabidopsis transcriptome fasta file.

*Solution: available from ensembl plants (http://plants.ensembl.org/) and other websites*

4. How many transcripts in this file?

*Solution: you could calculate this on the command line e.g.*

```bash
grep ">" my_transcript_reference.fa | wc -l
```

5. Choose one of the Arabidopsis accessions to test out the different aligners

6. Use BLAST to find the locations of the Col-0 transcripts in the genome of your chosen accession (this will take some time, move onto step 9 while this runs!)
		TIP: As we are BLASTing many sequences, output to a file and in tabular format
		TIP: Can you get BLAST to just give you the top hit?
    
 *Solution: use the command line to run blastn*
 
 *First, make a blast database:*
 
 ```bash
 makeblastdb -in my_genome_fasta.fa -dbtype nucl
 ```
 Then BLAST against the database, use `-outfmt 6` to get a tabular output, `-max_target_seqs 1` and `-max_hsps 1` should give the top hit
```bash
blastn -query my_transcript_reference.fa -db ref_genome.fa -outfmt 6 -max_target_seqs 1 -max_hsps 1 -out transcripts_top_hits_table.txt
```

7. What is the average sequence identity?

*Solution: calculate from the % sequence identity column in the output table of your blast*

8. On average, what percentage of the gene length is covered in these BLAST hits?

*Solution: try out the qcovs option in blast, or you could calculate yourself using the alignment length column in the blast table and the transcript lengths*

9. Use gmap to find the locations of the Col-0 transcripts in your chosen accession.
	HINT: You will need to use the gmap manual to find out how to use gmap.
	
*Solution: for gmap, you first need to build a gmap genome from your reference sequence and then run gmap. Manual here: https://github.com/juliangehring/GMAP-GSNAP*

```bash
#-D specifies the location of the gmap genome, -d the name of the gmap genome
gmap_build -D /desired/path/to/location/of/gmap/genome/ -d gmap_genome_name ref_genome.fa

#--min-identity=0.9 filters for minimum of 90 % identity, -f=2 asks for gff3 gene format
gmap -d gmap_genome_name my_transcript_reference.fa --min-identity=0.9 -f=2
```
10. Play around with the different output options and parameters in gmap.

11. Which aligner is better suited to finding the locations of the transcripts? Why?
*Solution: It depends exactly what you want! But gmap is well designed for finding cDNAs in genomes - you might have issues with blasts because each intron will produce a separate hit on the genomic dna*

Further exercises:
12. Use gmap to identify the locations of transcripts in all the different Arabidopsis accessions.
*Solution: use gmap as shown before, you will need to make a gmap genome for each accession*

13. Extract the gDNA of the transcripts from the Arabidopsis accessions.
*Solution: you can generate a bed file from the gff output of the accessions, and then use bedtools getfasta to extract fasta sequences. Remember - in a bed file, the coordinates can't be larger than the max size of each chromosome and the start (col2) must be smaller than the end (col3)*

```bash
cat transcripts.bed
chr1	284	2731	transcript_1
chr1	8572	9829	transcript_2

bedtools getfasta ```

14. Choose a random set of 20 transcripts and write a script to perform multiple alignments of each gene using the gDNA from the different Arabidopsis accessions.
