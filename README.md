# whatGene
Tools to search for genes in special chromosomes like sexual or B chromosomes.

## Installation
- Copy script to your binaries folder.
- Dependencies:
  * RepeatMasker [http://www.repeatmasker.org/RMDownload.html](http://www.repeatmasker.org/RMDownload.html)
  * seqTK [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)
  * BLAT [http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/)
  * Trimmomatic [http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
  * R (ggplot2)
  * pysamstats [https://github.com/alimanfoo/pysamstats](https://github.com/alimanfoo/pysamstats)
  * Python (Matplotlib, Pylab and Biopython)
  * RAxML [https://cme.h-its.org/exelixis/web/software/raxml/index.html](https://cme.h-its.org/exelixis/web/software/raxml/index.html) 
  * MAFFT [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/)

## Protocol 0. De novo assembly, CDS prediction and clustering

### 0.1 Assembly (Trinity)

```
$ /usr/local/lib/trinityrnaseq-Trinity-v2.5.0/Trinity --seqType fq --max_memory 50G --left Paut-0B-Ind11_1.fq.gz,Paut-2B-Ind19_1.fq.gz --right Paut-0B-Ind11_2.fq.gz,Paut-2B-Ind19_2.fq.gz --CPU 12 --trimmomatic --full_cleanup
```

```
#  --normalize_max_read_cov <int>	defaults to 50
#  --normalize_by_read_set		run normalization separate for each pair of fastq files,
#					then one final normalization that combines the individual
#					normalized reads.
#					Consider using this if RAM limitations are a consideration.
```

output - assembly: Trinity.fasta

### 0.2 CDS extraction (Transdecoder)


```
$ TransDecoder.LongOrfs -t Trinity.fasta
```

output: emon_zb_trinity.fasta.transdecoder_dir/longest_orfs.cds

### 0.3 Clustering (CD-HIT)

Slow:

```
$ cd-hit-est -T 8 -r 1 -i longest_orfs.cds -o Trinity_cds_nr80.fasta -M 0 -aS 0.8 -c 0.8 -G 0 -g 1
```

Alternatively (fast):

```
$ cd-hit-est -r 1 -i longest_orfs.cds -T 12 -c 0.80 -o longest_orfs.cds.nr80 -M 0
```

output: longest_orfs.cds.nr80

## Protocol 1. Coverage graphics

### 1.1 Map reads against a transcriptomic reference

```
$ ssaha2_run_multi.py list.txt reference.fasta 12
```
or
```
$ ssaha2_run_multi_pe_se.py list.txt reference.fasta 12
```

### 1.2 Generate coverage file

```
$ bam_coverage_join.py FastaReference ListOfBams 999999]
```

### 1.3 Generate graphics

You need to prepare a samples files with libraries sorted like in the latter command. This “samples.txt” file should have this structure for each library:

Name    TypeOfLibrary        LibrarySIzeInNucleotides    GenomeSizeInNucleotides

For each reference sequence, the script will generate a pdf page. Looking at the type of library, the script will create a line chart for each different term before the “_” and for each different term after the “_” will include a line in the chart, calculating the mean+/-stdev if you have two or more libraries with the same term.

An example:

```
01-LMIG-PAT-LEG-M01_checked_1_mapped    gDNA_SLCa0B     20493125000     6300000000
04-LMIG-PDH-LEG-M12_checked_1_mapped    gDNA_SLCa0B     22651377750     6300000000
02-LMIG-PAT-LEG-M04_checked_1_mapped    gDNA_SLCaPB     20184562750     6300000000
03-LMIG-PDH-LEG-M04_checked_1_mapped    gDNA_SLCaPB     18438224500     6300000000
05-LMIG-PDH-LEG-M14_checked_1_mapped    gDNA_SLCaPB     20354837250     6300000000
06-LMIG-PDH-LEG-M15_checked_1_mapped    gDNA_SLCaPB     21084833750     6300000000
SRR764583_1_mapped      gDNA_NLW        277474348600    6300000000
01_LMIG_PAT_TES_M01_L12_1_mapped        RNAtestis_SLCa0B        31.497002       1
04_LMIG_PDH_TES_M12_L12_1_mapped        RNAtestis_SLCa0B        35.184289       1
02_LMIG_PAT_TES_M04_L12_1_mapped        RNAtestis_SLCaPB        35.837611       1
03_LMIG_PDH_TES_M04_L12_1_mapped        RNAtestis_SLCaPB        23.952951       1
05_LMIG_PDH_TES_M14_L12_1_mapped        RNAtestis_SLCaPB        48.226476       1
06_LMIG_PDH_TES_M15_L12_1_mapped        RNAtestis_SLCaPB        33.787235       1
07_LMIG_PAT_LEG_M01_L12_1_mapped        RNAleg_SLCa0B   32.427158       1
10_LMIG_PDH_LEG_M12_L12_1_mapped        RNAleg_SLCa0B   38.024871       1
08_LMIG_PAT_LEG_M04_L12_1_mapped        RNAleg_SLCaPB   33.235481       1
09_LMIG_PDH_LEG_M04_L12_1_mapped        RNAleg_SLCaPB   38.558194       1
11_LMIG_PDH_LEG_M14_L12_1_mapped        RNAleg_SLCaPB   36.287354       1
12_LMIG_PDH_LEG_M15_L12_1_mapped        RNAleg_SLCaPB   38.215901       1
```
Then, we can run:

```
$ coverage_graphics.py CoverageFile SamplesFile FastaFile [PDF/SVG/NOPLOT]
```

## Protocol 2. SNP calling

### 2.1 Join mappings by library type

Perform mapping for each library separately an then join in groups like these: gdna_zerob, gdna_plusb, rna_zerob, rna_plusb.

```
$ samtools merge -u - gdna_zerob_ind1.bam gdna_zerob_ind2.bam | samtools sort - gdna_zerob
$ samtools index gdna_zerob.bam
```

For the new Samtools version, the command are: 

```
$ samtools merge -u - *bam | samtools sort -T aln.sorted - -o project_h.bam
$ samtools index gdna_zerob.bam
```

### 2.2 Join counts

Get counts from the joined bam files and join in a single table: 

```
$ bam_var_join.py FastaReference ListOfBams
```

Output: toico3.txt

### 2.3 snp_calling_bchr.py: SNP calling

In a text file, indicate the type of library. At least, you should include gdna_zero (library withouth the chromosome) and gdna_plus (library with the chromosome). Please, use the same order like in the following example:

```
gdna_zerob.bam    gdna_zero
gdna_plusb.bam    gdna_plus
rna_zerob.bam    rna1_zero
rna_plusb.bam        rna1_plus
rna_zerob_leg.bam    rna2_zero
rna_plusb_leg.bam    rna2_plus
```
Then we perform the SNP calling with this new file and the previously generated "toico3.txt" file:

```
$ snp_calling_bchr.py baba.txt toico3.txt
```

Outputs: snps_selected2.txt: list with the selected SNPs
         ref_alt2.txt: list with the variant Ref (fixed in gdna_zerob) and variant Alt (private of gdna_plusb)

### 2.4 Get variation per library

This protocol works to get counts for Ref and Alt variants in the selected SNPs per individual library (not per joined libraries).

Firstly, we get a table with the countings for each nucleotide for all the libraries:

```
$ bam_var_join.py FastaReference ListOfBams
```

We create a file "sel.txt" with selected positions. We can addiontally apply additional filters. The "sel.txt" files looks like this:

```
seq1    3
seq1    84
seq2    456
seq3    123
```

Then we extract the positions from the ref_alt2.txt file:

```
$ grep -wFf sel.txt ref_alt2.txt > ref_alt3.txt
```

And from the  toico3.txt got from individual library:

```
$ grep -wFf sel.txt toico3.txt > toico4.txt
```

Finally, we get the counts.

```
$ get_var_library.py # modify line 20 with the number of libraries
```

### 2.5 Get Reference and alternative sequences 

```
$ sequence_ref_alt.py FastaFile ref_alt3.txt
```

You can also use "ref_alt2.txt" file, if you do not apply filters for SNPs selection.

## Protocol 3. Phylogenies

This protocol works for both transcriptomic and genomic libraries if we map libraries

```
$ bam_coverage_join.py FastaReference ListOfBams
```

Get majority consensus from BAM files for all the libraries of interest:

```
$ bam_consensus.py toico3.txt
```

Then we perform RAxML phylogenies with all genes in a FASTA file for a list of genes:

```
$ massive_phylogeny.py sequences_all_247genes_mod.fasta.fam gene_list.txt
```

Draw trees from the new trees:

```
$ massive_phylogeny_figure.py NewickList Outgroup
```
