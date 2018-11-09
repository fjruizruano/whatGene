# whatGene
Tools to search for genes in special chromosomes like sexual or B chromosomes.

Repository under construction!!!!

## Installation
- Copy script to your binaries folder.
- Dependencies:
  * BioPython [http://biopython.org/wiki/Main_Page](http://biopython.org/wiki/Main_Page)
  * RepeatMasker [http://www.repeatmasker.org/RMDownload.html](http://www.repeatmasker.org/RMDownload.html)
  * seqTK [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)
  * DeconSeq [http://deconseq.sourceforge.net](http://deconseq.sourceforge.net)
  * BLAT [http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/)
  * Trimmomatic [http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

## 0. De novo assembly, CDS prediction and clustering

Trinity

Transdecoder

CD-HIT

## 1. Map reads against a transcriptomic reference

```
$ ssaha2_run_multi.py list.txt reference.fasta 12
```
or
```
$ ssaha2_run_multi_pe_se.py list.txt reference.fasta 12
```

## 2. Generate coverage file

```
$ bam_coverage_join.py FastaReference ListOfBams 999999]
```

## 3. Generate graphics

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


```
$ coverage_graphics.py CoverageFile SamplesFile FastaFile [PDF/SVG/NOPLOT]
```

