# Su_Xingyuan_May_to_August_2023
Bioinformatics project on Rhizobium leguminosarum 

# Step 1 - Genome Assembly using SPAdes <br>
version: SPAdes genome assembler v3.15.2

## Information about Files (All files on Info server)
- **rhizo_ee**, directory, pathname ``/home/xingyuan/rhizo_ee``: 
   - directory for this project, rhizobium experimental evolution; environmental variable ``$RHIZO`` <br>
- **raw_reads**, directory, pathname ``/home/xingyuan/rhizo_ee/raw_reads``: 
   - 767 symbolic links linking to the raw reads in ``/scratch/batstonelab/RltEE2020-PE_reads``(an example of the file format is ``9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R2_001.fastq.gz@``); 
   - MultiQC report of the raw reads (``Project_Heath_363_gDNA_multiqc_report.html``)

## Before Assembly
### 1. Running FastQC (Practice on 2018 strains)
FastQC v0.11.5

Input file format is ``.fastq``. An example is ``GSF2234-101A_S1_R1_001.fastq``.

Running FastQC: ``nohup fastqc -o /home/xingyuan/2018_strains/fastQC_raw_reads GSF2234-295A_S28_R2_001.fastq &`` 

``fastqc -o /home/xingyuan/2018_strains/fastQC_trimmed_reads *gz``

Output files are ``fastqc.html`` and ``fastqc.zip``. An example is ``GSF2234-101A_S1_R1_001_fastqc.html`` and ``GSF2234-101A_S1_R1_001_fastqc.zip``.

Commands are from https://home.cc.umanitoba.ca/~psgendb/doc/fastqc.help. 

### 2. Running MultiQC (Practice on 2018 strains) 
multiqc, version 1.9

Running MultiQC: ``multiqc .`` in the directory with the FastQC reports 

### 3. Running Trimmomatic (Practice on 2018 strains)
version: trimmomatic-0.39.jar

http://www.usadellab.org/cms/?page=trimmomatic

Modify this command: ``java -jar /path/to/trimmomatic.jar PE R1_001.fastq \ R2_001.fastq R1_P.fq.gz R1_UP.fq.gz R2_P.fq.gz R2_UP.fq.gz \ ILLUMINACLIP:/path/to/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE``

``java -jar /usr/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE /home/xingyuan/2018_strains/raw_reads/GSF2234-101A_S1_R1_001.fastq /home/xingyuan/2018_strains/raw_reads/GSF2234-101A_S1_R2_001.fastq /home/xingyuan/2018_strains/trimmed_reads/GSF2234-101A_S1_R1_P_001.fq.gz /home/xingyuan/2018_strains/trimmed_reads/GSF2234-101A_S1_R1_UP_001.fq.gz /home/xingyuan/2018_strains/trimmed_reads/GSF2234-101A_S1_R2_P_001.fq.gz /home/xingyuan/2018_strains/trimmed_reads/GSF2234-101A_S1_R2_UP_001.fq.gz ILLUMINACLIP:/usr/local/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE``

Running for all files using shell script: <br>
``#!/bin/bash`` <br>
``for f in $(ls *.fastq | sed 's/?_001.fastq//' | sort -u)`` <br>
``do`` <br>
``java -jar /usr/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE /home/xingyuan/2018_strains/raw_reads/${f}1_001.fastq /home/xingyuan/2018_strains/raw_reads/${f}2_001.fastq /home/xingyuan/2018_strains/trimmed_reads/${f}1_P_001.fq.gz /home/xingyuan/2018_strains/trimmed_reads/${f}1_UP_001.fq.gz /home/xingyuan/2018_strains/trimmed_reads/${f}2_P_001.fq.gz /home/xingyuan/2018_strains/trimmed_reads/${f}2_UP_001.fq.gz ILLUMINACLIP:/usr/local/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE`` <br>
``done``

### 4. Repeat 1 and 2 (Practice on 2018 strains)




## During Assembly 
https://www.melbournebioinformatics.org.au/tutorials/tutorials/assembly/spades/

## After Assembly 

**Output** <br>
``contigs.fasta``: you can check the number of contigs by counting the ``>`` symbol, length of each contig, k-mer coverage of the largest k-value used (k-mer coverage is always lower than read coverage). See [Contigs and scaffolds format](https://github.com/ablab/spades#contigs-and-scaffolds-format). 

## Progress 
**May 12, 2023** <br>
- strain 9_7_9 was assembled correctly. Type ``cd /home/xingyuan/rhizo_ee/raw_reads/9_7_9-spades`` to see results. <br>
     - Methods: ``Command line: /usr/local/spades/version.3.15.2/bin/spades.py --pe1-1 /home/xingyuan/rhizo_ee/raw_reads/9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R1_001.fastq.gz --pe1-2 /home/xingyuan/rhizo_ee/raw_reads/9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R2_001.fastq.gz -o   /home/xingyuan/rhizo_ee/raw_reads/9_7_9-spades``, ``System information: SPAdes version: 3.15.2 Python version: 3.5.1 OS: Linux-2.6.32-754.30.2.el6.x86_64-x86_64-with-centos-6.10-Final``, 
     - Results: 335 contigs, 333 scaffolds. 

### Pathnames
For original raw reads: ``/2/scratch/batstonelab/RltEE2020-PE_reads``

### Codes
``nohup spades.py --pe1-1 9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R1_001.fastq.gz --pe1-2 9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R2_001.fastq.gz -o 9_7_9-spades &``

``Bandage image assembly_graph_with_scaffolds.gfa assembly_graph_with_scaffolds.jpg``

``quast.py contigs.fasta``

``scp xingyuan@info.mcmaster.ca:/home/xingyuan/2018_strains/fastQC_raw_reads/GSF2234-101A_S1_R1_001_fastqc.html /Users/xingyuansu/Desktop``

## Error codes 

