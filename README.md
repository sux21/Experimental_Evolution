# Su_Xingyuan_May_to_August_2023
Bioinformatics project on Rhizobium leguminosarum 

# Problems/Questions 

# Step 1 - Genome Assembly using SPAdes <br>
version: SPAdes genome assembler v3.15.2

https://github.com/ablab/spades

## Information about Files (All files on Info server)
- **rhizo_ee**, directory, pathname ``/home/xingyuan/rhizo_ee``: 
   - directory for this project, rhizobium experimental evolution; environmental variable ``$RHIZO`` <br>
- **raw_reads**, directory, pathname ``/home/xingyuan/rhizo_ee/raw_reads``: 
   - 767 symbolic links linking to the raw reads in ``/scratch/batstonelab/RltEE2020-PE_reads``(an example of the file format is ``9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R2_001.fastq.gz@``); 
   - MultiQC report of the raw reads (``Project_Heath_363_gDNA_multiqc_report.html``)

## Before Assembly

## During Assembly 

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

``scp xingyuan@info.mcmaster.ca:/home/xingyuan/rhizo_ee/raw_reads/9_7_9-spades/quast_results/results_2023_05_16_10_24_26/report.html /Users/xingyuansu/Desktop``


## Error codes 

