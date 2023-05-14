# Su_Xingyuan_Summer_2023
Bioinformatics project on Rhizobium leguminosarum 

# Problems/Questions 

# Step 1 - Genome Assembly using SPAdes <br>
https://github.com/ablab/spades

## Information about Files (All files on Info server)
- **rhizo_ee**, directory, pathname ``/home/xingyuan/rhizo_ee``: 
   - directory for this project, rhizobium experimental evolution; environmental variable ``$RHIZO`` <br>
- **raw_reads**, directory, pathname ``/home/xingyuan/rhizo_ee/raw_reads``: 
   - 767 symbolic links linking to the raw reads in ``/scratch/batstonelab/RltEE2020-PE_reads``(an example of the file format is ``9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R2_001.fastq.gz@``); 
   - MultiQC report of the raw reads (``Project_Heath_363_gDNA_multiqc_report.html``)


## Progress 
**May 12, 2023** <br>
- strain 9_7_9 was assembled correctly. Type ``cd /home/xingyuan/rhizo_ee/raw_reads/9_7_9-spades`` to see results. Type ``less contigs.fasta``; result is ``>NODE_1_length_728405_cov_41.557482``, meaning 1 contig, 728405 nucleotides sequence length, and 41.557482 k-mer coverage for the last (largest) k value used.  Type ``less scaffolds.fasta``; result is ``NODE_1_length_728405_cov_41.557482``, same result as contig. 

### Pathnames
For original raw reads: ``/2/scratch/batstonelab/RltEE2020-PE_reads``

### Codes
``nohup spades.py --pe1-1 9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R1_001.fastq.gz --pe1-2 9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R2_001.fastq.gz -o 9_7_9-spades &``

### PATH
``.:/home/xingyuan/bin:/usr/lib64/qt-3.3/bin:/usr/lib64/openmpi/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/usr/local/python3/Python-3.5.1:/2/scratch/batstonelab/bin/SPAdes-3.15.5-Linux/bin:/home/xingyuan/.local/bin:/home/xingyuan/bin``

## Error codes 


