# Su_Xingyuan_Summer_2023
Bioinformatics project on Rhizobium leguminosarum 

# Questions
- The pathname of the original raw reads is ``/2/scratch/batstonelab/RltEE2020-PE_reads``, different from before. Change it if the SPAdes fails.

# Step 1 - Genome Assembly using SPAdes <br>
https://github.com/ablab/spades

## Information about Files (All files on Info server)
- **rhizo_ee**, directory, pathname ``/home/xingyuan/rhizo_ee``: 
   - directory for this project, rhizobium experimental evolution; environmental variable ``$RHIZO`` <br>
- **raw_reads**, directory, pathname ``/home/xingyuan/rhizo_ee/raw_reads``: 
   - 767 symbolic links linking to the raw reads in ``/scratch/batstonelab/RltEE2020-PE_reads``(an example of the file format is ``9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R2_001.fastq.gz@``); 
   - MultiQC report of the raw reads (``Project_Heath_363_gDNA_multiqc_report.html``)


## Progress 

### Codes
$SPADES --pe1-1 9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R1_001.fastq.gz --pe1-2 9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R2_001.fastq.gz -o 9_7_9-spades

/usr/local/python3/Python-3.5.1
/usr/local/spades/version.3.15.5/assembler/spades.py 
