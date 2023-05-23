# Su_Xingyuan_May_to_August_2023
Bioinformatics project on Rhizobium leguminosarum 

# Step 0 - Background Information of Illumina Sequencing
Paired-End Sequencing: https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html.

# Step 1 - Genome Assembly using SPAdes <br>

## Before Assembly
### 1. Running FastQC for raw data (Practice on 2018 strains)
FastQC v0.11.5

Running FastQC: 
```
nohup fastqc -o /home/xingyuan/2018_strains/fastQC_raw_reads GSF2234-295A_S28_R2_001.fastq &
```

Input files in /home/xingyuan/2018_strains/raw_reads. <br>
Output files in /home/xingyuan/2018_strains/fastQC_raw_reads.

### 2. Running MultiQC for raw data (Practice on 2018 strains) 
multiqc, version 1.9

Running MultiQC in the directory with the FastQC reports: 
```
multiqc .
```

Input files in /home/xingyuan/2018_strains/fastQC_raw_reads. <br>
Output files in /home/xingyuan/2018_strains/fastQC_raw_reads.

### 3. Running Trimmomatic (Practice on 2018 strains)
version: trimmomatic-0.39.jar

Running for all files using for loop in shell script: <br>
```
#!/bin/bash 
for R1 in *R1* 
do 
R2=${R1//R1_001.fastq/R2_001.fastq} 
R1_P=${R1//001.fastq/P_001.fq.gz} 
R1_UP=${R1//001.fastq/UP_001.fq.gz} 
R2_P=${R2//001.fastq/P_001.fq.gz} 
R2_UP=${R2//001.fastq/UP_001.fq.gz} 

java -jar /usr/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE /home/xingyuan/2018_strains/raw_reads/$R1 /home/xingyuan/2018_strains/raw_reads/$R2 /home/xingyuan/2018_strains/trim_2nd_attempt/$R1_P /home/xingyuan/2018_strains/trim_2nd_attempt/$R1_UP /home/xingyuan/2018_strains/trim_2nd_attempt/$R2_P /home/xingyuan/2018_strains/trim_2nd_attempt/$R2_UP ILLUMINACLIP:/usr/local/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE HEADCROP:15 
done
```

Input files in /home/xingyuan/2018_strains/fastQC_raw_reads. <br>
Output files in /home/xingyuan/2018_strains/trimmed_reads.

### 4. Run FastQC again for trimmed data (Practice on 2018 strains)
Run FastQC: 
```
nohup fastqc -o /home/xingyuan/2018_strains/fastQC_trim_2nd_attempt *_P_* &
```

Input files in /home/xingyuan/2018_strains/trimmed_reads. <br>
Output files in /home/xingyuan/2018_strains/fastQC_trimmed_reads.

### 5. Run MultiQC again for trimmed data (Practice on 2018 strains)
Run MultiQC in the directory with the FastQC reports: 
```
multiqc -f . 
```
Input files in /home/xingyuan/2018_strains/fastQC_trimmed_reads. <br>
Output files in /home/xingyuan/2018_strains/fastQC_trimmed_reads.

### 6. Steps 3, 4, 5 were repeated one time to remove the first 15bp of the read.

## During Assembly 
### 1. Running SPAdes (Practice on 2018 strains)
version: SPAdes genome assembler v3.15.2
https://www.melbournebioinformatics.org.au/tutorials/tutorials/assembly/spades/

#### For 101A <br>
**Run SPAdes on trimmed reads**: 
```
nohup spades.py --pe1-1 GSF2234-101A_S1_R1_P_001.fq.gz --pe1-2 GSF2234-101A_S1_R2_P_001.fq.gz -o spades-101A &
```

```
Quast output:
Assembly                   scaffolds
# contigs (>= 0 bp)        146      
# contigs (>= 1000 bp)     68       
Total length (>= 0 bp)     7138120  
Total length (>= 1000 bp)  7114449  
# contigs                  81       
Largest contig             553567   
Total length               7123407  
GC (%)                     60.85    
N50                        190187   
N75                        104963   
L50                        12       
L75                        25       
# N's per 100 kbp          8.70 
```

**Run SPAdes on trimmed reads - using --isolate**:
```
spades.py -1 GSF2234-101A_S1_R1_P_001.fq.gz -2 GSF2234-101A_S1_R2_P_001.fq.gz --isolate -o spades-101A-test
```

```
Quast results:
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        210      
# contigs (>= 1000 bp)     49       
Total length (>= 0 bp)     7155610  
Total length (>= 1000 bp)  7128413  
# contigs                  52       
Largest contig             564345   
Total length               7130378  
GC (%)                     60.85    
N50                        226578   
N75                        127550   
L50                        10       
L75                        20       
# N's per 100 kbp          8.55  
```
**Test the following command and compare with the above to see if there are differences**
```
spades.py --pe1-1 GSF2234-101A_S1_R1_P_001.fq.gz --pe1-2 GSF2234-101A_S1_R2_P_001.fq.gz --isolate -o spades-101A-test2 
```

Do the same for these strains: <br>
051: `` `` <br>
125: ``nohup spades.py --pe1-1 GSF2234-101A_S1_R1_P_001.fq.gz --pe1-2 GSF2234-101A_S1_R2_P_001.fq.gz -o spades-125A &`` <br>
144D: ``nohup spades.py --pe1-1  --pe1-2  -o spades- &``
153C: ``nohup spades.py --pe1-1 GSF2234-153C_S7_R1_P_001.fq.gz --pe1-2 GSF2234-153C_S7_R1_UP_001.fq.gz -o spades-153C &`` <br>
154B: ``nohup spades.py --pe1-1 GSF2234-154B_S8_R1_P_001.fq.gz --pe1-2 GSF2234-154B_S8_R2_P_001.fq.gz -o spades-154B &`` <br>
164A: ``nohup spades.py --pe1-1 GSF2234-164A_S9_R1_P_001.fq.gz --pe1-2 GSF2234-164A_S9_R2_P_001.fq.gz -o spades-164A &`` <br>
177: ``nohup spades.py --pe1-1 GSF2234-177C_S11_R1_P_001.fq.gz --pe1-2 GSF2234-177C_S11_R2_P_001.fq.gz -o spades-177C &`` <br>
214C: ``nohup spades.py --pe1-1 GSF2234-214C_S15_R1_P_001.fq.gz --pe1-2 GSF2234-214C_S15_R2_P_001.fq.gz -o spades-214C &`` <br>
218A: ``nohup spades.py --pe1-1 GSF2234-218A_S16_R1_P_001.fq.gz --pe1-2 GSF2234-218A_S16_R2_P_001.fq.gz -o spades-218A &`` <br>
225A: ``nohup spades.py --pe1-1 GSF2234-225A_S17_R1_P_001.fq.gz --pe1-2 GSF2234-225A_S17_R2_P_001.fq.gz -o spades-225A &`` <br>
272A: ``nohup spades.py --pe1-1 GSF2234-272A_S22_R1_P_001.fq.gz --pe1-2 GSF2234-272A_S22_R2_P_001.fq.gz -o spades-272A &`` <br>
295A: ``nohup spades.py --pe1-1 GSF2234-295A_S28_R1_P_001.fq.gz --pe1-2 GSF2234-295A_S28_R2_P_001.fq.gz -o spades-295A &`` <br>
298A: ``nohup spades.py --pe1-1 GSF2234-298A_S30_R1_P_001.fq.gz --pe1-2 GSF2234-298A_S30_R2_P_001.fq.gz -o spades-298A &`` <br>
301D: ``nohup spades.py --pe1-1  --pe1-2  -o spades- &`` <br>
336: ``nohup spades.py --pe1-1 GSF2234-336A_S39_R1_P_001.fq.gz --pe1-2 GSF2234-336A_S39_R2_P_001.fq.gz -o spades-336A &`` <br>
338A: ``nohup spades.py --pe1-1  --pe1-2  -o spades- &`` <br>
372: ``nohup spades.py --pe1-1  --pe1-2  -o spades- &`` <br>
377: ``nohup spades.py --pe1-1  --pe1-2  -o spades- &`` <br>
377A: ``nohup spades.py --pe1-1  --pe1-2  -o spades- &`` <br>
391: ``nohup spades.py --pe1-1  --pe1-2  -o spades- &`` <br>
524D: ``nohup spades.py --pe1-1  --pe1-2  -o spades- &`` <br>

## After Assembly 

### Codes

``spades.py -1 GSF2234-101A_S1_R1_P_001.fq.gz -2 GSF2234-101A_S1_R2_P_001.fq.gz --careful -o spades-101A-test``

``Bandage image assembly_graph_with_scaffolds.gfa assembly_graph_with_scaffolds.jpg``

``quast.py contigs.fasta``

``scp xingyuan@info.mcmaster.ca:/home/xingyuan/2018_strains/fastQC_trimmed_reads/GSF2234-105A_S2_R2_P_001_fastqc.html /Users/xingyuansu/Desktop``

## Error codes 

