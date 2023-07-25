
## After Assembly 
### 1. Run BWA
Version: 0.7.17-r1188 

#### 214C, 218A, 272A, 295A
1. Index the contigs/scaffolds
```
bwa index /home/xingyuan/2018_strains/trim_2nd_attempt/spades-295A/scaffolds.fasta     
```
2. align reads back to contigs/scaffolds
```
bwa mem -t 2 /home/xingyuan/2018_strains/trim_2nd_attempt/spades-295A/scaffolds.fasta /home/xingyuan/2018_strains/trim_2nd_attempt/GSF2234-295A_S28_R1_P_001.fq.gz /home/xingyuan/2018_strains/trim_2nd_attempt/GSF2234-295A_S28_R2_P_001.fq.gz > 295A_scaffolds.sam 
```
3. convert SAM file to BAM file (see commands [here](http://www.htslib.org/doc/samtools-view.html))
```
samtools view -bS 295A_scaffolds.sam > 295A_scaffolds.bam 
```
4. sort the BAM file (see commands [here](http://www.htslib.org/doc/samtools-sort.html))
```
samtools sort -o 295A_scaffolds.sorted.bam 295A_scaffolds.bam 
```
5. index the BAM file (see commands [here](http://www.htslib.org/doc/samtools-index.html))
```
samtools index 295A_scaffolds.sorted.bam 
```
6. obtain summary statistics (see commands [here](http://www.htslib.org/doc/samtools-flagstat.html))
```
samtools flagstat 295A_scaffolds.sorted.bam > mapping_summary 
```
7. Run qualimap for more information (see commands [here](http://qualimap.conesalab.org/doc_html/analysis.html))
```
qualimap bamqc -outdir bamqc -bam 295A_scaffolds.sorted.bam 
```
