# Genome assembly
## Quast outputs
Quast outputs the followings information in a ``report.txt`` file:
- [# contigs (>= ...bp)](https://quast.sourceforge.net/docs/manual.html#sec3)
- [Total length (>= ...bp)](https://quast.sourceforge.net/docs/manual.html#sec3)
- [Largest contig](https://quast.sourceforge.net/docs/manual.html#sec3)
- [Total length](https://quast.sourceforge.net/docs/manual.html#sec3)
- [GC (%)](https://quast.sourceforge.net/docs/manual.html#sec3)
- [N50, N90](https://quast.sourceforge.net/docs/manual.html#sec3)
- [auN](https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity)
- [L50, L90](https://quast.sourceforge.net/docs/manual.html#sec3)
- [# N's per 100 kbp](https://quast.sourceforge.net/docs/manual.html#sec3)
<img width="348" alt="Screenshot 2023-06-02 at 4 24 50 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/a29b3cbd-13c3-40a1-8d21-b304fd7a7ef3">

## Qualimap 
BAM QC outputs basic statistic for an alignment in a ``qualimapReport.html`` file. Explanations of these outputs below are copied from http://qualimap.conesalab.org/doc_html/analysis.html#output. 
- Globals: This section contains information about the total number of reads, number of mapped reads, paired-end mapping performance, read length distribution, number of clipped reads and duplication rate (estimated from the start positions of read alignments). <br>
<img width="381" alt="Screenshot 2023-06-03 at 11 24 39 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/9e03f341-04cc-45b3-af9c-f9110c8ee4ad"> <br>
- ACGT Content: Nucleotide content and GC percentage in the mapped reads. <br>
<img width="356" alt="Screenshot 2023-06-03 at 11 26 29 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/2c2bcefe-867d-4645-9e07-e5b6131e5c37"> <br>
- Coverage: Mean and standard deviation of the coverage depth. <br>
<img width="226" alt="Screenshot 2023-06-03 at 11 26 47 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/4f4d0885-299c-4e89-bafd-ac35bd06b09e"> <br>
- Mapping Quality: Mean mapping quality of the mapped reads. <br>
<img width="218" alt="Screenshot 2023-06-03 at 11 27 06 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/14840b8e-c037-484d-8b20-e59db70fe57b"> <br>
- Insert size: Mean, standard deviation and percentiles of the insert size distribution if applicable. The features are computed based on the TLEN field of the SAM file. <br>
<img width="256" alt="Screenshot 2023-06-03 at 11 27 24 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/55031796-e1b8-4a56-831d-0db84d6625d7"> <br>
- Mismatches and indels: The section reports general alignment error rate (computed as a ratio of total collected edit distance to the number of mapped bases), total number of mismatches and total number of indels (computed from the CIGAR values). Additionally fraction of the homopolymer indels among total indels is provided. Note, the error rate and mismatches metrics are based on optional fields of a SAM record (NM for edit distance, MD for mismatches). The features are not reported if these fields are missing in the SAM file. <br>
<img width="360" alt="Screenshot 2023-06-03 at 11 27 45 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/d59cb375-193c-4009-b50c-4fa79b710e53"> <br>
- Chromosome stats: Number of mapped bases, mean and standard deviation of the coverage depth for each chromosome as defined by the header of the SAM file. <br>
<img width="400" alt="Screenshot 2023-06-03 at 11 28 25 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/6089b6b0-9c14-4703-be3b-e7afed1fd326"> <br>
- Coverage across reference: This plot consists of two figures. The upper figure provides the coverage distribution (red line) and coverage deviation across the reference sequence. The coverage is measured in X [1]. The lower figure shows GC content across reference (black line) together with its average value (red dotted line). <br>
[1]	Example for the meaning of X: If one genomic region has a coverage of 10X, it means that, on average, 10 different reads are mapped to each nucleotide of the region. <br>
<img width="400" alt="Screenshot 2023-06-03 at 11 30 31 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/e9b30040-623d-45eb-b9ce-2dfef7c5f56d"> <br>
- Coverage Histogram: Histogram of the number of genomic locations having a given coverage rate. The bins of the x-axis are conveniently scaled by aggregating some coverage values in order to produce a representative histogram also in presence of the usual NGS peaks of coverage. <br> 
<img width="400" alt="Screenshot 2023-06-03 at 11 34 21 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/c2806e7d-803c-4247-82c4-b813baf52f6a"> <br>
- Coverage Histogram (0-50X): Histogram of the number of genomic locations having a given coverage rate. In this graph genome locations with a coverage greater than 50X are grouped into the last bin. By doing so a higher resolution of the most common values for the coverage rate is obtained. <br>
<img width="400" alt="Screenshot 2023-06-03 at 11 37 38 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/0e0e5bdb-5993-4604-946b-aece10fd57c2"> <br>
- Genome Fraction Coverage: Provides a visual way of knowing how much reference has been sequenced with at least a given coverage rate. This graph should be interpreted as in this example: <br>
If one aims a coverage rate of at least 25X (x-axis), how much of reference (y-axis) will be considered? The answer to this question in the case of the whole-genome sequencing provided example (copied below) is ~83%. <br>
<img width="400" alt="Screenshot 2023-06-03 at 11 41 13 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/d9b6fff2-bbfd-470f-8a57-5c5b66a6ff36"> <br>
- Duplication Rate Histogram: This plot shows the distribution of duplicated read starts. Due to several factors (e.g. amount of starting material, sample preparation, etc) it is possible that the same fragments are sequenced several times. For some experiments where enrichment is used (e.g. ChIP-seq ) this is expected at some low rate. If most of the reads share the exact same genomic positions there is very likely an associated bias. <br>
<img width="400" alt="Screenshot 2023-06-03 at 11 45 18 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/2c57f651-5ac4-4e54-a685-77a6cd896c45"> <br>
- Mapped Reads Nucleotide Content: This plot shows the nucleotide content per position of the mapped reads. <br>
<img width="400" alt="Screenshot 2023-06-03 at 11 47 45 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/172411ef-27d0-4377-a562-464b2e5641d4"> <br>
- Mapped Reads GC-content Distribution: This graph shows the distribution of GC content per mapped read. If compared with a precomputed genome distribution, this plot allows to check if there is a shift in the GC content. <br>
<img width="400" alt="Screenshot 2023-06-03 at 11 49 13 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/738d08c3-7e84-41b0-ac04-f50358e78392"> <br>
- Mapped Reads Clipping Profile: Represents the percentage of clipped bases across the reads. The clipping is detected via SAM format CIGAR codes ‘H’ (hard clipping) and ‘S’ (soft clipping). In addition, the total number of clipped reads can be found in the report Summary. The plot is not shown if there are no clipped-reads are found. Total number of clipped reads can be found in Summary. <br>
<img width="400" alt="Screenshot 2023-06-03 at 11 52 03 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/acf2bfb9-f8a0-4836-8ec6-e310a8148be1"> <br>
- Homopolymer Indels: This bar plot shows separately the number of indels that are within a homopolymer of A’s, C’s, G’s or T’s together with the number of indels that are not within a homopolymer. Large numbers of homopolymer indels may indicate a problem in a sequencing process. An indel is considered homopolymeric if it is found within a homopolymer (defined as at least 5 equal consecutive bases). Owing to the fact that Qualimap works directly from BAM files (and not from reference genomes), we make use of the CIGAR code from the corresponding read for this task. Indel statistics cam be found in a dedicated section of the Summary report. <br>
This chart is not shown if the sample doesn’t contain any indels. <br>
<img width="400" alt="Screenshot 2023-06-03 at 11 54 58 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/80a971b6-0870-4c92-acf8-367716ea0e80"> <br>
- Mapping Quality Across Reference: This plot provides the mapping quality distribution across the reference. To construct the plot mean mapping quality is computed for each window. <br>
<img width="400" alt="Screenshot 2023-06-03 at 11 56 36 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/bb05a380-c3bf-43df-a3d6-aa0b63cce2b4"> <br>
- Mapping Quality Histogram: Histogram of the number of genomic locations having a given mapping quality. To construct the histogram mean mapping quality is computed at each genome position with non-zero coverage and collected. According to Specification of the SAM format the range for the mapping quality is [0-255]. [SAMv1.pdf](https://github.com/sux21/Su_Xingyuan_Summer_2023/files/11644726/SAMv1.pdf) <br>
<img width="400" alt="Screenshot 2023-06-04 at 12 02 23 AM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/50621f0a-0d00-4b8f-b8d1-74dc65bfb71c"> <br>
- Insert Size Across Reference: This plot provides the insert size distribution across the reference. Insert size is collected from the SAM alignment field TLEN. Only positive values are taken into account. To construct the plot mean insert size is computed for each window. <br>
<img width="400" alt="Screenshot 2023-06-04 at 12 06 10 AM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/ccf58f76-a31e-4b4d-be35-60f06721e03d"> 
<img width="310" alt="Screenshot 2023-06-04 at 12 34 22 AM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/3b9ec814-7017-4799-9ce4-2e08c694f40f"> <br>
- Insert Size Histogram: Histogram of insert size distribution. To construct the histogram all collected insert size values (number of read alignments of a certain insert size) are seprated in bins. The default number of bins is 50. The number of reads in the bin will be the sum from insert size of the bins. The detailed values for insert are reported in raw data report. <br>
<img width="400" alt="Screenshot 2023-06-04 at 12 09 59 AM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/6c5a5faa-220f-4120-abc1-08b2322ebec2"> 
<img width="200" alt="Screenshot 2023-06-04 at 12 32 31 AM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/04846ff3-ac33-4b6b-b988-fea232d61b60">


