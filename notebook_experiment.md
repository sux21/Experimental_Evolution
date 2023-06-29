# Quality after trimming
## Unusual per base sequence content

### 19_1_9 has more AT than GC
<img width="700" alt="Screenshot 2023-06-26 at 11 34 31 AM" src="https://github.com/sux21/2020_Experimental_Evoluntion/assets/132287930/575a7592-2d5e-4180-979d-b79e98daacb7"> <br>
<img width="700" alt="Screenshot 2023-06-26 at 11 35 29 AM" src="https://github.com/sux21/2020_Experimental_Evoluntion/assets/132287930/e1b26e4c-03db-48f4-82ce-aa95deeb8f57">

### 19_4_7 has more similar amount of AT and GC
<img width="700" alt="Screenshot 2023-06-26 at 11 38 24 AM" src="https://github.com/sux21/2020_Experimental_Evoluntion/assets/132287930/58445e3d-6170-469e-a75f-aabbc4b660d0"> <br>
<img width="700" alt="Screenshot 2023-06-26 at 11 38 42 AM" src="https://github.com/sux21/2020_Experimental_Evoluntion/assets/132287930/b8584e8b-2f28-4051-8d49-845a590ba41e">

### as5_2_4 has more similar amount of AT and GC
<img width="700" alt="Screenshot 2023-06-26 at 11 40 55 AM" src="https://github.com/sux21/2020_Experimental_Evoluntion/assets/132287930/a74ddee3-1444-498b-b948-9916aa826291"> <br>
<img width="700" alt="Screenshot 2023-06-26 at 11 41 16 AM" src="https://github.com/sux21/2020_Experimental_Evoluntion/assets/132287930/fd8e1115-d260-4e10-8eb5-d7fef60b58eb">

## Failed per sequence gc content
See this website for explanation of per base gc content: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html. 

### 19_1_9
<img width="700" alt="Screenshot 2023-06-26 at 11 35 57 AM" src="https://github.com/sux21/2020_Experimental_Evoluntion/assets/132287930/4547f8e3-5018-4042-868b-9f0fdd0f7a53"> <br>
<img width="700" alt="Screenshot 2023-06-26 at 11 36 18 AM" src="https://github.com/sux21/2020_Experimental_Evoluntion/assets/132287930/6fbb1454-8a96-4f1e-84cf-b3b54a40db6a">

### 19_4_7
<img width="700" alt="Screenshot 2023-06-26 at 11 39 57 AM" src="https://github.com/sux21/2020_Experimental_Evoluntion/assets/132287930/aaafab26-49ca-4da2-b849-d49c492a4690"> <br>
<img width="700" alt="Screenshot 2023-06-26 at 11 40 16 AM" src="https://github.com/sux21/2020_Experimental_Evoluntion/assets/132287930/e300edc7-13e1-4c24-8d91-46413216f368">

## Failed Kmer content
See this website for explanation of Kmer content: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/11%20Kmer%20Content.html. 

### 7_1_5

### 20_6_1

### 19_6_4

### 15_2_2

# Pairwise comparison of average nucleotide identity results
Reference to reference (comparing 56 original strains to itself)
- Comparing same strain to same strain: 
- Comparing Rht_173_C to Rht_209_N: IQ-Tree shows they are identical, but ANI value is 99.9999. Don't know why they are not 100.
- Comparing as_5_2_4 to all: no ANI values using FastANI default settings. Top hits of blast search for the first five contigs are *Paenibacilus*. 






