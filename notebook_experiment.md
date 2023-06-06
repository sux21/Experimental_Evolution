# Quality control by FastQC, MultiQC
## Unusual per base sequence content

### 19_1_9 has more AT than GC
<img width="600" alt="Screenshot 2023-06-01 at 4 56 26 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/82c07678-5004-4f19-9d31-d4665e9c2d1a">
<img width="600" alt="Screenshot 2023-06-01 at 4 56 48 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/31ccd674-0b64-4ec5-ade6-d42d0087a5b8">

### 19_4_7 has more similar amount of AT and GC
<img width="600" alt="Screenshot 2023-06-01 at 4 59 28 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/62f520d0-52e4-46dd-ba20-f3cfe9953aa3">
<img width="600" alt="Screenshot 2023-06-01 at 4 59 41 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/57bdbd42-76d8-4410-9f96-1d830f41bff9">

### as5_2_4 has more similar amount of AT and GC
<img width="600" alt="Screenshot 2023-06-01 at 5 07 46 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/0e5664b8-994f-4271-bf45-8bacc955abad">
<img width="600" alt="Screenshot 2023-06-01 at 5 08 07 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/8208affc-2244-4828-a11b-90d556c7e5d5">

## Failed per sequence gc content
See this website for explanation of per base gc content: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html. 

### 19_1_9
<img width="600" alt="Screenshot 2023-06-01 at 5 17 03 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/0a9ec26b-5c9a-4df0-8fc9-4397555bd411">
<img width="600" alt="Screenshot 2023-06-01 at 5 17 37 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/76780695-5ed8-49af-9eae-ca8532374fe2">

### 19_4_7
<img width="600" alt="Screenshot 2023-06-01 at 5 19 29 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/a6b269db-3c95-40d1-8eeb-57bae70f48a7">
<img width="600" alt="Screenshot 2023-06-01 at 5 19 48 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/22ae3924-4643-4e6e-9423-3b53eae82262">

## Failed Kmer content
See this website for explanation of Kmer content: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/11%20Kmer%20Content.html. 

### 7_1_5
<img width="600" alt="Screenshot 2023-06-01 at 5 25 15 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/0a15678c-8fc8-49d7-bc2d-b9112a52ef17">
<img width="600" alt="Screenshot 2023-06-01 at 5 26 00 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/c5a2ef12-b53d-4381-ba60-8d4f054717b4">

### 20_6_1
<img width="600" alt="Screenshot 2023-06-01 at 5 27 39 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/be9a43de-0704-41f5-bf86-aa5e16dcf751">
<img width="600" alt="Screenshot 2023-06-01 at 5 28 05 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/c8296a73-5156-4657-a1e8-3c8c04f6a1fc">

### 19_6_4
<img width="600" alt="Screenshot 2023-06-01 at 5 33 27 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/f5737bc9-50be-4d49-bedc-ae036243df9d">
<img width="600" alt="Screenshot 2023-06-01 at 5 33 48 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/1c1557af-9aea-411c-a635-0beb053c9a32">


### 15_2_2
<img width="600" alt="Screenshot 2023-06-01 at 5 28 50 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/2bb7eb05-5ee2-4f80-8fe2-33820f010768">
<img width="600" alt="Screenshot 2023-06-01 at 5 29 10 PM" src="https://github.com/sux21/Su_Xingyuan_Summer_2023/assets/132287930/e3c03b40-b403-4df6-a119-499ca7055b38">

# 19Genome assembly 
## Quast
### 52 samples from from 2020 strains in Rhizobium_leguminosarum_EE2021-Single_strain_experiment Google sheets
|#| sample | # scaffolds | total length | largest scaffold | GC (%) | N50 | N90 | auN | L50 | L90 | # N's per 100 kbp  |
|:--:| :----: | :----:    | :----:       | :---:          | :--:   | :--:|:---:|:--: |:--: |:--: | :--:               |
|1| 10_1_8 | 269 |  7698930 | 928458 | 60.61 | 400806 | 116845 | 436794.2 | 7 | 21 | 5.03 |
|2| 10_1_9 | 119 | 7377880  | 1009279 | 60.75 | 430333 | 173037 | 461108.1 | 6 | 18 |  7.94  |
|3| 10_7_6 | 229 | 7561496 | 966116 | 60.71 | 447903 |  86934 | 445840.7 | 7 | 19 |  5.20 | 
|4| 11_4_2 | 227 | 7490241 | 856458 | 60.69 | 334660 | 65380 |  364679.9 | 8 | 24 | 5.19 |
|5| 11_4_4 | 113 | 7526111 | 1446462 | 60.60 | 654969 | 159952 | 678186.6 | 4 | 14 |  6.52  |
|6| 11_5_6 | 136 | 7531401 | 1117631 | 60.60 | 590789 | 205988 | 654709.3  |  5  | 12 | 7.81 | 
|7| 13_4_1 | 216 | 7296160 | 1167888 | 60.75 | 500937 | 120765 | 535202.8 | 5 | 17 |  9.42 |
|8| 14_4_6 (has warning)| 150 | 7534545 | 1101709 | 60.60 | 654969 | 228788 | 677489.1 | 5 | 11 |  7.78 |
|9| 14_5_3 | 217 | 7299266 | 1165591 | 60.76 | 526898 | 120765 | 559835.2 | 5 | 15 | 9.43 |
|10| 15_4_4 | 274 | 7699117 | 970156 | 60.61 | 409389 |  105773  | 486538.2 | 6 | 19 | 6.34 | 
|11| 15_4_6 | 253 | 7654724 | 869297 | 60.62 | 445323 | 114888 | 453915.5 | 6 |  19 | 3.85 | 
|12| 16_1_6 | 350 | 7619799 | 899177 | 60.96 | 291640 | 50655 | 396217.0 | 8 | 28 | 2.73 |
|13| 16_1_7 | 270 | 7696264 | 970290 | 60.62 | 414354 | 105773 | 455952.9 | 6 | 21 | 7.68 |
|14| 16_1_8 | 286 | 7316283 | 1492913 | 60.70 | 517446 | 61754 | 582815.2 | 5 | 20 | 9.38 |
|15| 16_4_2 | 216 | 7561806 | 883671 | 60.70 | 459276 | 74856 | 424792.2 | 7 | 21 | 6.39 | 
|16| 16_4_3 | 425 | 7357352 | 1509109 | 60.69 | 508829 | 69558 | 586323.5 | 5 | 20 | 9.35 |
|17| 16_6_6 (has warning) | 171 | 6593588 | 627707 | 60.73 | 492317 | 74856 | 409345.3 | 6 | 17 | 4.41 | 
|18| 17_2_1 (has warning) | 215 | 7339141 |  979257 | 60.77 | 562328 | 147161 | 534068.6 | 6 | 14 | 6.69 |
|19| 17_2_8 | 270 | 7697265 | 928458 | 60.61 | 409389 | 105773 | 459457.4 | 6 | 20 | 5.09 |
|20| 17_2_9 | 229 | 7697644 | 1046772 | 60.61 | 436676 | 116844 | 517397.3 | 6 | 18 | 6.39 |
|21| 18_1_4 | 214 | 7341027 | 728361 | 60.77 | 573721 | 147161 | 485908.3 | 6 | 15 | 6.77 | 
|22| 18_1_5 (has warning) | 129 | 7162478 | 676748 | 60.79 | 389809 | 152413 | 397830.8 | 7 | 18 | 6.92 |
|23| 19_1_1 | 222 | 7334731 | 728361 | 60.77 | 534564 | 147161 | 465215.6 | 6 | 16 | 5.47 |
|24| 19_5_8 (has warning) | 200 | 7340052 | 728361 | 60.76 | 575221 | 147161 | 503422.5 | 6 | 14 | 7.97 |
|25| 2_2_5 | 259 | 7694619 | 928458 | 60.61 | 468807 | 89372 | 487128.9 | 6 | 21 | 6.34 | 
|26| 2_5_2 | 274 | 7314867 | 1101302 | 60.70 | 453203 | 68938 | 460247.1 | 6 | 20 | 10.69 |
|27| 2_6_4 | 298 | 7703285 | 970854 | 60.60 | 354826 | 116844 | 463980.2 | 6 | 21 | 5.05 | 
|28| 3_1_5 | 234 | 7313912 | 1493013 | 60.70 | 462326 | 85152 | 582249.4 | 6 | 19 | 10.84 |
|29| 3_2_1 | 213 | 7333973 | 728361 | 60.77 | 562328 | 147161 | 476943.0 | 6 | 15 | 6.68 | 
|30| 3_2_3 | 301 | 7097471 | 1046074 | 60.61 | 372426 | 99732 | 485678.2 | 6 | 19 | 4.17 |
|31| 3_2_6 | 187 | 7332725 | 728361 | 60.77 | 562175 | 147161 | 483739.7 | 6 | 15 | 8.11 |
|32| 3_2_7 | 175 | 7333153 | 980003 | 60.77 | 534902 | 147161 | 521483.3 | 6 | 15 | 8.14 |
|33| 3_3_5 | 226 | 7291290 | 1165711 | 60.76 | 456720 | 89373 | 516143.8 | 5 | 19 | 10.79 |
|34| 3_3_7 | 177 | 7332380 | 728361 | 60.78 | 575221 | 147161 | 514321.2 | 6 | 14 | 6.71 |
|35| 3_3_9 | 218 | 7487662 | 864148 | 60.69 | 334675 | 88616 | 366581.0 | 8 | 24 | 3.87 |
|36| 4_1_2 | 214 | 7313784 | 1492913 | 60.70 | 517446 | 90004 | 594910.5 | 5 | 18 | 9.43 |
|37| 4_1_4 | 206 | 7311227 | 1493013 | 60.70 | 453119 | 86038 | 569765.1 | 6 | 20 | 10.75 |
|38| 4_2_1 | 186 | 7292087 | 1165711 | 60.75 | 526898 | 120765 | 554707.6 | 5 | 16 | 9.43 |
|39| 6_4_5 | 250 | 7320356 | 1492674 | 60.71 | 517398 | 89372 | 589651.3 | 5 | 18 | 9.48 |
|40| 6_4_7 | 295 | 7698123 | 970156 | 60.61 | 409389 | 105773 | 492640.7 | 6 | 20 | 7.82 |
|41| 6_7_5 | 246 | 7298635 | 1165491 | 60.76 | 433102 | 104189 | 527662.5 | 6 | 17 | 8.07 |
|42| 7_1_2 | 163 | 7280581 | 979296 | 60.80 | 562328 | 147420 | 538716.8 | 6 | 14 | 8.21 |
|43| 7_1_5 | 215 | 7340092 | 805435 | 60.78 | 575220 | 95675 | 505635.1 | 6 | 15 | 6.73 |
|44| 7_6_3 |
