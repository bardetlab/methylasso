# MethyLasso
A segmentation approach to analyze DNA methylation patterns and identify differentially methylation regions from whole-genome datasets.

MethyLasso models DNA methylation data using a nonparametric regression framework known as a Generalized Additive Model. It relies on the fused lasso method to segment the genome by estimating regions in which methylation is constant to study DNA methylation patterns in a single condition or capture dynamic changes of DNA methylation levels across conditions.  
MethyLasso identifies low-methylated regions (LMRs), unmethylated regions (UMRs), DNA methylation valleys (DMVs) and partially methylated domains (PMDs) in a single condition as well as differentially methylated regions (DMRs) between two conditions.

By Delphine Balaramane, Yannick G Spill, Michaël Weber & Anaïs F Bardet.  

Publication: https://www.biorxiv.org/content/10.1101/2023.07.27.550791v1

## Important fix for bismark input files
MethyLasso version 1.0.0 from December 12th 2024 should be used

## 1. &emsp;Installation
Dependencies:
- Requires R 3.6 or greater
- Requires R packages: Rcpp, RcppEigen, methods,scales, Matrix, matrixStats, data.table, ggplot2, foreach, doParallel, stringr, R.utils
- Requires conda (if you do the conda installation)

mkdir MethyLasso  
unzip methylasso-main.zip  
R CMD INSTALL --preclean methylasso-main  

Conda installation:  

unzip methylasso-main.zip  
conda env create -f methylasso-main/conda_env.yml  
conda activate MethyLasso  
R CMD INSTALL --preclean methylasso-main
  
  
## 2. &emsp;Usage
Using data from two different conditions:  
Rscript MethyLasso.R --n1 [name1] --c1 [condition1_replicate1,condition1_replicate2,..] --n2 [name2] --c2 [condition2_replicate1,condition2_replicate2,..] [options]

Using data from only one condition (not calling DMRs):  
Rscript MethyLasso.R --n1 [name1] --c1 [condition1_replicate1,condition1_replicate2,..] [options]
  
  
## 3. &emsp;Arguments

|**argument**|**description**|
| ---------- | --- |
| --n1|Name of condition 1 (used in output file names)   |
| --c1|Input file(s) for condition 1 (separated by comma without space for replicates) |
| --n2|Name of condition 2 (used in output file names and corresponding to reference condition for DMR identification) | 
| --c2|Input file(s) for condition 2 (separated by comma without space for replicates) |  
  
## 4. &emsp;Input format

By Default, Methylasso takes a tab-delimited file containing for each CpG:  
chr / start / end / percent_methylation / count_methylated / count_unmethylated 
This file can be generated within the Bismark pipeline with the command Coverage2Cytosine with option --mergeCpG.  

Alternatively, files with different formats can be provided by specifying the columns. For each CpG, their genomic coordinates should be in the first three columns (chr / start / end) and columns containing the counts of methylated and unmethylated Cs OR total count and methylation percentage should be specified.  

|**argument**|**description**|
| ---------- | --- |
| --mC | Column number containing count of methylated Cs  |
| --uC | Column number containing count of unmethylated Cs   |
OR
| --cov | Column number containing total count |
| --meth | Column number containing percentage of methylation [0-1] |

    
## 5. &emsp;Options  

|**argument**|**description**|
| ---------- | --- |
|-c | Minimum coverage per CpG (default 5)|
|-n | Minimum number of CpGs in LMRs, UMRS, DMVs and DMRs (default 4)|
|-m | Maximum mean methylation for PMDs (default 0.7)|
|-d | Minimum methylation difference for DMRs [0-1] (default 0.1)|
|-p | P-value threshold for DMRs (default 0.05)|
|-q | Q-value (FDR) threshold for DMRs (to be used instead of p-value threshold)|
|-r | Replicate coverage score for CpGs in DMRs (default 0.7)|
|-s | Skip LMR, UMR, DMV and PMD identification (default FALSE)|
|-t | Number of threads to use (default 1)|
|-o | Output directory (default current directory)|
|-f | Output figures in pdf format (default TRUE)|
|--quiet | Do not print processing information (default FALSE)|
|--version | Print version|
|--help | Print help|

## 6. &emsp;Additional options for advances configurations 

| **Argument**              | **Description**                                                                     |
|---------------------------|-------------------------------------------------------------------------------------|
| --umrlmr_lambda2          | Segmentation value for UMR/LMR identification (default 25)                          |
| --tol_val                 | Tolerance value for optimization (default 0.01)                                     |
| --pmd_lambda2             | Segmentation value for PMD identification (default 1000)                            |
| --pmd_std_threshold       | Standard deviation threshold for PMD (default 0.15)                                 |
| --umr_std_threshold       | Standard deviation threshold for UMR (default 0.1)                                  |
| --umr_max_beta            | Maximum methylation for UMR (default 0.1)                                           |
| --lmr_max_beta            | Maximum methylation for LMR (default 0.5)                                           |
| --valley_max_beta         | Maximum methylation for valleys (default 0.1)                                       |
| --pmd_valley_min_width    | Minimum width for PMD valleys (default 5000)                                        |
| --max_distance            | Maximum distance for merging segments (default 500)                                 |
| --min_width               | Minimum width for segments (default 30)                                             |
| --umr_large_width         | Large width threshold for UMRs (default 500)                                        |
| --umr_large_density       | Density threshold for large UMRs (default 0.03)                                     |
| --dmr_lambda2             | Segmentation value for DMR identification (default 25)                              |
| --min.overlap.pc          | Minimum overlap percentage for DMRs (default 10)                                    |
| --flank.width             | Minimum width of flanking segments (default 300)                                    |
| --flank.dist.max          | Maximum distance of UMR/LMR to flanking segment (default 5000)                      |
| --flank.num.cpgs          | Consider as small UMR/LMRs those which have num.cpgs smaller than this threshold (default 10) |
| --flank.beta.min          | For small UMR/LMRs, require that flanking segments have mean methylation above this threshold (default 0.5) |
| --split.pmds              | Whether to split PMDs which contain UMRs (default TRUE)                             |
| --merge.pmds              | Whether to merge adjacent PMDs into a single region (default TRUE)                  |
| --drop                    | Whether to drop verbose columns and unannotated segments (default TRUE)             |


## 7. &emsp;Output files and plots  

Two histograms are generated for each replicate in each condition in a pdf file name_replicate_plots.pdf containing: a distribution of the sequencing depth and distribution of the methylation of CpGs obtained after applying the minimum coverage per CpG (option -m).

For the LMR/UMR/DMV/PMD part, two tab-delimited files will be generated for each condition name_lmr_umr_dmv.tsv and name_pmd.tsv containing:
chr / start / end / number of CpGs / methylation / standard deviation / category

After PMDs identification, a smooth scatter plot containing the standart deviation and mean methylation of each PMDs is generated in a pdf file for each condition name1_pmds_plot.pdf and name2_pmds_plot.pdf.

For the DMR part, one tab-delimited file will be generated name1_vs_name2_dmrs.tsv containing:
chr / start / end / number of CpGs condition 1 / number of CpGs condition 2 / replicate coverage score / mean methylation 1 / mean methylation 2 / methylation difference / p-value / q-value / annotation
If the first part is not run (option -s), the last annotation column will not be added.

After DMRs identification, a smooth scatter plot containing the mean methylation 1 and mean methylation 2 of each DMRs is generated in a pdf file name1_vs_name2_dmrs_plot.pdf.

    
## 8. &emsp;Example

An example of command line to run MethyLasso is provided using files from Hansen et al. 2011 comparing healthy and colon cancer conditions with 3 replicates for each ([healthy_1](https://g-948214.d2cf88.03c0.data.globus.org/input_files/WGBS_colon-primary-normal_1_meth.bed.gz), [healthy_2](https://g-948214.d2cf88.03c0.data.globus.org/input_files/WGBS_colon-primary-normal_2_meth.bed.gz), [healthy_3](https://g-948214.d2cf88.03c0.data.globus.org/input_files/WGBS_colon-primary-normal_3_meth.bed.gz), [cancer_1](https://g-948214.d2cf88.03c0.data.globus.org/input_files/WGBS_colon-primary-tumor_1_meth.bed.gz), [cancer_2](https://g-948214.d2cf88.03c0.data.globus.org/input_files/WGBS_colon-primary-tumor_2_meth.bed.gz), [cancer_3](https://g-948214.d2cf88.03c0.data.globus.org/input_files/WGBS_colon-primary-tumor_3_meth.bed.gz)). The input files contain for each CpG , their genomic coordinates in the first three columns (chr / start / end), the total Count of Cs in columns four (--cov 4) and the percentage of methylation in columns 5 (--meth 5). DMRs will be identified with a q-value threshold 0.05 (-q 0.05).

Rscript MethyLasso.R --n1 healthy --c1 WGBS_healthy_1_meth.bed.gz,WGBS_healthy_2_meth.bed.gz,WGBS_healthy_3_meth.bed.gz --n2 cancer --c2 WGBS_cancer_1_meth.bed.gz,WGBS_cancer_2_meth.bed.gz,WGBS_cancer_3_meth.bed.gz --cov 4 --meth 5 -q 0.05

Output files obtained by running this command are LMRs/UMRs/DMVs file [healthy](https://g-948214.d2cf88.03c0.data.globus.org/output_files/normal_lmr_umr_dmv.tsv) and [cancer](https://g-948214.d2cf88.03c0.data.globus.org/output_files/cancer_lmr_umr_dmv.tsv) , PMDs file for [healthy](https://g-948214.d2cf88.03c0.data.globus.org/output_files/normal_pmd.tsv) and [cancer](https://g-948214.d2cf88.03c0.data.globus.org/output_files/cancer_pmd.tsv) and DMRs file [cancer vs healthy](https://g-948214.d2cf88.03c0.data.globus.org/output_files/cancer_vs_normal_dmrs.tsv). 


Also plots containing information on each input files are generated for healthy condition:  
<img src="https://github.com/abardet/methylasso/blob/0040e082382a6c32eb51fc7f3d557d57aa491595/output_plots/healthy_1_plots.png" alt="healthy1" width=30% height=40%> <img src="https://github.com/abardet/methylasso/blob/0040e082382a6c32eb51fc7f3d557d57aa491595/output_plots/healthy_2_plots.png" alt="healthy2" width=30% height=40%> <img src="https://github.com/abardet/methylasso/blob/0040e082382a6c32eb51fc7f3d557d57aa491595/output_plots/healthy_3_plots.png" alt="healthy3" width=30% height=40%>

and for cancer condition:   
<img src="https://github.com/abardet/methylasso/blob/013e3fae5dc76578d158a269d1429f89da65bc57/output_plots/cancer_1_plots.png" alt="cancer1" width=30% height=40%> <img src="https://github.com/abardet/methylasso/blob/013e3fae5dc76578d158a269d1429f89da65bc57/output_plots/cancer_2_plots.png" alt="cancer2" width=30% height=40%> <img src="https://github.com/abardet/methylasso/blob/013e3fae5dc76578d158a269d1429f89da65bc57/output_plots/cancer_3_plots.png" alt="cancer3" width=30% height=40%>

A plot on PMDs for each condition, where the user can see what threshold to take if he wants to recover more methylated PMDs and change the -m option.   
<img src="https://github.com/abardet/methylasso/blob/0040e082382a6c32eb51fc7f3d557d57aa491595/output_plots/healthy_pmds_plot.png" alt="pmd_healthy" width=40% height=40%> <img src="https://github.com/abardet/methylasso/blob/0040e082382a6c32eb51fc7f3d557d57aa491595/output_plots/cancer_pmds_plot.png" alt="pmd_cancer" width=40% height=40%>

And a plot on DMRs methylation in cancer and healthy conditions.   
<img src="https://github.com/abardet/methylasso/blob/0040e082382a6c32eb51fc7f3d557d57aa491595/output_plots/cancer_vs_healthy_dmrs_plot.png" alt="dmrs_cancer_vs_healthy" width=40% height=40%>


