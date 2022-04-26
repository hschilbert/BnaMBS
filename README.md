# BnaMBS
scripts belonging to the GSL, SPC, and SOC mapping-by-sequencing analyses.


### Please get in touch if you need help running the scripts on your own data sets: [Hanna Schilbert (email)](mailto:hschilbe@cebitec.uni-bielefeld.de?subject=[GitHub]BnaMBS_scripts_request) ###

## Usage

### General recommendation

Full paths should be used to specify input and output files and folders. Sequence names should not contain white space characters like spaces and TABs. Underscores can be used to replace spaces. All files should be provided in tab separated format.

## 2.5 Generation of the gold standard 

The scripts belonging to the generation of the gold standard are marked with a 1., they should be applied in the listed order.

### 2.5.1 filter_parent_variants.py
This script filters for high quality variants by applying coverage filters. 

```
Usage:
  python filter_parent_variants.py --vcf <FILE> --cov <FILE> --out <FILE> --min_cov_cov <INT> --max_cov_cov <INT> --min_cov_vcf <INT> --max_cov_vcf <INT> --black_vcf <FILE>

  Mandatory:
  
  Inputs 
  --vcf           STR   vcf file containing variants derived from the mapping of parental genotype 1 against B. napus Darmor-bzh
  --cov           STR   coverage file containing read coverage information per position derived from the mapping of parental genotype 2 against B. napus Darmor-bzh
          
  Output 
  --out           STR   output file
  
  Optional:
  --min_cov_cov   INT   minimal read coverage at one position in the bam file [10]
  --max_cov_cov   INT   maximal read coverage at one position in the bam file [100]
  --min_cov_vcf   INT   minimal read coverage at one position in the vcf file [10]
  --max_cov_vcf   INT   maximal read coverage at one position in the vcf file [100]
  --black_vcf     STR   vcf file BOAS

```

`--vcf` vcf file containing variants derived from the mapping of parental genotype 1 against B. napus Darmor-bzh.

`--cov` coverage file containing read coverage information per position derived from the mapping of parental genotype 2 against B. napus Darmor-bzh.

`--out` specify the output file where the filtered variants will be stored.

`--min_cov_cov`, `--max_cov_cov`, `--min_cov_vcf`, `--max_cov_vcf` integers should be given for these parameters determining the filter threshold based on the files given in `--cov` and `--vcf`, respectively.

`--black_vcf` vcf file containing BOAS


### 2.5.2 combine_homo_VCFs_vs_Bn41.py
This script filters for homozygous variants which are unique per parental genotype. Triallelic variants and variants present in both parents were excluded from further analyses as these are not contrasting between the pools.

```
Usage:
  python combine_homo_VCFs_vs_Bn41.py --P1_VCF <FILE> --P2_VCF <FILE> --out <DIR> 

  Mandatory:
  
  Inputs 
  --P1_VCF        STR   vcf file containing the variants filtered for coverage of parental genotype 1 (see 1.1)
  --P2_VCF        STR   vcf file containing the variants filtered for coverage of parental genotype 2 (see 1.1)
          
  Output 
  --out           STR   output directory
```

### 2.5.3 filter_vcf_F1.py
This script uses the set of homozygous SNVs in the parents and further screens the set for heterozygosity in a reconstituted F1 population. The reconstituted F1 variant set comprises variants derived from all analysed genomic sequencing data of our study. Heterozygous variants are defined as having an allele frequency between 0.2 0.8 against the B. napus reference genome sequence.

```
Usage:
  python combine_homo_VCFs_vs_Bn41.py --F1_VCF <FILE> --P1_VCF <FILE> --P2_VCF <FILE> --out <DIR> 

  Mandatory:
  
  Inputs 
  --F1_VCF        STR   vcf file containing reconstituted F1 variant set
  --P1_VCF        STR   vcf file containing the variants filtered for homozygosity of parental genotype 1 (see 1.2)
  --P2_VCF        STR   vcf file containing the variants filtered for homozygosity of parental genotype 2 (see 1.2)
          
  Output 
  --out           STR   output directory
```

### 2.5.4 merge_vcfs.py
This script combines the homozygous SNVs of the parents (1.2) with the variants identified to be heterozygous in step 1.3 to generate the final gold standard. 

```
Usage:
  python merge_vcfs.py --in <FILE> --fasta <FILE> --out <FILE> --sort_script <FILE>

  Mandatory:
  
  Inputs  
  --in            STR   path to input folder where vcfs to be combined are located
  --fasta         STR   path to reference genome sequence fasta file
          
  Output 
  --out           STR   output file
  
  Optional:
  --sort_script   STR   path to sort script sort_vcf_by_fasta.py 
  
```

`--sort_script` full path to sort_vcf_by_fasta.py can be provided. If this parameter is not used sort_vcf_by_fasta.py should be located in the working directory.

### 2.5.5 sort_vcf_by_fasta.py
This script sorts a vcf file based on a reference genome sequence. 

```
Usage:
  python sort_vcf_by_fasta.py --vcf <FILE> --fasta <FILE> --output <FILE>

  Mandatory:
  
  Inputs  
  --vcf           STR   vcf file to be sorted by chromosomal position
  --fasta         STR   path to reference genome sequence fasta file
          
  Output 
  --out           STR   output file, which is a sorted vcf file
```

## 2.6 Filter raw variants per pool for delta allele frequency calculation

### 2.6.1 filter_pools_vcfs_for_gold_standard.py
This script filters the pool vcf files for the variants in the gold standard. 

```
Usage:
  python filter_pools_vcfs_for_gold_standard.py --out <DIR>
	
  Mandatory:
  
  Inputs
  Full paths to the to be filtered pool vcf files need to be specified in the script via their full paths and save in the variables vcf_hp (VCF high pool) and vcf_lp (VCF low pool). The vcf files containing the gold standard variants from the parents need to be specified via their full paths in the script and save in the variables vcf_J and vcf_L.
  
  Output file
  --out           STR   path to output folder
```

### 2.6.2 combine_single_VCFs.py
This script combines vcf files located in the same folder. 

```
Usage:
  python combine_single_VCFs.py --in <DIR> --out <FILE>

  Mandatory:
  
  Input  
  --in            STR   path to vcf file folder, which contains the vcf files which are going to be merged
          
  Output 
  --out           STR   path to output file, which is a merged vcf file
```			
					
## 2.7 Interval detection

### 2.7.1 fisher_exact_test_corrects_for_multiples_testing.py
This script identifies “statistically meaningful differential Allele specific Read Counts” (dARCs). For the interval detection Fisher’s exact test was applied on the raw SNVs of the pools to yield variants with a significant delta allele frequency. A p-value cut-off of 0.05 was applied after correction for multiple testing. The passing SNVs are called “statistically meaningful differential Allele specific Read Counts” (dARCs).

```
Usage:
  python fisher_exact_test_corrects_for_multiples_testing.py --in <FILE> --sig <INTEGER> --pool1 <COMMA_SEPARATED_LIST_OF_SAMPLES> --pool2 <COMMA_SEPARATED_LIST_OF_SAMPLES> --out <FILE>

  Mandatory:
  
  Input  
  --in          STR   path to vcf file
  --sig		INT   integer setting the significant level alpha [0.05] 
  --pool1       STR   list of samples, comma seperated
  --pool2       STR   list of samples, comma seperated
          
  Output 
  --out       	STR   path to output file
```
					
### 2.7.2 get_intervals_based_on_dARCs_Bn41_v4.py
This script identifies genomic intervals based on dARCs were used to identify genomic intervals associated with the analyzed traits 


### 2.7.3 PAV_finder.py 
This script identifies I) Zero coverage regions (ZCRs) by using the coverage information of both pools and applying a genome wide screening with a window size of 200 bp per chromosome. ZCRs are considered during the interval detection, as they are often responsible for the splitting of genomic intervals into parts. Where an interval is missing in both pools compared to the Darmor bzh reference genome sequence, no variants and hence no dARCs can be detected. 

This script identifies II) presence absence variants (PAVs). First, the average coverage per gene region per pool was calculated by calculating the mean of all coverage values per gene region. Next, a minimum coverage cut off was applied (-mincov 10) by ensuring that the sum of the average coverage for each gene region of both pools need to be greater than 10. This step was done to insure that the gene region is present in at least one pool. If no coverage was detected for one gene region of one pool, the coverage was set to 0.01. The average coverage was then normalized to the overall coverage of a sample by dividing through the median of all mean coverage values of a sample. Next, the log2 of the normalized coverage values of one gene region derived from both pools was calculated. Then high quality PAVs were extracted by filtering with an absolute value of log2(normalized coverage of gene region X of pool 1 / normalized coverage of gene region X of pool 2) > 1 and an absolute z score of > 1.5 and the normalized coverage value of at least one pool must be < 0.4 to insure a very low coverage alias absence of this gene region (github PAV_parser_genes.py)

```
Usage:
  python PAV_finder.py --cov1 <FILE> --cov2 <FILE> --out <DIR>

  Mandatory:
  
  Input  
  --in          STR   path to vcf file

  --cov1        STR   path to coverage file derived from one pool mapping 
  --cov2        STR   path to coverage file derived from the other pool mapping 
  
  Optional:
  --mode        STR   mode of PAV detection (gene|genomic|zcr)>[genomic]
  --gff         STR   gff annotation file, necessary if gene mode is chosen [none]
  --anno        STR   functional annotation file [none]
  --mincov      STR   minimal combined coverage of both samples per gene [-1]
  --blocksize   INT   size for genomic PAC or ZCR detection [3000]
  --maxrelcov   INT   relative coverage cutoff of ZCR detection [0.1]
	
  Output 
  --out       	STR   path to output directory
					
```

## 2.7 Generation of delta allele frequency plots

## 2.8 Generation of delta allele frequency plots

### 2.8.1 sophisticated_cov_plot.py
This script generates coverage and delta allele frequency plots. 

```
Usage:
  python sophisticated_cov_plot.py --input_vcf <FILE> --input_vcf_sig_SNP <FILE> --in_merged_ori_vcf <FILE> --reference_file <FILE> --high_pool <sample name in VCF; multiple samples names can be provided comma-seperated> --low_pool <sample name in VCF; multiple samples names can be provided comma-seperated> --output_dir <DIR> 

  Mandatory:
  
  Input  
  --input_vcf <FILE> 
  --input_vcf_sig_SNP <FILE> 
  --in_merged_ori_vcf <FILE> 
  --reference_file <FILE> 
  --high_pool <sample name in VCF; multiple samples names can be provided comma-seperated> 
  --low_pool <sample name in VCF; multiple samples names can be provided comma-seperated>   
  
  Output 
  --output_dir           STR   path to output directory [will be generated if required]
```

## 2.9 Functional annotation and candidate genes

### 2.9.1 fetch_gene_IDs_from_gff3_file.py
### 2.9.2 map_mean_exp_to_cand_genes_in_reg.py
### 2.9.3 map_PAVs_to_genes_in_regs.py

## 3.2 RNA-Seq

The scripts belonging to the RNA-Seq analysis are presented here.

### 3.2.1 parse_STAR_log_file_create_mapping_statistic.py

```
Usage:
  python parse_STAR_log_file_create_mapping_statistic.py --path_to_log_files <DIR> --out <FILE>

  Mandatory:
  
  Inputs  
  --path_to_log_files           STR   path to directory containing STAR log files

  Output 
  --out           STR   output file
```

### 3.2.2 generate_figures_only_mean_expression_calc.py

```
Usage:
  python generate_figures_only_mean_expression_calc.py --in <FILE> --genes <FILE> --samples <FILE> --out <DIR>

  Mandatory:
  
  Inputs  
  --in           STR   path to count table with samples as column names and gene IDs as row names
  --genes        STR   path to gene ID file, one column with one gene ID per row
  --samples      STR   path to sample ID file, one column with one sample ID per row. Sample ID must be present in count table
  
  Output 
  --out           STR   output directory
```

### 3.2.3 map_mean_exp_to_cand_genes_in_reg.py

```
Usage:
  python map_mean_exp_to_cand_genes_in_reg.py --mean_exp_table <FILE> --anno_vars_w_reg <FILE> --out <FILE>

  Mandatory:
  
  Inputs  
  --mean_exp_table        STR   path to count table with samples as column names and gene IDs as row names
  --anno_vars_w_reg       STR   path to file, containing annotated variants in genomic regions
  
  Output 
  --out           STR   output file
```

