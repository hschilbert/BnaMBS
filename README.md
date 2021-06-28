# BnaMBS
scripts belonging to the GSL, SPC, and SOC mapping-by-sequencing analyses.


### Please get in touch if you need help running the scripts on your own data sets: [Hanna Schilbert (email)](mailto:hschilbe@cebitec.uni-bielefeld.de?subject=[GitHub]BnaMBS_scripts_request) ###

## Usage

### General recommendation

Full paths should be used to specify input and output files and folders. Sequence names should not contain white space characters like spaces and TABs. Underscores can be used to replace spaces. All files should be provided in tab separated format.

### 1 Gold standard

The scripts belonging to the generation of the gold standard are marked with a 1., they should be applied in the listed order.

## 1.1 filter_parent_variants.py
This script filters for high quality variants by applying coverage filters. 

```
Usage:
  python filter_parent_variants.py --vcf <FILE> --cov <FILE> --out <FILE> --min_cov_cov <INT> --max_cov_cov <INT> --min_cov_vcf <INT> --max_cov_vcf <INT> --black_vcf <FILE>

  Mandatory:
  
  Inputs 
  --vcf           STR   vcf file containing variants derived from the mapping of parental genotype 1 against B. napus Darmor-bzh
  --cov           STR   coverage file containing read coverage information per position derived from the mapping of parental genotype 2 against B. napus Darmor-bzh
          
  Output file
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


## 1.2 combine_homo_VCFs_vs_Bn41.py
This script filters for homozygous variants which are unique per parental genotype. Triallelic variants and variants present in both parents were excluded from further analyses as these are not contrasting between the pools.

```
Usage:
  python combine_homo_VCFs_vs_Bn41.py --P1_VCF <FILE> --P2_VCF <FILE> --out <DIR> 

  Mandatory:
  
  Inputs 
  --P1_VCF        STR   vcf file containing the variants filtered for coverage of parental genotype 1 (see 1.1)
  --P2_VCF        STR   vcf file containing the variants filtered for coverage of parental genotype 2 (see 1.1)
          
  Output directory
  --out           STR   output directory
```

## 1.3 filter_vcf_F1.py
This script uses the set of homozygous SNVs in the parents and further screens the set for heterozygosity in a reconstituted F1 population. The reconstituted F1 variant set comprises variants derived from all analysed genomic sequencing data of our study. Heterozygous variants are defined as having an allele frequency between 0.2 0.8 against the B. napus reference genome sequence.

```
Usage:
  python combine_homo_VCFs_vs_Bn41.py --F1_VCF <FILE> --P1_VCF <FILE> --P2_VCF <FILE> --out <DIR> 

  Mandatory:
  
  Inputs 
  --F1_VCF        STR   vcf file containing reconstituted F1 variant set
  --P1_VCF        STR   vcf file containing the variants filtered for homozygosity of parental genotype 1 (see 1.2)
  --P2_VCF        STR   vcf file containing the variants filtered for homozygosity of parental genotype 2 (see 1.2)
          
  Output directory
  --out           STR   output directory
```

## 1.4 merge_vcfs.py
This script combines the homozygous SNVs of the parents (1.2) with the variants identified to be heterozygous in step 1.3 to generate the final gold standard. 

```
Usage:
  python merge_vcfs.py --in <FILE> --fasta <FILE> --out <FILE> --sort_script <FILE>

  Mandatory:
  
  Inputs  
  --in            STR   path to input folder where vcfs to be combined are located
  --fasta         STR   path to reference genome sequence fasta file
          
  Output file
  --out           STR   output file
  
  Optional:
  --sort_script   STR   path to sort script sort_vcf_by_fasta.py 
  
```

`--sort_script` full path to sort_vcf_by_fasta.py can be provided. If this parameter is not used sort_vcf_by_fasta.py should be located in the working directory.

## 1.4.1 sort_vcf_by_fasta.py
This script sorts a vcf file based on a reference genome sequence. 

```
Usage:
  python sort_vcf_by_fasta.py --vcf <FILE> --fasta <FILE> --output <FILE>

  Mandatory:
  
  Inputs  
  --vcf           STR   vcf file to be sorted by chromosomal position
  --fasta         STR   path to reference genome sequence fasta file
          
  Output file
  --out           STR   output file, which is a sorted vcf file
```
