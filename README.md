# Male - Female Fst re-analysis

In this folder, we compare male female Fst of published datasets on a SNP per SNP basis to a null expected distribution. This README aims at guiding the interested reader quickly to key functions and data of the project. We would be happy to provide further help with this folder, do not hesitate to contact dutoit.ludovic@gmail.com.

This folder does not contain the raw or intermediate data files but it contains all the code from the raw data as well as the vectors used to create the figures in the file [vectors_reanalysis.RData](vectors_reanalysis.RData). See below for more detailed information. 

## Datasets

This re-analyis include three datasets	:


* 1. [Dutoit et al. 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14789) using 95,974 sites sequenced for 47 males and 47 females. All the sites are in coding sequences (no missing data).Minor allele frequency(i.e. MAF)>0.05

* 2. [Flanagan and Jones 2017](https://onlinelibrary.wiley.com/doi/full/10.1111/evo.13173) using 44,773 sites RAD-Seq generated SNPs, 57 females and 167 males. Up to 50% missing data. Number of sequenced alleles saved for each sex at each site. MAF>0.05

* 3. 1000 Human genomes. 1233 males and 1271 females. 7,477 sites. All the sites are in coding sequences. MAF>0.05

* 4. 1000 human genomes, but no MAF filtering. 1233 males and 1271 females. ll the sites are in coding sequences.

## Analysis 

### Generate allele-frequency files

Several different allele frequencies files were generated for each species. [generate_freq_pipefish.md](generate_freq_pipefish.md) and [generate_freq_flycatcher.md](generate_freq_flycatcher.md) were used to generate frequency files for pipefish and flycatcher. [downloadhuman.md](downloadhuman.md) was followed by [generate_freq_human.md](generate_freq_human.md) to generate frequency files for humans. The first part of those files is just generating clean VCF files with only the sites of interest as described above. While those files are a clean record of what has been done from the raw data to the alleles frequency files, they might be a bit long and tedious and the reader might be interested in the summary below on how to get from a vcf to alleles_frequency (permuted/non-permuated):

```bash
#!/bin/sh
##Parameters
NMALES=47 # 47 in the flycatcher example
NFEMALES=47 # 47 in the flycatcher example

#randomise individuals starting from a list of male codes and females codes in males.txt and females.txt :
cat males.txt females.txt > indstokeep.txt
shuf indstokeep.txt > randomisedindstokeep.txt 

#Split reandomised males and females list files
head -n $NFEMALES randomisedindstokeep.txt > females_randomised.txt
tail -n $NMALES randomisedindstokeep.txt > males_randomised.txt

##Allele frequency calculation. For real data, just use males and females.txt instead of the randomised above created files
vcftools --vcf out.recode.vcf --keep males_randomised --freq
mv out.frq males_frequency.frq
vcftools --vcf out.recode.vcf --keep females_randomised.txt --freq
mv out.frq females_frequency.frq

### Paste males and females frequency to create one allele_frequency file
cut -f 5 males_frequency.frq | cut -f 2  -d ':' > tempmale_p.txt #Male allele_freq
cut -f 4  males_frequency.frq > male_allele_N.txt #Male allele number
cut -f 5 females_frequency.frq | cut -f 2  -d ':' > tempfemale_p.txt # Female allele_freq
cut -f 1-2 males_frequency.frq > temp_posinfo.txt #Position information
cut -f 4  females_frequency.frq > female_allele_N.txt  # Female alle;e number
#paste everything together and replace the header as below too
paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt | tail -n " +str(nsites) + " | cat header.txt - >  permuted_allele_frequency.txt
#The header can then be replaced by:
scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq
```

### Creating_vectors

From the allele-frequency files, all the useful functions are stored in [Creating_vectors.md](Creating_Vectors.md), that also create all the data useful for this article in the R object.

The functions are summarised below:

```r
nei_fst_func(pm,pf) # calculatte fst for a snp or a vector of snps each SNP based on inputted allele_frequencies
#Parameter of the chi-square distribution given nm and nf 
exp_chisq(nm,nf) # calculate the parameter of the chi-square distribution to create the null theoretical distribution based on inputted number of ALLES for males and females
quantile_from_expdistrib(pm,pf,nm,nf) # Calculate the quantile of the expected null distribution in which a SNP fall, given nm and nf and the observed allele frequencies. Useful to compare SNP with different number of males and females
howmuchabove(pm,pf,nm,nf,percent) # Compute the percentage of data above a given quantile (i.e)

#theoretical_null function  draws a theoretical null by taking one point from the chi-square null for every SNP, allow to deal with missing data when different SNPs have different Null because different numbers of alleles, it requires draw_one
draw_one(x)
theoretical_null(nm,nf)
```


They then use these functions to generate all the data necessary for the figures, in the one object  [vectors_reanalysis.RData](vectors_reanalysis.RData)

it can be loaded in R as follow:

```
load("vectors_reanalysis.RData")
README()
```
**Overall FST and chi-square**

We also compared overall Fst to 1000 permutation per dataset created using [1000permutations.py](1000permutations.py). That was done in the script: [overallFST.md](overallFST.md)

We compared the tail of the observed distribution to the tail of the expected distributions and used chi-square to test for tail enrichment at 5% and 1% in [chisquare.md](chisquare.md) 
### Figures

The figure are produced in the file [plots.R](plots.R)



### Raw data storage

While all the data is published, a copy of all the raw files is stored on the High Capacity Storage at Otago University under: 

```
/Volumes/sci-bioinformatics-project-archive/sexantagonismfst
```

Contact dutoit.ludovic@gmail.com  to ask for access.
