## generate_freq_pipefish

first I create files with only males and females 

```python
#indstokeep
allmales = ["sample_NPM005_align","sample_NPM006_align","sample_NPM007_align","sample_NPM008_align","sample_NPM010_align","sample_NPM011_align","sample_NPM012_align","sample_NPM1128_align","sample_PRM001_align","sample_PRM002_align","sample_PRM003_align","sample_PRM005_align","sample_PRM006_align","sample_PRM007_align","sample_PRM009_align","sample_PRM010_align","sample_PRM011_align","sample_PRM012_align","sample_PRM013_align","sample_PRM014_align","sample_PRM015_align","sample_PRM016_align","sample_PRM017_align","sample_PRM018_align","sample_PRM019_align","sample_PRM022_align","sample_PRM023_align","sample_PRM025_align","sample_PRM026_align","sample_PRM028_align","sample_PRM041_align","sample_PRM042_align","sample_PRM043_align","sample_PRM044_align","sample_PRM045_align","sample_PRM046_align","sample_PRM047_align","sample_PRM048_align","sample_PRM049_align","sample_PRM050_align","sample_PRM051_align","sample_PRM053_align","sample_PRM054_align","sample_PRM055_align","sample_PRM056_align","sample_PRM057_align","sample_PRM058_align","sample_PRM059_align","sample_PRM060_align","sample_PRM061_align","sample_PRM062_align","sample_PRM064_align","sample_PRM065_align","sample_PRM066_align","sample_PRM067_align","sample_PRM068_align","sample_PRM069_align","sample_PRM070_align","sample_PRM071_align","sample_PRM072_align","sample_PRM073_align","sample_PRM074_align","sample_PRM076_align","sample_PRM077_align","sample_PRM078_align","sample_PRM079_align","sample_PRM080_align","sample_PRM081_align","sample_PRM082_align","sample_PRM083_align","sample_PRM084_align","sample_PRM085_align","sample_PRM086-23_align","sample_PRM086R_align","sample_PRM087_align","sample_PRM088_align","sample_PRM089_align","sample_PRM090_align","sample_PRM091_align","sample_PRM092_align","sample_PRM093_align","sample_PRM096_align","sample_PRM097_align","sample_PRM098_align","sample_PRM099_align","sample_PRM100_align","sample_PRM101_align","sample_PRM102_align","sample_PRM103_align","sample_PRM104_align","sample_PRM105_align","sample_PRM106_align","sample_PRM107_align","sample_PRM108_align","sample_PRM109_align","sample_PRM110_align","sample_PRM111_align","sample_PRM112_align","sample_PRM113_align","sample_PRM114_align","sample_PRM115_align","sample_PRM117_align","sample_PRM118_align","sample_PRM119_align","sample_PRM120_align","sample_PRM121_align","sample_PRM122_align","sample_PRM124_align","sample_PRM125_align","sample_PRM126_align","sample_PRM127_align","sample_PRM128_align","sample_PRM129_align","sample_PRM130_align","sample_PRM131_align","sample_PRM132_align","sample_PRM133_align","sample_PRM134_align","sample_PRM135_align","sample_PRM136_align","sample_PRM137_align","sample_PRM139_align","sample_PRM140_align","sample_PRM143_align","sample_PRM144_align","sample_PRM145_align","sample_PRM146_align","sample_PRM147_align","sample_PRM148_align","sample_PRM149_align","sample_PRM150_align","sample_PRM151_align","sample_PRM152_align","sample_PRM153_align","sample_PRM154_align","sample_PRM155_align","sample_PRM156_align","sample_PRM159_align","sample_PRM160_align","sample_PRM161_align","sample_PRM162_align","sample_PRM163_align","sample_PRM164_align","sample_PRM165_align","sample_PRM166_align","sample_PRM167_align","sample_PRM168_align","sample_PRM169_align","sample_PRM171_align","sample_PRM172_align","sample_PRM173_align","sample_PRM174_align","sample_PRM175_align","sample_PRM176_align","sample_PRM177-1_align","sample_PRM177_align","sample_PRM178_align","sample_PRM179_align","sample_PRM181_align","sample_PRM182_align","sample_PRM183_align","sample_PRM184_align","sample_PRM185_align","sample_PRM186_align","sample_PRM187_align","sample_PRM188_align","sample_PRM189_align"]

allfemales =  ["sample_FEM001_align","sample_FEM002_align","sample_FEM004_align","sample_FEM005_align","sample_FEM006_align","sample_FEM008_align","sample_FEM009_align","sample_FEM010_align","sample_FEM011_align","sample_FEM013_align","sample_FEM014_align","sample_FEM015_align","sample_FEM016_align","sample_FEM017_align","sample_FEM018_align","sample_FEM019_align","sample_FEM020_align","sample_FEM021_align","sample_FEM022_align","sample_FEM023_align","sample_FEM024_align","sample_FEM025_align","sample_FEM026_align","sample_FEM027_align","sample_FEM028_align","sample_FEM029_align","sample_FEM030_align","sample_FEM031_align","sample_FEM032_align","sample_FEM033_align","sample_FEM034_align","sample_FEM035_align","sample_FEM036_align","sample_FEM037_align","sample_FEM039_align","sample_FEM041_align","sample_FEM046_align","sample_FEM047_align","sample_FEM051_align","sample_FEM052_align","sample_FEM054-1_align","sample_FEM056_align","sample_FEM058_align","sample_FEM064_align","sample_FEM066_align","sample_FEM067_align","sample_FEM072_align","sample_FEM074_align","sample_FEM076_align","sample_FEM077_align","sample_FEM078_align","sample_FEM080_align","sample_FEM081_align","sample_FEM082_align","sample_FEM083_align","sample_FEM084_align","sample_FEM087_align"]

output =  open("indstokeep.txt", "w")

for ind in allmales+allfemales:
	output.write(ind+"\n")


output.close()
os.system("wc -l indstokeep.txt")

output =  open("males.txt", "w")
for ind in allmales :
	output.write(ind+"\n")

output.close()

output =  open("females.txt", "w")
for ind in allfemales :
	output.write(ind+"\n")

output.close()

```

keep only those individuals and sites with less than 50% data, calculate allele frequencies and then small text processing to get it into one file


```bash
module load VCFtools
vcftools --vcf batch_1.vcf --keep indstokeep.txt --max-missing-count 112 --recode --min-alleles 2 --max-alleles 2

#separate males and fgemales in two vcf to get frequency fast
vcftools --vcf out.recode.vcf --keep males.txt --freq
mv out.frq males_frequency.frq
vcftools --vcf out.recode.vcf --keep females.txt --freq
mv out.frq females_frequency.frq

### small test
cut -f 5 males_frequency.frq | cut -f 2  -d ":" > tempmale_p.txt
cut -f 4  males_frequency.frq > male_allele_N.txt

cut -f 5 females_frequency.frq | cut -f 2  -d ":" > tempfemale_p.txt
cut -f 1-2 males_frequency.frq > temp_posinfo.txt
cut -f 4  females_frequency.frq > female_allele_N.txt


paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt > clean_frequencies_pipefish.txt
## replace manually the first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq

```



```
###summary info:
45939 sites RAD SNPs, 57 females, up to 167 males, but missing data. Number of sequenced alleles saved for each sex at each site.
```

### MAF

```bash
 module load VCFtools
vcftools --vcf batch_1.vcf --keep indstokeep.txt --max-missing-count 112 --recode --min-alleles 2 --maf 0.05

#After filtering, kept After filtering, kept 44773 out of a possible 69109 Sites
#separate males and fgemales in two vcf to get frequency fast!!!!!
vcftools --vcf out.recode.vcf --keep males.txt --freq
mv out.frq males_frequency.frq
vcftools --vcf out.recode.vcf --keep females.txt --freq
mv out.frq females_frequency.frq

### small test
cut -f 5 males_frequency.frq | cut -f 2  -d ":" > tempmale_p.txt
cut -f 4  males_frequency.frq > male_allele_N.txt

cut -f 5 females_frequency.frq | cut -f 2  -d ":" > tempfemale_p.txt
cut -f 1-2 males_frequency.frq > temp_posinfo.txt
cut -f 4  females_frequency.frq > female_allele_N.txt


paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt >  clean_frequencies_pipefishMAFabove005.txt
## replace manually the first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq

```

below 0.05

```bash
module load VCFtools
vcftools --vcf batch_1.vcf --keep indstokeep.txt --max-missing-count 112 --recode --min-alleles 2 --max-maf 0.05
#After filtering, kept 1216 out of a possible 69109 Sites

#separate males and fgemales in two vcf to get frequency fast!!!!!
vcftools --vcf out.recode.vcf --keep males.txt --freq
mv out.frq males_frequency.frq
vcftools --vcf out.recode.vcf --keep females.txt --freq
mv out.frq females_frequency.frq

### small test
cut -f 5 males_frequency.frq | cut -f 2  -d ":" > tempmale_p.txt
cut -f 4  males_frequency.frq > male_allele_N.txt

cut -f 5 females_frequency.frq | cut -f 2  -d ":" > tempfemale_p.txt
cut -f 1-2 males_frequency.frq > temp_posinfo.txt
cut -f 4  females_frequency.frq > female_allele_N.txt


paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt >  clean_frequencies_pipefishMAFbelow005.txt
## replace manually the first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq
```

Those MAF filtered file are in[freq_files/clean_frequencies_pipefishMAFbelow005.txt](freq_files/clean_frequencies_pipefishMAFbelow005.txt) and 
frequencies_pipefishMAFabove005.txt](freq_files/clean_frequencies_pipefishMAFabove005.txt) 


### Finally I create the global SFS
```bash
module load VCFtools
vcftools --vcf batch_1.vcf --keep indstokeep.txt --max-missing-count 112 --recode --min-alleles 2 --max-alleles 2
vcftools --vcf out.recode.vcf --keep indstokeep.txt --freq
cut -d ":" -f 3 out.frq  > alternate_allele_pipefish.freq
```


```r 
maf<-function(x){return(min(x,1-x))}
data<-read.table("alternate_allele_pipefish.freq")
freq<-as.numeric(as.character(data[,1]))
minor<-unlist(lapply(freq,maf))
pdf("hist_pipefish.pdf")
hist(minor,xlim=c(0,0.5),breaks=39,main="pipefish")
dev.off()
```

[freq_files/hist_pipefish.pdf](freq_files/hist_pipefish.pdf)


## Generate 100 bootstraps pipefish

In this case, I create 100 samples where I assign males and females randomly but keep number of males and females equal to what they are originally.

```bash
module load VCFtools
vcftools --vcf batch_1.vcf --keep indstokeep.txt --max-missing-count 112 --recode --min-alleles 2 --max-alleles 2
echo -e "scaf\tpos\tn_males_allele_covered\tmale_freq\tn_females_allele_covered\tfemale_freq" > header.txt

``` 

 The below code does the bootstrapping from the vcf (note it is a bit gly as it is essentially bash code wrapped into python)
```python
import os
#os.mkdir("100bootstraps")
nmales=167
nfemales=57
nsites= 45938
for i in range(1,101):
	print("Bootstrap",i)
	##randomise the indtokeepfile
	os.system("shuf indstokeep.txt > randomisedindstokeep.txt")
	#grab males and females randomly
	os.system("head -n "+str(nfemales)+" randomisedindstokeep.txt > females.txt")
	os.system("tail -n "+str(nmales)+" randomisedindstokeep.txt > males.txt")
	##go through the frequency calculation
	os.system("vcftools --vcf out.recode.vcf --keep males.txt --freq")
	os.system("mv out.frq males_frequency.frq")
	os.system("vcftools --vcf out.recode.vcf --keep females.txt --freq")
	os.system("mv out.frq females_frequency.frq")

	### small test
	os.system("cut -f 5 males_frequency.frq | cut -f 2  -d ':' > tempmale_p.txt")
	os.system("cut -f 4  males_frequency.frq > male_allele_N.txt")
	os.system("cut -f 5 females_frequency.frq | cut -f 2  -d ':' > tempfemale_p.txt")
	os.system("cut -f 1-2 males_frequency.frq > temp_posinfo.txt")
	os.system("cut -f 4  females_frequency.frq > female_allele_N.txt")
	#paste everything together and replace the header as below too
	os.system("paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt | tail -n " +str(nsites) + " | cat header.txt - >  100bootstraps/boot_"+str(i)+"pipefish.txt")
## first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq
```

###Permutation maf 005

```bash
 module load VCFtools
vcftools --vcf batch_1.vcf --keep indstokeep.txt --max-missing-count 112 --recode --min-alleles 2 --maf 0.05
```

The below code does the bootstrapping from the vcf (note it is a bit gly as it is essentially bash code wrapped into python)
```python
import os
#os.mkdir("100bootstraps")
nmales=167
nfemales=57
nsites= 44773
for i in range(1,2):
	print("Bootstrap",i)
	##randomise the indtokeepfile
	os.system("shuf indstokeep.txt > randomisedindstokeep.txt")
	#grab males and females randomly
	os.system("head -n "+str(nfemales)+" randomisedindstokeep.txt > females.txt")
	os.system("tail -n "+str(nmales)+" randomisedindstokeep.txt > males.txt")
	##go through the frequency calculation
	os.system("vcftools --vcf out.recode.vcf --keep males.txt --freq")
	os.system("mv out.frq males_frequency.frq")
	os.system("vcftools --vcf out.recode.vcf --keep females.txt --freq")
	os.system("mv out.frq females_frequency.frq")

	### small test
	os.system("cut -f 5 males_frequency.frq | cut -f 2  -d ':' > tempmale_p.txt")
	os.system("cut -f 4  males_frequency.frq > male_allele_N.txt")
	os.system("cut -f 5 females_frequency.frq | cut -f 2  -d ':' > tempfemale_p.txt")
	os.system("cut -f 1-2 males_frequency.frq > temp_posinfo.txt")
	os.system("cut -f 4  females_frequency.frq > female_allele_N.txt")
	#paste everything together and replace the header as below too
	os.system("paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt | tail -n " +str(nsites) + " | cat header.txt - >  permutation_mafabove005_"+str(i)+"pipefish.txt")
## first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq
```
