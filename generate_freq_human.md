# Generate freq files human

## generate a list of males and females
[igsr_samples.tsv](igsr_samples.tsv) is the metadata including the sex. We use to create males and females list files in python.

the individuals we have:

```python
dict_sex = {}
#metadata
with open("igsr_samples.tsv") as f :
	for line in f:
		dict_sex[line.split()[0]] =line.split()[1] 

#our individuals
with open("1000humans_onlycds_clean.vcf") as f:
	for line in f:
		if line.startswith("#C"): break
inds = line.split()[9:]
allmales = [ind for ind in inds if dict_sex[ind]=="male"]
allfemales = [ind for ind in inds if dict_sex[ind]=="female"]

len(allmales) + len (allfemales) == len(inds)
```


```python
import os
	
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
We now have males.txt and females.txt.

We can now filter the vcf to keep only those individuals and sites without msising data, calculate allele frequencies and then small text processing to get it into one frquency file.


### frequency file without MAF filtering


```bash
module load VCFtools
vcftools --vcf 1000humans_onlycds_clean.vcf --keep indstokeep.txt --max-missing-count 0 --recode --min-alleles 2 --max-alleles 2
#After filtering, kept 121487 out of a possible 121866 Sites'
#separate males and females in two vcf to get allele frequencies.
vcftools --vcf out.recode.vcf --keep males.txt --freq
mv out.frq males_frequency.frq
vcftools --vcf out.recode.vcf --keep females.txt --freq
mv out.frq females_frequency.frq

### combine file
cut -f 5 males_frequency.frq | cut -f 2  -d ":" > tempmale_p.txt
cut -f 4  males_frequency.frq > male_allele_N.txt

cut -f 5 females_frequency.frq | cut -f 2  -d ":" > tempfemale_p.txt
cut -f 1-2 males_frequency.frq > temp_posinfo.txt
cut -f 4  females_frequency.frq > female_allele_N.txt


paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt >  clean_frequencies_humansCDS.txt
## replace manually the first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq

```

We realised at this point that a few positions were invariable, we remove them here in R. This is not a problem for any other dataset as they are fildered for maf.

```r
data<-read.table("clean_frequencies_humansCDS.txt",h=T)
dim(data[which(data$male_freq==data$female_freq),]) # 333 positions and none of them are at p=0.5, all p= 1 or 0, which would generate Na fst.
data<-data[-which(data$male_freq==data$female_freq),]
write.table(data,"clean_frequencies_humansCDS.txt",sep="\t",quote=F,row.names=F)


That last file [freq_files/clean_frequencies_humansCDS.txt](freq_files/clean_frequencies_humansCDS.txt) is the same I'll obtain for every dataset

###summary info:

121148 sites in coding sequences of non-overlapping genes
all site are sequenced for every male (2466 alleles) and every female (2542 females).
```
### Permutation no maf filtering

The below code does the bootstrapping from the vcf (note it is a bit ugly as it is essentially bash code wrapped into python)

```python
import os
nmales=1233
nfemales=1271
nsites= 121148
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

### combine file
os.system("cut -f 5 males_frequency.frq | cut -f 2  -d ':' > tempmale_p.txt")
os.system("cut -f 4  males_frequency.frq > male_allele_N.txt")
os.system("cut -f 5 females_frequency.frq | cut -f 2  -d ':' > tempfemale_p.txt")
os.system("cut -f 1-2 males_frequency.frq > temp_posinfo.txt")
os.system("cut -f 4  females_frequency.frq > female_allele_N.txt")
#paste everything together and replace the header as below too
os.system("paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt | tail -n " +str(nsites) + " | cat header.txt - >  permutation_nomaffilter.txt")
## first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq
```
### CDS MAF>0.05

above  0.05 . Relatively few snps but cannot really go at much lower threshold because I don't have that many flycatchers and pipefish.

```bash
module load VCFtools
vcftools --vcf 1000humans_onlycds_clean.vcf --keep indstokeep.txt --max-missing-count 0 --recode --min-alleles 2 --max-alleles 2 --maf 0.05
#After filtering, kept 7477 out of a possible 121866 Sites
#separate males and fgemales in two vcf to get frequency fast!!!!!
vcftools --vcf out.recode.vcf --keep males.txt --freq
mv out.frq males_frequency.frq
vcftools --vcf out.recode.vcf --keep females.txt --freq
mv out.frq females_frequency.frq

### combine file
cut -f 5 males_frequency.frq | cut -f 2  -d ":" > tempmale_p.txt
cut -f 4  males_frequency.frq > male_allele_N.txt

cut -f 5 females_frequency.frq | cut -f 2  -d ":" > tempfemale_p.txt
cut -f 1-2 males_frequency.frq > temp_posinfo.txt
cut -f 4  females_frequency.frq > female_allele_N.txt


paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt >  clean_frequencies_humansCDSMAFabove005.txt
## replace manually the first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq

```



### CDS MAF>0.05 permuted sex

above  0.05 . Relatively few snps but cannot really go at much lower threshold because I don't have that many flycatchers and pipefish.

```bash
module load VCFtools
vcftools --vcf 1000humans_onlycds_clean.vcf --keep indstokeep.txt --max-missing-count 0 --recode --min-alleles 2 --max-alleles 2 --maf 0.05
```

The below code does the bootstrapping from the vcf (note it is a bit ugly as it is essentially bash code wrapped into python)

```python
import os
#os.mkdir("100bootstraps")
nmales=1233
nfemales=1271
nsites= 7477
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

### combine file
os.system("cut -f 5 males_frequency.frq | cut -f 2  -d ':' > tempmale_p.txt")
os.system("cut -f 4  males_frequency.frq > male_allele_N.txt")
os.system("cut -f 5 females_frequency.frq | cut -f 2  -d ':' > tempfemale_p.txt")
os.system("cut -f 1-2 males_frequency.frq > temp_posinfo.txt")
os.system("cut -f 4  females_frequency.frq > female_allele_N.txt")
#paste everything together and replace the header as below too
os.system("paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt | tail -n " +str(nsites) + " | cat header.txt - >  permutation_mafabove005_humans.txt")
## first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq
```
