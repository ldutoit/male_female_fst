# Re analyse flycatcher data


## Creating a VCF of CDS only

First, I get CDS only based on code [here](https://github.com/ldutoit/SBE/blob/master/createCDSonly.py (permisson needed for access)) and the flycatcher annnotation [Ficedula_albicollis.fAlb15.e73.gtf]( Ficedula_albicollis.fAlb15.e73.gtf). I reproduce what has been done for Dutoit et al. 2018.


```python
#module load Python/2.7.16-gimkl-2018b BEDTools
print "starting"
import windows_tools as wt  # third party module from https://github.com/ldutoit/personal_libraries
import os
os.mkdir("beds")
os.mkdir("beds/CDSbeds")
os.chdir("beds") 
!cut -f 9 ../Ficedula_albicollis.fAlb15.e73.gtf | sed -E 's/gene_id\s\"(ENSFALG[0-9]+).*/\1/g' > gene_names.bed
os.system("cut -f 1,3,4,5 ../Ficedula_albicollis.fAlb15.e73.gtf  | paste - gene_names.bed  > Ficedula_albicollis.fAlb15.e73_clean.gtf") # not really a gtf anymore
#avoid having starts before end
output = open("Ficedula_albicollis.fAlb15.e73_cleanswapped.gtf","w")
with open("Ficedula_albicollis.fAlb15.e73_clean.gtf") as f:	
	for line in f:
		start = str(min([int(x) for x in line.split("\t")[2:4]]))
		end = str(max([int(x) for x in line.split("\t")[2:4]]))
		output.write("\t".join(line.split()[:2]+[start]+[end]+line.split()[4:])+"\n")
output.close()

os.system("awk '{print $1,$3-1,$4,$5,$2}' Ficedula_albicollis.fAlb15.e73_cleanswapped.gtf >  Ficedula_albicollis.fAlb15.e73_clean.bed") # gtf to bed
os.system("sed -i 's/\s/\t/g' Ficedula_albicollis.fAlb15.e73_clean.bed") # remove spaces for tab
os.system("sort -k 4  -o  Ficedula_albicollis.fAlb15.e73_clean_sorted.bed Ficedula_albicollis.fAlb15.e73_clean.bed") # sort by gene names
os.system("rm Ficedula_albicollis.fAlb15.e73_cleanswapped.gtf Ficedula_albicollis.fAlb15.e73_clean.gtf Ficedula_albicollis.fAlb15.e73_clean.bed  Ficedula_albicollis.fAlb15.e73.gtf  gene_names.bed")	##CAREFUL THERE IS OVERLAP#
os.chdir("..") 

Annotated_file = "beds/Ficedula_albicollis.fAlb15.e73_clean_sorted.bed"


#get a CDS only file
os.system(" grep CDS "+Annotated_file+" > beds/Ficedula_albicollis.fAlb15.e73_clean_sortedCDSonly.bed")

## create a dictionnary of genes as keys with CDS as values
annot = wt.Bed("beds/Ficedula_albicollis.fAlb15.e73_clean_sortedCDSonly.bed", 0, 1, 2).windows

dict_genes = {}
ngene=0
with open(Annotated_file) as f:
	for line in f:
		ngene+=1
		print ngene
		if "CDS" in line:
			geneID = line.split()[3]
			if geneID in dict_genes:
				dict_genes[geneID].append(line)
			else:
				dict_genes[geneID]= [line]
ngene=0
for gene in dict_genes.keys():
	ngene+=1
	print ngene
	output=open("beds/CDSbeds/"+gene+".bed","w")	
	output.write("".join(dict_genes[gene]))
	output.close()	#

ngene
#Out[11]: 15290


##merge CDS if needed 
i = 0
for bed in os.listdir("beds/CDSbeds/"):
	i+=1
	print i,"beds/CDSbeds/"+bed
	os.system("sortBed -i beds/CDSbeds/" + bed + " | mergeBed -i - -d -1 > temp")
	output=open("beds/CDSbeds/" + bed,"w")
	with open("temp") as f:
		for line in f:
			output.write(line.strip()+"\tCDS\n")#
	output.close()
	os.system("rm temp")

###remove genes over many scaffolds 

gene="ENSFALG00000000001"
manyscaf = []
i=0
with open("beds/Ficedula_albicollis.fAlb15.e73_clean_sorted.bed") as f:
	wingenes= []
	for line in f:
		if "CDS" in line and gene==line.split()[3]:
			wingenes.append(line.split()[0])
			#print wingenes
		else:
			if len(set(wingenes))>1:manyscaf.append(line.split()[3])
			i+=1
			print i
			wingenes= []
			gene= line.split()[3]


manyscaf

#"['ENSFALG00000000003',
#" 'ENSFALG00000000807',
#" 'ENSFALG00000005247',
#" 'ENSFALG00000009191',
#" 'ENSFALG00000010625',
#" 'ENSFALG00000010684',
#" 'ENSFALG00000011070',
#" 'ENSFALG00000012157']

for gene in manyscaf:
	if os.path.exists("beds/CDSbeds/"+gene+".bed"):
		os.remove("beds/CDSbeds/"+gene+".bed")

####remnove all genes overlapping each other
# One I create a tempfile with all windows

output = open("temp","w")
i=0
for filename in os.listdir("beds/CDSbeds/"):
	i+=1
	print i
	with open("beds/CDSbeds/"+filename) as f:
		for line in f:
			output.write(line.strip()+"\t"+filename.split(".bed")[0]+"\n")
output.close()

#I load it in memory and split it in lists by scaffolds within a dictionnary
cdswindows  = wt.Bed("temp",0,1,2)
window_dict= {}
for win in cdswindows.windows:
	if win.seq in window_dict.keys():
		window_dict[win.seq].append(win)
	else:
		window_dict[win.seq] = [win]
 
#I store all the overlaps
overlaps = []
i=0
for item in cdswindows.windows:
	i+=1
	print i
	a = item.overlap_listwindows(window_dict[item.seq],out="list")
	if len(a)>1: 
		for overlap in a:
			if overlap.allcols[1] != item.allcols[1]:
				overlaps.append([item,overlap])

#I look at how many genes are involved
geneset = set()
for ov in overlaps:
 	geneset.add(ov[0].allcols[1])
 	geneset.add(ov[1].allcols[1])

i=0
for gene in geneset:
	if not "upward" in gene:
		if os.path.exists("beds/CDSbeds/"+gene+".bed"):
			i+=1
			os.remove("beds/CDSbeds/"+gene+".bed")
print i

len(os.listdir("beds/CDSbeds/"))
```

2 more genes that Dutoit at al. 2018, not sure which ones but no issue out of 14915.

## Creating males and females files
What is happening now is specific to this project. I will go grab all the SNPs for each of the positions.


```python
#module load SAMtools

############PSEUDOCOCDE
#for gene in the cds ..
#	open the file :
#		tabix each region to a new vcf#

#add headers
#sort the new vcf
#remove duplicates

import os
i=0
os.system("rm allcds_toclean.vcf")
for filename in os.listdir("beds/CDSbeds/"):
	i+=1
	print i, filename
	with open("beds/CDSbeds/"+filename) as f:
		for line in f:
			scaf,start,end = line.split()[:3]
			os.system("tabix GotlandsSNps95.vcf.gz "+scaf+":"+str(int(start)+1)+"-"+end+" >> allcds_toclean.vcf")# the plus one is bed 0 based to vcf 1 based

add headers, sortit

os.system(`zcat GotlandsSNps95.vcf.gz | head -n 100000 | grep "^#" | cat - allcds_toclean.vcf | bcftools sort - | uniq > GotlandsSNps95_onlycds_clean.vcf`)
```

*The rest is focusing on selecting the right individuals*

I keep only the right males and females ( infrp from dutoit et al suppinfo)

```python
allmales =  [ "15M153", "15M155", "15M158", "15M160", "15M161", "15M162", "15M163", "15M201", "15M202", "15M203", "15M204", "15M207", "15M468", "15M469", "15M475", "15M477", "15M49", "15M537", "15M568", "15M571", "15M573", "15M589", "15M684", "15M724", "93M25", "93M27", "93M28", "93M29", "93M36", "93M38", "93M39", "93M40", "93M41", "93M46", "93M53", "93M55", "93M58", "93M71", "93M72", "93M73", "93M78", "93M79", "93M80", "93M81", "93M83", "93M84", "93M86"]
allfemales = [ "15F129", "15F130", "15F131", "15F135", "15F142", "15F143", "15F145", "15F149", "15F151", "15F17", "15F18", "15F21", "15F22", "15F23", "15F24", "15F25", "15F29", "15F447", "15F448", "15F450", "15F453", "15F457", "15F459", "15F460", "93F24", "93F26", "93F30", "93F32", "93F34", "93F35", "93F42", "93F44", "93F45", "93F47", "93F54", "93F56", "93F59", "93F74", "93F75", "93F77", "93F82", "93F88", "93F89", "93F90", "93F92", "93F93", "93F94"]
	
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

keep only those individuals and sites without msising data, calculate allele frequencies and then small text processing to get it into one file


```bash
module load VCFtools
vcftools --vcf GotlandsSNps95_onlycds_clean.vcf --keep indstokeep.txt --max-missing-count 0 --recode --min-alleles 2 --max-alleles 2

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


paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt >  clean_frequencies_flycatcher.txt
## replace manually the first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq

```


That last file [freq_files/clean_frequencies_flycatcher.txt](freq_files/clean_frequencies_flycatcher.txt) is the same I'll obtain for every dataset

```
###summary info:

137202 sites in coding sequences of non-overlapping genes
all site are sequenced for every male (N=47 > 94 alleles) and every female (N=47 > 94 females).
```


### CDS MAF

above  0.05 (A lot more above 0.05 than humans, highlighting the different SFS)
```bash
module load VCFtools
vcftools --vcf GotlandsSNps95_onlycds_clean.vcf --keep indstokeep.txt --max-missing-count 0 --recode --min-alleles 2 --max-alleles 2 --maf 0.05
#After filtering, kept 95974 out of a possible 162675 Sites
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


paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt >  clean_frequencies_flycatcherMAFabove005.txt
## replace manually the first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq

```



## Generate 1 bootstrap flycatcher

I create 100 samples where I assign males and females randomly but keep number of males and females equal.

```bash
module load VCFtools
vcftools --vcf GotlandsSNps95_onlycds_clean.vcf --keep indstokeep.txt --max-missing-count 0 --recode --min-alleles 2 --max-alleles 2
echo -e "scaf\tpos\tn_males_allele_covered\tmale_freq\tn_females_allele_covered\tfemale_freq" > header.txt

``` 

The below code does the bootstrapping from the vcf (note it is a bit gly as it is essentially bash code wrapped into python)
```python
import os
os.mkdir("100bootstraps")
nmales=47
nfemales=47
nsites= 137202

print("Bootstrap",i)
##randomise the indtokeepfile
os.system(" shuf indstokeep.txt > randomisedindstokeep.txt")
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
os.system("paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt | tail -n " +str(nsites) + " | cat header.txt - >  100bootstraps/boot_flycatcher.txt")
## first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq
```

