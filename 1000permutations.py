#!/usr/bin/env python
#This file creates 1000 permutations for each dataset to compare overall Fst in the data tom Fst from 1000 bootstraped samples

#From 1000 genomes analysis folder
import os
os.mkdir("1000bootstraps")
nmales=1233
nfemales=1271
nsites= 44527
for i in range(1,1001):
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
	os.system("paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt | tail -n " +str(nsites) + " | cat header.txt - >  1000bootstraps/permutation_mafabove005_humans_"+str(i)+".txt")
	## first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq




#From the pipefish analysis folder

import os
os.mkdir("1000bootstraps")
nmales=167
nfemales=57
nsites= 44773
for i in range(1,1001):
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
	os.system("paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt | tail -n " +str(nsites) + " | cat header.txt - >  1000bootstraps/permutation_mafabove005_"+str(i)+"pipefish.txt")
## first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq




#From the flycatcher analysis folder


import os
os.mkdir("1000bootstraps")
nmales=47
nfemales=47
nsites= 95974
for i in range(1,1001):
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
	os.system("paste temp_posinfo.txt male_allele_N.txt  tempmale_p.txt  female_allele_N.txt tempfemale_p.txt | tail -n " +str(nsites) + " | cat header.txt - >  1000bootstraps/permutation_mafabove005_"+str(i)+"flycatcher.txt")
	## first line to scaf	pos	n_males_allele_covered	male_freq	n_females_allele_covered	female_freq

