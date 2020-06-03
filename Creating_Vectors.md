
Creating Vectors
================

This short script aims at creating vectors to play with. It
saves everything in one file with a README() function.

Running README() is explaining everything.

## Functions

The functions below define the way we infer and draw the theoretical
null distribution, Nei’s FST, and how we obtain quantiles.

``` r
library("ggplot2")
library("lemon")
```

   

``` r
#Nei's Fst
nei_fst_func <- function(pm,pf){((pm-pf)^2)/(4*((pm+pf)/2)*(1-((pm+pf)/2)))}

#Parameter of the chi-square distribution given nm and nf 
exp_chisq <- function(nm,nf){return(1 / ( 2 * ( 2 / ((1/nm)+(1/nf)) ) ))}

# Small function that calculate the quantile of the expected null distribution in which a SNP fall, given nm and nf and the observed allele frequencies. Usefull to compare SNP with different number of males and females
quantile_from_expdistrib<-function(pm,pf,nm,nf){
  exp_snp<-exp_chisq(nm,nf)
  quantilesnp <- pchisq(nei_fst_func(pm,pf)/exp_snp,df=1)
  return(quantilesnp)
}

#Small function that compute the percentage of data above a given quantile (i.e)
howmuchabove<-function(pm,pf,nm,nf,percent,output="empty",verbose=T){
  exp_param<-exp_chisq(nm,nf)
  n_above = length(which(nei_fst_func(pm,pf)> qchisq(p=percent,df=1)*exp_chisq(nm,nf)))
  perc_above = (n_above/length(pm))*100
  if (verbose==T){
    print(paste(n_above,"SNPs(", round(perc_above,4),"%) above percentile",percent*100,"% in the observed distribution of",length(pm)," SNPs",sep=" "))
  }
  ### Output ( can be  "empty", "number","perc")
  if (output == "empty"){}
  else if (output == "number"){
    return(n_above)}
  else if (output == "perc"){
     return(perc_above)}
  
  else{
    print(c("output should be 'number' or 'perc', not returning anything"))
  }
  }

#Small function that draws a theoretical null by taking one point from the chi-square null of every SNP
draw_one<-function(x){return(rchisq(n = 1, df =1 )*x)}

theoretical_null<-function(nm,nf){
  exp_snp<-exp_chisq(nm,nf)
  return(unlist(lapply(exp_snp,draw_one)))
}
```

The function below is the analysis pipeline applied to any dataset.
Refer to the README for information on how the allele frequencies files
where generated.

``` r
output_vectors<-function(filename){
  ###Read data
  data_observed<-read.table(filename,h=T)
  
  
  
  ####################REAL data######################################
  
  #calculate FST from allele frequencies
  fst.obs<-nei_fst_func(data_observed$male_freq,data_observed$female_freq)
  
  #calculate in which quantile each snp falls respective to its own theoretical distribution
  quantile_observed<-quantile_from_expdistrib(data_observed$male_freq,data_observed$female_freq,data_observed$n_males_allele_covered,data_observed$n_females_allele_covered)
  
  #draw one value of the theoretical distribution of each SNP as a fair null (only matters when number of sampled alleles varies SNP to SNP)
  theoretical_null_draw<-theoretical_null(data_observed$n_males_allele_covered,data_observed$n_females_allele_covered)
  
  #fill a vector with the number of obseerved SNPs that fall in each quantile
  quantiles.perc<-c()
  for (i in 1:100){
    quantiles.perc[i] <- length(which(quantile_observed<i/100. & quantile_observed>=(i-1)/100.))/length(quantile_observed)
  } 
  return(list(fst.obs,quantile_observed,quantiles.perc))
}
```

### Create Vectors

``` r
##Flycatcher
file_observed_data<-"freq_files/clean_frequencies_flycatcherMAFabove005.txt"
flydata <-output_vectors(filename = file_observed_data)
fly_fst_observed<-as.numeric(flydata[[1]])
fly_quantile_observed<-as.numeric(flydata[[2]])
fly_percquantile_observed<-as.numeric(flydata[[3]])

file_permuted_data<-"freq_files/permutation_mafabove005_1flycatcher.txt"
flydatapermuted <-output_vectors(filename = file_permuted_data)
fly_fst_permuted<-as.numeric(flydatapermuted[[1]])
fly_quantile_permuted<-as.numeric(flydatapermuted[[2]])
fly_percquantile_permuted<-as.numeric(flydatapermuted[[3]])


##pipefish
file_observed_data<-"freq_files/clean_frequencies_pipefishMAFabove005.txt"
fishdata <-output_vectors(filename = file_observed_data)
fish_fst_observed<-as.numeric(fishdata[[1]])
fish_quantile_observed<-as.numeric(fishdata[[2]])
fish_percquantile_observed<-as.numeric(fishdata[[3]])

#theoretical null of pipefish with varying 
data<-read.table(file_observed_data,h=T)
fish_theoretical_null<-theoretical_null(data$n_males_allele_covered,data$n_females_allele_covered)


file_permuted_data<-"freq_files/permutation_mafabove005_1pipefish.txt"
fishdatapermuted <-output_vectors(filename = file_permuted_data)
fish_fst_permuted<-as.numeric(fishdatapermuted[[1]])
fish_quantile_permuted<-as.numeric(fishdatapermuted[[2]])
fish_percquantile_permuted<-as.numeric(fishdatapermuted[[3]])

##human
file_observed_data<-"freq_files/clean_frequencies_humansCDSMAFabove005.txt"
humandata <-output_vectors(filename = file_observed_data)
human_fst_observed<-as.numeric(humandata[[1]])
human_quantile_observed<-as.numeric(humandata[[2]])
human_percquantile_observed<-as.numeric(humandata[[3]])

file_permuted_data<-"freq_files/permutation_mafabove005_1humans.txt"
humandatapermuted <-output_vectors(filename = file_permuted_data)
human_fst_permuted<-as.vector(humandatapermuted[[1]])
human_quantile_permuted<-as.vector(humandatapermuted[[2]])
human_percquantile_permuted<-as.vector(humandatapermuted[[3]])


## Human without maf filtering

file_observed_data<-"freq_files/clean_frequencies_humansCDS.txt"
humandata_nomafFilter <-output_vectors(filename = file_observed_data)
human_fst_observed_nomafFilter<-as.numeric(humandata_nomafFilter[[1]])
human_quantile_observed_nomafFilter<-as.numeric(humandata_nomafFilter[[2]])
human_percquantile_observed_nomafFilter<-as.numeric(humandata_nomafFilter[[3]])

file_permuted_data<-"freq_files/boot_1human.txt"
humandatapermuted_nomafFilter <-output_vectors(filename = file_permuted_data)
human_fst_permuted_nomafFilter<-as.vector(humandatapermuted_nomafFilter[[1]])
human_quantile_permuted_nomafFilter<-as.vector(humandatapermuted_nomafFilter[[2]])
human_perc_quantile_permuted_nomafFilter<-as.vector(humandatapermuted_nomafFilter[[3]])



```

### Saving object

``` r
README<-function(){
  cat("

####### SUMMARY

NAME STRUCTURE: SPECIES_DATA_SOURCE

SPECIES: fly | fish | human
DATA: fst | quantile | percquantile
SOURCE: observed | permuted

example:

fly_fst_permuted contains the fst value for each SNP in the permuted flycatcher data.

####### Description

There are  19 vectors:

for each dataset i.e.(flycatcher, pipefish, human) observed or permuted there is: 

-- Fst
-- Quantiles of  data in theoretical Fst (i.e. the quantile value for each SNP in its own theoretical distribution)
-- Percentage of data falling into each quantile ( 100 data points)

For fish, the theoretical is a bit special since there is missing data, that is the 19th vector: fish_theoretical_null

######### ALL vectors

fly_fst_observed: flycatcher fst observed
fly_quantile_observed: flycatcher quantile observed
fly_perc_quantile_observed: flycatcher % observed

fly_fst_permuted: flycatcher fst permuted
fly_quantile_permuted: flycatcher quantile permuted
fly_perc_quantile_permuted: flycatcher % permuted

fish_fst_observed: pipefish fst observed
fish_quantile_observed: pipefish quantile observed
fish_perc_quantile_observed: pipefish % observed

fish_fst_permuted: pipefish fst permuted
fish_quantile_permuted: pipefish quantile permuted
fish_perc_quantile_permuted: pipefish % permuted
fish_theoretical_null : Fish theoretical distribution accounting for missing data varying at each position

human_fst_observed: human fst  observed
human_quantile_observed: human quantile observed 
human_perc_quantile_observed: human % observed 

human_fst_permuted_nomafFilter: human fst permuted no maf filter for suppmat
human_quantile_permuted_nomafFilter: human quantile permuted no maf filter for suppmat
human_perc_quantile_permuted_nomafFilter: human % permuted no maf filter for suppmat
")
}

save(README,fly_fst_observed,fly_quantile_observed,fly_percquantile_observed,fly_fst_permuted,fly_quantile_permuted,fly_percquantile_permuted, fish_fst_observed, fish_quantile_observed, fish_percquantile_observed, fish_fst_permuted, fish_quantile_permuted, fish_percquantile_permuted, human_fst_observed, human_quantile_observed, human_percquantile_observed, human_fst_permuted, human_quantile_permuted, human_percquantile_permuted,fish_theoretical_null,human_fst_permuted_nomafFilter,human_quantile_permuted_nomafFilter,human_perc_quantile_permuted_nomafFilter , file="vectors_reanalysis.RData")
```
