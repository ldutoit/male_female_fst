Chi-square test
================
Ludovic Dutoit
3/16/2020

This document does the chi-square tests to look for enrichment in the tails of the observed distribution compared to expected null.

Functions
---------

The functions below define the way we infer and draw the theoretical null distribution, Nei's FST, and how we obtain quantiles.

``` r
library("ggplot2")
```

    ## Warning: package 'ggplot2' was built under R version 3.5.2

``` r
library("lemon")
```

    ## Warning: package 'lemon' was built under R version 3.5.2

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

The function below is the analysis pipeline applied to any dataset. Refer to the README for information on how the allele frequencies files where generated.

``` r
reanalysis<-function(file_observed_data,file_permutated_data,prefix="prefix"){
  ###Read data
  data_permuted<-read.table(file_permutated_data,h=T)
  data_observed<-read.table(file_observed_data,h=T)
  
  
  
  ####################REAL data######################################
  
  #calculate FST from allele frequencies
  fst.obsreal<-nei_fst_func(data_observed$male_freq,data_observed$female_freq)
  
  
  
  ########################################## Visualisations###################

 
 #####chisquare
  theory_quantile99 <- 1.659*((1/data_observed$n_males_allele_covered)+(1/data_observed$n_females_allele_covered))
  top_1_obs <- length(which(nei_fst_func(data_observed$male_freq,data_observed$female_freq)>=theory_quantile99))
  top_1_exp <- length(nei_fst_func(data_observed$male_freq,data_observed$female_freq))/100.
  
  #Observed and expected for bottom 99% quantiles
  bottom_99_obs<-length(which(nei_fst_func(data_observed$male_freq,data_observed$female_freq)<=theory_quantile99))
  bottom_99_exp<-  length(nei_fst_func(data_observed$male_freq,data_observed$female_freq))*0.99
  matrix1<-cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp))
  rownames(matrix1)<-c("top_1_obs","top_1_exp")
  colnames(matrix1)<-c("bottom_99_obs","bottom_99_exp")
  print("matrix 1%")
  print (matrix1)
  print(chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp))))
  
  #####chisquare
  theory_quantile95 <- qchisq(p=0.95,df=1)*exp_chisq(data_observed$n_males_allele_covered,data_observed$n_females_allele_covered)
  top_5_obs <- length(which(nei_fst_func(data_observed$male_freq,data_observed$female_freq)>=theory_quantile95))
  top_5_exp <- length(nei_fst_func(data_observed$male_freq,data_observed$female_freq))*0.05
  
  #Observed and expected for bottom 99% quantiles
  bottom_95_obs<-length(which(nei_fst_func(data_observed$male_freq,data_observed$female_freq)<=theory_quantile95))
  bottom_95_exp<-  length(nei_fst_func(data_observed$male_freq,data_observed$female_freq))*0.95
  matrix5<-cbind(c(top_5_obs,top_5_exp),c(bottom_95_obs,bottom_95_exp))
  rownames(matrix5)<-c("top_5_obs","top_5_exp")
  colnames(matrix5)<-c("bottom_95_obs","bottom_95_exp")
  print("matrix 5%")
  print (matrix5)
  print(chisq.test(cbind(c(top_5_obs,top_5_exp),c(bottom_95_obs,bottom_95_exp))))
  }
```

### Humans (&gt;maf.0.05)

1233 males, 1271 females, 44,527 sites

``` r
file_permutated_data<-"freq_files/permutation_mafabove005_1humans.txt"
file_observed_data<-"freq_files/clean_frequencies_humansCDSMAFabove005.txt"
reanalysis(file_observed_data,file_permutated_data,prefix="Humans")
```

    ## [1] "matrix 1%"
    ##           bottom_99_obs bottom_99_exp
    ## top_1_obs        405.00      44122.00
    ## top_1_exp        445.27      44081.73
    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  cbind(c(top_1_obs, top_1_exp), c(bottom_99_obs, bottom_99_exp))
    ## X-squared = 1.8312, df = 1, p-value = 0.176
    ## 
    ## [1] "matrix 5%"
    ##           bottom_95_obs bottom_95_exp
    ## top_5_obs       2162.00      42365.00
    ## top_5_exp       2226.35      42300.65
    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  cbind(c(top_5_obs, top_5_exp), c(bottom_95_obs, bottom_95_exp))
    ## X-squared = 0.96192, df = 1, p-value = 0.3267

### Flycatcher (maf&gt;0.05)

47 males, 47 females, 95974 sites

``` r
file_permutated_data<-"freq_files/permutation_mafabove005_1flycatcher.txt"
file_observed_data<-"freq_files/clean_frequencies_flycatcherMAFabove005.txt"
reanalysis(file_observed_data,file_permutated_data,prefix="Flycatcher")
```

    ## [1] "matrix 1%"
    ##           bottom_99_obs bottom_99_exp
    ## top_1_obs       1745.00      94229.00
    ## top_1_exp        959.74      95014.26
    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  cbind(c(top_1_obs, top_1_exp), c(bottom_99_obs, bottom_99_exp))
    ## X-squared = 230.65, df = 1, p-value < 2.2e-16
    ## 
    ## [1] "matrix 5%"
    ##           bottom_95_obs bottom_95_exp
    ## top_5_obs        5376.0       90598.0
    ## top_5_exp        4798.7       91175.3
    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  cbind(c(top_5_obs, top_5_exp), c(bottom_95_obs, bottom_95_exp))
    ## X-squared = 34.469, df = 1, p-value = 4.331e-09

### Pipefish (maf&gt;0.05)

Pipefish has significantly more missing data. We therefore kept any site with less than 50% missing data.

167 males, 57 femailes, 44'773 sites

``` r
file_permutated_data<-"freq_files/permutation_mafabove005_1pipefish.txt"
file_observed_data<-"freq_files/clean_frequencies_pipefishMAFabove005.txt"
reanalysis(file_observed_data,file_permutated_data,prefix="Pipefish")
```

    ## [1] "matrix 1%"
    ##           bottom_99_obs bottom_99_exp
    ## top_1_obs       1125.00      43648.00
    ## top_1_exp        447.73      44325.27
    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  cbind(c(top_1_obs, top_1_exp), c(bottom_99_obs, bottom_99_exp))
    ## X-squared = 295.99, df = 1, p-value < 2.2e-16
    ## 
    ## [1] "matrix 5%"
    ##           bottom_95_obs bottom_95_exp
    ## top_5_obs       3774.00      40999.00
    ## top_5_exp       2238.65      42534.35
    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  cbind(c(top_5_obs, top_5_exp), c(bottom_95_obs, bottom_95_exp))
    ## X-squared = 419.73, df = 1, p-value < 2.2e-16
