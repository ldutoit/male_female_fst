# Overall FST

First we extract  Fst from each permutations and store mean and medians, before inferring the 95%CI from quantiles.
```r
nei_fst_func <- function(pm,pf){((pm-pf)^2)/(4*((pm+pf)/2)*(1-((pm+pf)/2)))}


paths<-c("1000genomes/1000bootstraps/","pipefish/1000bootstraps/", "flycatcher/1000bootstraps/")

for (path in paths){
	n_perm = 0
	means<-c()
	medians<-c()
	for (permutation in dir(path)){
		n_perm = n_perm +1
		if (n_perm%%10==0){print(n_perm)}
		#print(c(path,permutation))
		perm<-read.table(paste(path,permutation,sep=""),h=T)
		fst<-nei_fst_func(as.numeric(as.character(perm$male_freq)),as.numeric(as.character(perm$female_freq)))

		means<-c(means, mean(fst))
		medians<-c(medians, median(fst))
	}
		print(path)
		print("means")
		print(quantile(means,c(0.025,0.975)))
		print("medians")
		print(quantile(medians,c(0.025,0.975)))
}		#

#[1] "1000genomes/1000bootstraps/"
#[1] "means"
#        2.5%        97.5%
#0.0001745306 0.0003184714
#[1] "medians"
#        2.5%        97.5%
#7.761167e-05 1.520658e-04#
#

#[1] "pipefish/1000bootstraps/"
#[1] "means"
#       2.5%       97.5%
#0.003513656 0.003907291
#[1] "medians"
#       2.5%       97.5%
#0.001507340 0.001693004#

[1] "flycatcher/1000bootstraps/"
[1] "means"
       2.5%       97.5%
0.005099666 0.005797476
[1] "medians"
       2.5%       97.5%
0.002247276 0.002566306

```
#NOW the original values on observed data

```r
nei_fst_func <- function(pm,pf){((pm-pf)^2)/(4*((pm+pf)/2)*(1-((pm+pf)/2)))}

data<-read.table("freq_files/clean_frequencies_flycatcherMAFabove005.txt",h=T)
fst<-nei_fst_func(as.numeric(data$male_freq),as.numeric(data$female_freq))
mean(fst)
#[1] 0.007550288
median(fst)
#[1] 0.002318819


data<-read.table("freq_files/clean_frequencies_humansCDSMAFabove005.txt",h=T)
fst<-nei_fst_func(as.numeric(data$male_freq),as.numeric(data$female_freq))
mean(fst)
#[1] 0.0002004794
median(fst)
#[1]  9.578959e-05



data<-read.table("freq_files/clean_frequencies_pipefishMAFabove005.txt",h=T)
fst<-nei_fst_func(as.numeric(data$male_freq),as.numeric(data$female_freq))
mean(fst)
#[1] 0.004555063
median(fst)
#[1]  0.001875247
```

so in summary:


**Flycatcher**
mean observed Fst: 0.00756 (95%CI of the permutations: 0.005099666 -  0.005797476) SIGNIFICANT  
median observed Fst: 0.00232  (95%CI of the permutations: 0.00225 0.002576306 ) NOT SIGNIFICANT

**Pipefish**
mean observed Fst: 0.0046 (95%CI of the permutations: 0.0035 - 0.0039) SIGNIFICANT
median observed Fst: 0.0019  (95%CI of the permutations: 0.0015 0.0017#) SIGNIFICANT 


**Humans**

mean observed Fst:0.00020 (95%CI of the permutations: 0.00017 - 0.00031)  NOT SIGNIFICANT 
median observed Fst:  9.579e-05  (95%CI of the permutations: 7.761-05 1.5202-04) NOT SIGNIFICANT 



