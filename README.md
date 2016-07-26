
## -----------whiteOaksABC-----------


 This is a GitHub respository containing scripts and programs used to reconstruct the demographic history of white oaks using an approximate Bayesian computation framework (submitted paper)
 
 For more details or any help, please contact me: thibault.leroy_at_pierroton.inra.fr


### 1/ PROGRAMS

To compile mscalc, the calculator of summary statistics (Ross-Ibarra et al. 2008; 2009; Roux et al. 2011):

gcc *.c -lm -o mscalc

To compile msnsam, the coalescent simulator (Ross-Ibarra et al. 2008):

./clms

Priorgen.py generates prior for 4 different models of speciation (SI, AM, IM, SC)

priorgen.py -h




### 2/ DATASETS

 All datasets used to perform the ABC are available. For each pair of species, a bpfile, a spinput.txt and a target file are avaible. bpfile & spinput need to be in the directory where you use the coalescent simulator. 
 
 The target file contains the 19 summary statistics calculated on each real data set.




### 3/ EXAMPLE: Multilocus coalescent simulations:

- general information concerning the script

25000 multilocus simulations assuming an AM scenario between Q. robur & Q. petraea [i.e. number of SNPs (=3304) x number of simulations (=25000) = 82600000]

Number of SNPs: 2nd line of the spinput.txt file

Number of simulations: penultimate line of the spinput.txt file



- bash script (note that programs are assumed to be in your bin directory):

mknod myfifo p

mscalc < myfifo &

priorgen.py bpfile=bpfile n1=0 n1=100 n2=0 n2=100 nA=0 nA=100 tau=0 tau=100 M1=0 M1=100 M2=0 M2=100 shape1=0 shape1=100 shape2=0 shape2=500 model=AM nreps=25000 Nvariation=homo Mvariation=hetero symMig=asym parameters=priorfile | msnsam tbs 82600000 -s 1 -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs >myfifo




### 4/ EXAMPLE: ABC Model Choice (Rscript):

here for each model, we assume : 400 directories containing an ABCstat.txt file in which we have 25,000 lines of summary statistics =10,000,000 simulations/model

library(nnet)
- import source cv4ABC (Csillery et al. 2012)

source("cv4abc.R")

- import observed summary statistics

setwd("./robur-petraea/")

target=read.table("target_rob-pet.txt",skip=2,h=F)

ss=c(2:20)

- import each line of summary statistics corresponding to your simulations

M1heteroNhomoM=M1heteroNheteroM=M2heteroNhomoM=M2heteroNheteroM=M3heteroNhomoM=M3heteroNheteroM=M4heteroNhomoM=M4heteroNheteroM=NULL

for(i in 1:400){

  M1heteroNhomoM=rbind(M1heteroNhomoM, matrix(scan(paste("rob-pet_strictisolation_Nhetero_Mhomo_v1_",i,"/ABCstat.txt", sep=""), skip=2), byrow=T, ncol=20)[,ss])

  M2homoNheteroM=rbind(M2homoNheteroM, matrix(scan(paste("rob-pet_island_Nhomo_Mhetero_v1_",i,"/ABCstat.txt", sep=""), skip=2), byrow=T, ncol=20)[,ss])

  M3homoNheteroM=rbind(M3homoNheteroM, matrix(scan(paste("rob-pet_ancientmig_Nhomo_Mhetero_v1_",i,"/ABCstat.txt", sep=""), skip=2), byrow=T, ncol=20)[,ss])

  M4homoNheteroM=rbind(M4homoNheteroM, matrix(scan(paste("rob-pet_2ndarycont_Nhomo_Mhetero_v1_",i,"/ABCstat.txt", sep=""), skip=2), byrow=T, ncol=20)[,ss])

}

- replace all "NaN" [!!! The number of lines of summary statistics need to be the same for the 4 models !!!]

for(i in 1:ncol(M1heteroNhomoM)){

  M1heteroNhomoM[which(M1heteroNhomoM[,i]=="NaN"),i]=mean(M1heteroNhomoM[,i], na.rm=T)

  M2homoNheteroM[which(M2homoNheteroM[,i]=="NaN"),i]=mean(M2homoNheteroM[,i], na.rm=T)

  M3homoNheteroM[which(M3homoNheteroM[,i]=="NaN"),i]=mean(M3homoNheteroM[,i], na.rm=T)

  M4homoNheteroM[which(M4homoNheteroM[,i]=="NaN"),i]=mean(M4homoNheteroM[,i], na.rm=T)

}


- generate a long vector of numbers from the simulations (1=model1, 2=model2..., 4=model4), used as the dependent variable for the regression

x=c(rep(1:4, each=nrow(M1heteroNhomoM)))


- to perform several ABC analyses (here 100), duplicate your real dataset

obs=matrix(rep(target[ss],100), byrow=T, nrow=100)

- then perform your ABC analysis [note that tol is the most important paramters = required proportion of points nearest the target values (here 10,000/40,000,000 = 0.00025 best simulations)]

res=model_selection_abc_nnet(target=obs, x=x, sumstat=rbind(M1heteroNhomoM,M2homoNheteroM,M3homoNheteroM,M4homoNheteroM), tol=10000/(4*nrow(M1heteroNhomoM)), noweight=F, rejmethod=F, nb.nnet=20, size.nnet=8, output="OBS_rob-pet_SI_withHeteroNe_vs_IM_AM_SC_withHomoNeHeteroM_100ABC_tol10000_040116")




### 5/ EXAMPLE: generate posterior distribution under the best model (Rscript):

library(nnet)

- import source cv4estimations.R ((Csillery et al. 2012))

source("cv4estimations.R")


- import summary statistics (real data set)

setwd("./robur-petraea2/")

ss=c(2:20)

ncv=10 #number of validations to perform

tmp=as.numeric(read.table("target_rob-pet.txt",skip=2)[ss])

target=NULL

for(i in 1:ncv){target=rbind(target, tmp)}

- import all summary statistics (simulated data sets) [edit the number of directory (for i in 1:x..) to compile all your priorfile & ABCstat.txt files], please check the number of columns in your files!

prior=NULL

tmp.stat=NULL

tmp.parameters=NULL

for(i in 1:1){ 

  tmp.stat=matrix(scan("sumup_ABCstat_SC_homo_hetero_2.txt"), byrow=T, ncol=20)[,ss] 

  tmp.parameters=matrix(scan("sumup_priorfile_SC_homo_hetero.txt"), byrow=T, ncol=13) # ncol à modifier

  if(nrow(tmp.stat)==nrow(tmp.parameters)){

    prior.tmp=cbind(tmp.stat, tmp.parameters)

    prior=rbind(prior, prior.tmp)
  }

}

- then generate posteriors 

prior=na.omit(prior)

abc_nnet_multivar(target=target, x=prior[, -(1:length(ss))], sumstat=prior[, 1:length(ss)], tol=10000/(nrow(prior)), 
rejmethod=F, noweight=F, transf=rep("logit", ncol(prior)-length(ss)), bb=rbind(apply(prior[,-(1:length(ss))], MARGIN=2, FUN="min"), apply(prior[,-(1:length(ss))], MARGIN=2, FUN="max")), nb.nnet=25, size.nnet=10, trace=T, output="Posterior_10000_SC_MHETERO_v2_")
