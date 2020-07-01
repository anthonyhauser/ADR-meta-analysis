#Model_estimation
#selection of variables
#imputation of time on ART missing values, stepwise selection

#########################################################################################################################
library(rstan) #https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
library(gtools)
library(gdata)
library(bayesplot)
library("R.utils")
source('C:/Users/ahauser/Desktop/systematic_review_AH2/R/stan/stan_utility.R')
#Stan configuration
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
###############################################################################################################################################################################
###############################################################################################################################################################################

#Data management
load(file="C:/Users/ahauser/Desktop/systematic_review_AH2/R/meta_analysis/sr.RData")
#-------------------------------------------------------------------------------------------------------------
#Aggregate results from same study
sr=sr[,-c(which(colnames(sr) %in% c("inclusion","PID","RT_seq","PR_seq","Seq_Tech")),
          which(colnames(sr)=="nb_on_treatment"):which(colnames(sr)== "subtyp_c"))]
#first
f<-1:which(colnames(sr)=="repeated")
#weighting average
wa<-c(which(colnames(sr) %in% c("female_perc","Time_to_fail")),
      which(colnames(sr)=="ABC"):which(colnames(sr)=="PIdrug"))
#sum
#s<-which(colnames(sr)=="Genotyped"):which(colnames(sr)=="G190_AERS")
s<-which(colnames(sr)=="Genotyped"):which(colnames(sr)=="TAM")
double_r=lapply(as.list(sort(unique(sr$record_id))),function(x) which(sr$record_id==x))
double_r_id=which(sapply(double_r,function(x) length(x)>1))
double_r_loc=double_r[double_r_id]
nb_row_removed=0
for(i in 1:length(double_r_id)){
  mat=sr[(double_r_loc[[i]]-nb_row_removed),]
  v=mat[1,]
  v[wa]=apply(mat[,wa],2,function(x) sum(mat[,"Genotyped"]*x)/sum(mat[,"Genotyped"]))
  v[s]=apply(mat[,s],2,function(x) sum(x))
  v["Publication"]=substr(as.character(v["Publication"]),1,nchar(as.character(v["Publication"]))-1)
  sr=sr[-((double_r_loc[[i]])[-1]-nb_row_removed),]
  sr[(double_r_loc[[i]])[1]-nb_row_removed,]=v
  nb_row_removed=nb_row_removed+length(double_r_loc[[i]])-1
}

#statistics
#genotyped
which_na_nnrti<-apply(sr[,colnames(sr) %in% c("EFV","NVP")],1,function(x) sum(!is.na(x))==0)
which_na_nrti<-apply(sr[,colnames(sr) %in% c("FTC.3TC","ddI","d4T","TDF","ZDV")],1,function(x) sum(is.na(x))>0)
sum(!which_na_nrti)
sum(!which_na_nnrti)
sum(!which_na_nrti | !which_na_nnrti)

sum(sr[!which_na_nrti | !which_na_nnrti,]$Genotyped)
median(sr[!which_na_nrti | !which_na_nnrti,]$Genotyped)
min(sr[!which_na_nrti | !which_na_nnrti,]$Genotyped)
max(sr[!which_na_nrti | !which_na_nnrti,]$Genotyped)
#drugs
apply(sr[,colnames(sr) %in% c("FTC.3TC","ddI","d4T","TDF","ZDV","EFV","NVP")],2,function(x){
  sum(x[!is.na(x)]*sr[!which_na & !is.na(x),]$Genotyped)/sum(sr[!which_na & !is.na(x),]$Genotyped)
})
apply(sr[,colnames(sr) %in% c("FTC.3TC","d4T","TDF","ZDV")],2,function(x){
  sum(x[!is.na(x)]*sr[!which_na & !is.na(x),"FTC.3TC"]/10000*sr[!which_na & !is.na(x),]$Genotyped)/sum(sr[!which_na & !is.na(x),]$Genotyped)
})
apply(sr[!which_na,colnames(sr) %in% c("FTC.3TC","ddI","d4T","TDF","ZDV","EFV","NVP")],2,function(x) sum(x*sr[!which_na,]$Genotyped)/sum(sr[!which_na,]$Genotyped))

#-------------------------------------------------------------------------------------------------------------
#NRTI mutations data
#retrieve nrti drugs
#nrti_pos=which(colnames(sr)=="M41_L"):which(colnames(sr)=="K219_DENQR")
nrti_pos=c(which(colnames(sr)=="M41_L"):which(colnames(sr)=="K219_DENQR"),which(colnames(sr)=="TAM"))
nrti_names<-colnames(sr)[nrti_pos]
x_nrti=as.matrix(cbind(sr$Time_to_fail,sr$'FTC.3TC',sr$ddI,sr$d4T,sr$TDF,sr$ZDV))
colnames(x_nrti)<-c("Time_to_fail","FTC.3TC","ddI","d4T","TDF","ZDV")

#remove studies with NA drugs
na_drug_x=apply(x_nrti[,-1],1,function(x) sum(is.na(x))==0)
x_nrti=x_nrti[na_drug_x,]
sr_nrti=sr[na_drug_x,]

#save backtransformed data (time in months, drug in proportion)
x_nrti[,-grep("Time_to_fail",colnames(x_nrti))]=x_nrti[,-grep("Time_to_fail",colnames(x_nrti))]/100
x_nrti_bt=x_nrti

#k, n
n_nrti=sr_nrti$Genotyped
k_nrti=sr_nrti[,nrti_pos]
k_nrti=as.data.frame(k_nrti)
k_nrti[is.na(k_nrti)]=0 #NA is considered as 0

#Dimensions
N_nrti=length(n_nrti)
M_nrti=dim(k_nrti)[2]

#-------------------------------------------------------------------------------------------------------------
#Select data and variables
#Data
k=k_nrti[,-which(colnames(k_nrti)=="TAM")] #remove TAMs for now
nrti_names=colnames(k_nrti)
N=N_nrti #number of studies
M=M_nrti-1 #number of nrti mutations, -1 as we removed TAMs
n=n_nrti #number of genotyped by study

#Time
time_bt = x_nrti[,"Time_to_fail"] #time backtransformed (real time)
time_bt[(which(is.na(time_bt)))[2]] = 4.1 * 12 #ART duration is 4.1 years for people with VL>80, assumption duration same for people with VL>1000, see Labhart
time_nna = time_bt[!is.na(time_bt)] #times with no NA
mean_t = mean(time_nna) #mean time
var_t = var(time_nna) #var time
na_time_pos = which(is.na(time_bt)) #position of NA values
nna_time_pos = which(!is.na(time_bt)) #position of non NA values
time_tr=(time_bt-mean_t)/sqrt(var_t) #transformed times
time_tr_nna=time_tr[nna_time_pos] #transformed times with no NA
#y: other variables: NRTI drugs
y_nrti_bt = x_nrti_bt[,-1] #removing time
y_mean = apply(y_nrti_bt,2,mean) #mean by NRTI drugs
y_sd = apply(y_nrti_bt,2,sd) #sd by NRTI drugs
y_nrti = (y_nrti_bt - structure(rep(y_mean,each=dim(y_nrti_bt)[1]),dim=dim(y_nrti_bt)))/ #transformed NRTI drugs
  structure(rep(y_sd,each=dim(y_nrti_bt)[1]),dim=dim(y_nrti_bt))
#time_y: time binded with y
time_y_bt=cbind(time_bt,y_nrti_bt) #backtransformed time_y
time_y_bt[na_time_pos,1]=0 #Replace NA by 0 as Stan does not accept NA values
time_y=cbind(time_tr,y_nrti) #transformed time_y
time_y[na_time_pos,1]=0 #Replace NA by 0 as Stan does not accept NA values
x_mean=c(mean_t,y_mean) #Mean by variables
x_sd=c(sqrt(var_t),y_sd) #sd by variables



###############################################################################################################################################################################
###############################################################################################################################################################################
#Model the prevalence of NRTI mutations
mod_NRTI<-stan_model(file="C:/Users/ahauser/Desktop/systematic_review_AH2/R/stan/TAM/NRTI_stan3.stan") #does not calculate p_x_study
mod_NRTI<-stan_model(file="C:/Users/ahauser/Desktop/systematic_review_AH2/R/stan/TAM/NRTI_stan2.stan") #calculate p_x_study

#Starting setting
B_tot=dim(time_y_bt)[2]
B_list_assumed=rep(list(1:B_tot),M)
#B_list_assumed=list(c(1,2,6),c(1,3,4,5),c(1,3,4),c(1,3,5),c(1,2,4),c(1,2,6),c(1,3,4,6),c(1))
n_run=0
n_change_tot=1
last_method=FALSE
t_lim=800 #time limit before considering as a failure

#Stepwise method
while(n_run<15 & (n_change_tot>0 | !last_method)){
  y_list<-lapply(B_list_assumed,function(x){
    return(structure(time_y[,x],dim=c(dim(time_y)[1],length(x))))
  })
  y_bt_list<-lapply(B_list_assumed,function(x){
    return(structure(time_y_bt[,x],dim=c(dim(time_y_bt)[1],length(x))))
  })
  
  #Dimensions
  B=sapply(B_list_assumed,function(x) length(x))
  B_tot=dim(time_y_bt)[2]
  
  for(i in 1:M){
    assign(paste("x_m",i,sep=""),y_list[[i]])
    assign(paste("x_study",i,sep=""),y_bt_list[[i]])
    if(B[i]>1) prov=c(36,as.numeric((B_list_assumed[[i]])[-1] %in% c(2,5))) else prov=36
    assign(paste("x_FTC_TDF",i,sep=""),prov)
    if(B[i]>1) prov=c(36,as.numeric((B_list_assumed[[i]])[-1] %in% c(2,6))) else prov=36
    assign(paste("x_FTC_ZDV",i,sep=""),prov)
  }
  
  #Model estimation
  stan_NRTI=list(N=N,
                 M=M,
                 B=B,
                 B_tot=B_tot,
                 k=matrix(round(as.matrix(k)),nrow=N),
                 n=n,
                 tam_pos = setdiff(1:9,which(nrti_names %in% c("K65_R","M184_IV","TAM"))),
                 time_tr_nna=time_tr_nna, #ART time in non-missing position
                 mean_t=mean_t, #mean of ART time
                 var_t =var_t, #variance of ART time
                 na_time_pos=na_time_pos, #position of missing ART times
                 nna_time_pos=nna_time_pos, #position of non-missing ART times
                 x_m1=x_m1,
                 x_m2=x_m2,
                 x_m3=x_m3,
                 x_m4=x_m4,
                 x_m5=x_m5,
                 x_m6=x_m6,
                 x_m7=x_m7,
                 x_m8=x_m8,
                 x1_mean=structure(x_mean[B_list_assumed[[1]]],dim=B[1]),
                 x2_mean=structure(x_mean[B_list_assumed[[2]]],dim=B[2]),
                 x3_mean=structure(x_mean[B_list_assumed[[3]]],dim=B[3]),
                 x4_mean=structure(x_mean[B_list_assumed[[4]]],dim=B[4]),
                 x5_mean=structure(x_mean[B_list_assumed[[5]]],dim=B[5]),
                 x6_mean=structure(x_mean[B_list_assumed[[6]]],dim=B[6]),
                 x7_mean=structure(x_mean[B_list_assumed[[7]]],dim=B[7]),
                 x8_mean=structure(x_mean[B_list_assumed[[8]]],dim=B[8]),
                 x1_sd=structure(x_sd[B_list_assumed[[1]]],dim=B[1]),
                 x2_sd=structure(x_sd[B_list_assumed[[2]]],dim=B[2]),
                 x3_sd=structure(x_sd[B_list_assumed[[3]]],dim=B[3]),
                 x4_sd=structure(x_sd[B_list_assumed[[4]]],dim=B[4]),
                 x5_sd=structure(x_sd[B_list_assumed[[5]]],dim=B[5]),
                 x6_sd=structure(x_sd[B_list_assumed[[6]]],dim=B[6]),
                 x7_sd=structure(x_sd[B_list_assumed[[7]]],dim=B[7]),
                 x8_sd=structure(x_sd[B_list_assumed[[8]]],dim=B[8]),
                 x_FTC_TDF1=structure(x_FTC_TDF1,dim=B[1]),
                 x_FTC_TDF2=structure(x_FTC_TDF2,dim=B[2]),
                 x_FTC_TDF3=structure(x_FTC_TDF3,dim=B[3]),
                 x_FTC_TDF4=structure(x_FTC_TDF4,dim=B[4]),
                 x_FTC_TDF5=structure(x_FTC_TDF5,dim=B[5]),
                 x_FTC_TDF6=structure(x_FTC_TDF6,dim=B[6]),
                 x_FTC_TDF7=structure(x_FTC_TDF7,dim=B[7]),
                 x_FTC_TDF8=structure(x_FTC_TDF8,dim=B[8]),
                 x_FTC_ZDV1=structure(x_FTC_ZDV1,dim=B[1]),
                 x_FTC_ZDV2=structure(x_FTC_ZDV2,dim=B[2]),
                 x_FTC_ZDV3=structure(x_FTC_ZDV3,dim=B[3]),
                 x_FTC_ZDV4=structure(x_FTC_ZDV4,dim=B[4]),
                 x_FTC_ZDV5=structure(x_FTC_ZDV5,dim=B[5]),
                 x_FTC_ZDV6=structure(x_FTC_ZDV6,dim=B[6]),
                 x_FTC_ZDV7=structure(x_FTC_ZDV7,dim=B[7]),
                 x_FTC_ZDV8=structure(x_FTC_ZDV8,dim=B[8]),
                 init_res_1=rep(-4,M),
                 init_res_2=rep(4,M),
                 #init_add=structure(init_add,dim=M),
                 sigma_1=rep(1,M),
                 tau_1=1,
                 beta_m1_1=structure(rep(0,B[1]),dim=B[1]),
                 beta_m2_1=structure(rep(0,B[2]),dim=B[2]),
                 beta_m3_1=structure(rep(0,B[3]),dim=B[3]),
                 beta_m4_1=structure(rep(0,B[4]),dim=B[4]),
                 beta_m5_1=structure(rep(0,B[5]),dim=B[5]),
                 beta_m6_1=structure(rep(0,B[6]),dim=B[6]),
                 beta_m7_1=structure(rep(0,B[7]),dim=B[7]),
                 beta_m8_1=structure(rep(0,B[8]),dim=B[8]),
                 beta_m1_2=structure(rep(4,B[1]),dim=B[1]),
                 beta_m2_2=structure(rep(4,B[2]),dim=B[2]),
                 beta_m3_2=structure(rep(4,B[3]),dim=B[3]),
                 beta_m4_2=structure(rep(4,B[4]),dim=B[4]),
                 beta_m5_2=structure(rep(4,B[5]),dim=B[5]),
                 beta_m6_2=structure(rep(4,B[6]),dim=B[6]),
                 beta_m7_2=structure(rep(4,B[7]),dim=B[7]),
                 beta_m8_2=structure(rep(4,B[8]),dim=B[8]),
                 inference=1)
  
  fit_NRTI=NULL
  tryCatch({
    withTimeout({
      fit_NRTI<- sampling(object=mod_NRTI,
                          data = stan_NRTI,
                          warmup = 1000,
                          iter = 2500,
                          chains = 4,
                          cores = 4,
                          thin = 1,
                          control=list(adapt_delta=0.99,
                                       max_treedepth=15))
    }, timeout = t_lim)
  },error=function(e) e)
  
  if(is.null(fit_NRTI)){break}
  
  
  # mod1<-stan_model(file ="C:/Users/ahauser/Desktop/systematic_review_AH2/R/stan/stan_1.stan")
  # fit_NRTI2<- sampling(object=mod1,
  #                     data = stan_NRTI,
  #                     warmup = 1000,
  #                     iter = 2500,
  #                     chains = 4,
  #                     cores = 4,
  #                     thin = 1,
  #                     control=list(adapt_delta=0.99,
  #                                  max_treedepth=10))
  

  
  #Variable selection
  d_all<-summary(fit_NRTI)$summary
  d_all<-d_all[grep("beta_bt_m",rownames(d_all)),]
  d_matrix<-as.matrix(fit_NRTI)
  d_matrix<-d_matrix[,grep("beta_bt_m",colnames(d_matrix))]
  
  n_change=0
  print("Variables kept")
  if(n_run>=2 & n_change_tot==0){last_method=TRUE}else{last_method=FALSE}
  for(i in 1:M){
    print(i)
    d<-d_all[grep(paste("beta_bt_m",i,sep=""),rownames(d_all)),]
    if(sum(c(2,3) %in% B_list_assumed[[i]])>0) var1<-1:sum(c(2,3) %in% B_list_assumed[[i]])+1 else var1<-NULL
    if(sum(c(4:6) %in% B_list_assumed[[i]])>0) var2<-1:sum(c(4:6) %in% B_list_assumed[[i]])+1+length(var1) else var2<-NULL
    var_test1=NULL
    var_test2=NULL
    if(length(var1)>0){
      #var_test1<-d[var1,"25%"]+(d[var1,"50%"]-min(d[B_list_assumed[[i]],"50%"]))>0
      if(length(var1)>1){
        var_test1=rep(TRUE,length(var1))
        var_test1[which(d[var1,"50%"]==min(d[var1,"50%"]))]=FALSE
      }else{
        var_test1<-d[var1,"25%"]+(d[var1,"50%"]-min(d[2:length(B_list_assumed[[i]]),"50%"]))>0
      }
    }
    if(length(var2)>0){
      #var_test2<-d[var2,"25%"]+(d[var2,"50%"]-min(d[B_list_assumed[[i]],"50%"]))>0
      if(length(var2)>2){
        var_test2=rep(TRUE,length(var2))
        var_test2[which(d[var2,"50%"]==min(d[var2,"50%"]))]=FALSE
      }else{
        var_test2<-d[var2,"25%"]+(d[var2,"50%"]-min(d[2:length(B_list_assumed[[i]]),"50%"]))>0
      }
    }
    #if(n_run>=2 & n_change_tot==0){last_method=TRUE}else{last_method=FALSE}
    if(length(c(var1,var2))>0){
      print(c(TRUE,var_test1,var_test2))
      B_list_assumed[[i]]=(B_list_assumed[[i]])[c(TRUE,var_test1,var_test2)]
      n_change=n_change+sum(!c(var_test1,var_test2))
    }else{
      print("-------------")
    }
  }
  n_change_tot=n_change
  print(paste("number of changes in run", n_run))
  print(n_change_tot)
  print("Assumed variables")
  print(B_list_assumed)
  print("****************************************************************************")
  n_run<-n_run+1
}
B_list_assumed_NRTI <- B_list_assumed
d<-summary(fit_NRTI)$summary
save(B_list_assumed_NRTI,fit_NRTI,file="C:/Users/ahauser/Desktop/systematic_review_AH2/R/stan/TAM/fit_NRTI.RData")

###############################################################################################################################################################################
load("C:/Users/ahauser/Desktop/systematic_review_AH2/R/stan/TAM/fit_NRTI.RData")
d<-summary(fit_NRTI)$summary
d1<-d[setdiff(grep("p_init",rownames(d)),grep("p_init_pred",rownames(d))),]
beta<-d[grep("beta_bt",rownames(d)),]
score_comp = data.frame(beta_est=numeric(),
                        prev_change=numeric(),
                        beta_over0=logical(),
                        beta_score=numeric())

for(i in 1:8){
  beta_m=rep(0,5)
  beta_over0=rep(F,5)
  if(length(B_list_assumed_NRTI[[i]])>1){
    beta_m[(B_list_assumed_NRTI[[i]])[-1]-1] = (beta[grep(paste("beta_bt_m",i,sep=""),rownames(beta)),"50%"])[-1]
    beta_over0[(B_list_assumed_NRTI[[i]])[-1]-1] = ((beta[grep(paste("beta_bt_m",i,sep=""),rownames(beta)),"2.5%"])[-1])>0
    print((B_list_assumed_NRTI[[i]])[-1]-1)
  }
  score_comp[1:5+5*(i-1),"prev_change"]=inv.logit(logit(d1[i,"50%"])+beta_m)-d1[i,"50%"]
  score_comp[1:5+5*(i-1),"beta_est"]=beta_m
  score_comp[1:5+5*(i-1),"beta_over0"]=beta_over0
}

score_comp[,"beta_score"]<-as.vector(matrix(c(0,15,0,10,60,0,0,0,
                                              10,30,5,10,10,10,10,5,
                                              15,30,10,15,-10,15,30,10,
                                              5,30,5,15,-10,5,10,5,
                                              15,-10,10,-5,-10,15,30,10),nrow=5,byrow=TRUE))
plot(beta_score,beta_est)
plot(beta_score,prev_change)
boxplot(beta_score ~ beta_over0,data=score_comp)
cor(beta_score,prev_change)

p <- ggplot(score_comp, aes(x=beta_over0, y=beta_score)) + 
  geom_violin() + 
  theme_bw() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8)
###############################################################################################################################################################################
###############################################################################################################################################################################
#Model the prevalence of any TAM

d<-as.data.frame(fit_NRTI)
mu=array(NA,dim=c(sum(k_nrti$TAM>0),2))
cov=array(NA,dim=c(sum(k_nrti$TAM>0),2,2))
for(i in 1:sum(k_nrti$TAM>0)){
  v1<-logit(d[,paste("p_noTAM_min[",which(k_nrti$TAM>0)[i],"]",sep="")])
  v2<-logit(d[,paste("p_noTAM_max[",which(k_nrti$TAM>0)[i],"]",sep="")])
  mu[i,]<-c(mean(v1),mean(v2))
  cov[i,,]<-cov(cbind(v1,v2))
}
v1<-logit(d[,"p_noTAM_TDF_min"])
v2<-logit(d[,"p_noTAM_TDF_max"])
mu_TDF<-c(mean(v1),mean(v2))
cov_TDF<-cov(cbind(v1,v2))

v1<-logit(d[,"p_noTAM_ZDV_min"])
v2<-logit(d[,"p_noTAM_ZDV_max"])
mu_ZDV<-c(mean(v1),mean(v2))
cov_ZDV<-cov(cbind(v1,v2))

v1<-logit(d[,"p_noTAM_init_min"])
v2<-logit(d[,"p_noTAM_init_max"])
mu_init<-c(mean(v1),mean(v2))
cov_init<-cov(cbind(v1,v2))

# a=mvrnorm(1000,mu=c(1,1), Sigma=matrix(c(2,1,1,3),nrow=2))
# var(a[,1])
# cov(a[,1],a[,2])


mod_TAM<-stan_model(file="C:/Users/ahauser/Desktop/systematic_review_AH2/R/stan/TAM/stan_TAM3.stan")
stan_TAM=list(N=sum(k_nrti$TAM>0),
              M=6,
              k=round(k_nrti[k_nrti$TAM>0,"TAM"]),
              n=round(n_nrti[k_nrti$TAM>0]),
              mu=mu,
              cov = cov,
              mu_TDF= mu_TDF,
              cov_TDF = cov_TDF,
              mu_ZDV = mu_ZDV,
              cov_ZDV = cov_ZDV,
              mu_init = mu_init,
              cov_init = cov_init,
              inference=1) #max prevalence of no TAM, min of (1 - estimated prevalence with 100% FTC and TDF, after 3 years)

fit_TAM<- sampling(object=mod_TAM,
                   data = stan_TAM,
                   warmup = 1000,
                   iter = 2500,
                   chains = 4,
                   cores = 4,
                   thin = 1,
                   control=list(adapt_delta=0.99,
                                max_treedepth=10))
d4<-summary(fit_TAM)$summary
d4[grep("alpha",rownames(d4)),]
save(fit_TAM,file="C:/Users/ahauser/Desktop/systematic_review_AH2/R/stan/TAM/fit_TAM.RData")

load(file="C:/Users/ahauser/Desktop/systematic_review_AH2/R/stan/TAM/fit_TAM.RData")
###############################################################################################################################################################################
###############################################################################################################################################################################
#Figure
load("C:/Users/ahauser/Desktop/systematic_review_AH2/R/stan/TAM/fit_NRTI.RData")
d<-summary(fit_NRTI)$summary

d1<-d[setdiff(grep("p_init",rownames(d)),grep("p_init_pred",rownames(d))),]
x_which=list(NULL,c("ddI","TDF","d4T"),NULL,NULL,"FTC.3TC",NULL,c("d4T","ZDV"),NULL)
unlist(sapply(x_which,function(x) c("Time",x)))
d1<-data.frame(Mutation=nrti_names,
               Median=c(d1[,"50%"],d4["p_init","50%"]),
               lower=c(d1[,"2.5%"],d4["p_init","2.5%"]),
               upper=c(d1[,"97.5%"],d4["p_init","97.5%"]))

#p_x
d2<-d[grep("p_x_FTC_ZDV",rownames(d)),]
x_which=list(NULL,c("ddI","TDF","d4T"),NULL,NULL,"FTC.3TC",NULL,c("d4T","ZDV"),NULL)
unlist(sapply(x_which,function(x) c("Time",x)))
d2<-data.frame(Mutation=nrti_names,
               Median=c(d2[,"50%"],d4["p_FTC_ZDV","50%"]),
               lower=c(d2[,"2.5%"],d4["p_FTC_ZDV","2.5%"]),
               upper=c(d2[,"97.5%"],d4["p_FTC_ZDV","97.5%"]))

d3<-d[grep("p_x_FTC_TDF",rownames(d)),]
#d2<-d[grep("p_xstat_pred",rownames(d)),]
x_which=list(NULL,c("ddI","TDF","d4T"),NULL,NULL,"FTC.3TC",NULL,c("d4T","ZDV"),NULL)
unlist(sapply(x_which,function(x) c("Time",x)))
d3<-data.frame(Mutation=nrti_names,
               Median=c(d3[,"50%"],d4["p_FTC_TDF","50%"]),
               lower=c(d3[,"2.5%"],d4["p_FTC_TDF","2.5%"]),
               upper=c(d3[,"97.5%"],d4["p_FTC_TDF","97.5%"]))


d<-rbind(d1,d1,d1,d2,d3)


d$reg=factor(c(rep("T1",9),
               rep("Baseline",9),
               rep("T2",9),
               rep("After 3 years on FTC/3TC+ZDV",9),
               rep("After 3 years on FTC/3TC+TDF",9)),
             levels=c("T1","Baseline","T2",
                      "After 3 years on FTC/3TC+TDF",
                      "After 3 years on FTC/3TC+ZDV"))
# d$Mutation=factor(rep(c("M41","K65","D67","K70","M184","L210","T215","K219"),5),
#                   levels=c("M41","K65","D67","K70","M184","L210","T215","K219"))
d$Mutation=factor(rep(c("M41","K65","D67","K70","M184","L210","T215","K219","TAM"),5),
                  levels=c("K65","M184","M41","D67","K70","L210","T215","K219","TAM"))
# d<-rbind(d[d$reg=="Baseline",],d[d$reg=="Baseline",],d)
# d$trans=rep(c(1,0,1,0,0),each=dim(d)[1]/5)


#study prevalences
n=n_nrti
k=k_nrti
p=k/matrix(rep(n,9),ncol=9)
study_est=data.frame(est=(as.vector(as.matrix(p))),
                     Mutation=factor(rep(c("M41","K65","D67","K70","M184","L210","T215","K219","TAM"),each=dim(p)[1]),
                                     levels=c("K65","M184","M41","D67","K70","L210","T215","K219","TAM")),
                     n_study=rep(n,9))
study_est=study_est[-which(study_est$Mutation=="TAM" & study_est$est==0),]

#plot
print(ggplot(d, aes(x=Mutation, y=Median),group=interaction(Mutation,color)) +
        geom_point(aes(color=reg),shape=16,size=3, position = position_dodge(width = 0.6)) +
        geom_errorbar(width=.1, size=1,aes(ymin=lower, ymax=upper,color=reg),
                      position = position_dodge(width = 0.6)) +
        labs(x="Mutation",y= "Prevalence (%)",hjust=0,
             title="Figure : Prevalence of 8 NRTI DRMs by first-line regimen.",
             subtitle="Black open circles: single studies estimates. Points and vertical lines: median and 95% credibility intervals of baseline prevalence (black),
prevalence after 3 years on 3TC/FTC + TDF (red) or 3TC/FTC + ZDV (blue). TAM: M41 D67, K70, L210, T215, K219.") +
        labs(color="Prevalence")+
        scale_y_continuous(limits=c(0,1),labels = function(x) paste0(x*100, "%"))+
        theme_bw()+
        geom_point(data=study_est,aes(y=est,size=sqrt(n_study)),shape=1,alpha=0.5)+
        labs(size = "Study size")+
        scale_size(limits = c(3.5, 30),breaks = sqrt(c(20,100,500)),labels=c(20,100,500))+
        scale_colour_manual(values=c("white","black","white","darkred","blue","black"),
                            breaks=c(levels(d$reg)[c(2,4,5)])))

print(ggplot(d, aes(x=Mutation, y=Median),group=interaction(Mutation,color)) +
        geom_point(aes(color=reg),shape=16,size=3, position = position_dodge(width = 0.6)) +
        geom_errorbar(width=.1, size=1,aes(ymin=lower, ymax=upper,color=reg),
                      position = position_dodge(width = 0.6)) +
        labs(x="Mutation",y= "Prevalence (%)",hjust=0) +
        labs(color="Prevalence")+
        scale_y_continuous(limits=c(0,1),labels = function(x) paste0(x*100, "%"))+
        theme_bw()+
        geom_point(data=study_est,aes(y=est,size=sqrt(n_study)),shape=1,alpha=0.5)+
        labs(size = "Study size")+
        scale_size(limits = c(3.5, 30),breaks = sqrt(c(20,100,500)),labels=c(20,100,500))+
        scale_colour_manual(values=c("white","black","white","darkred","blue","black"),
                            breaks=c(levels(d$reg)[c(2,4,5)])))
