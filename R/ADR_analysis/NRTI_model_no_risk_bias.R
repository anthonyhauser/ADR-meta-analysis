#Model_estimation
#selection of variables
#imputation of time on ART missing values, stepwise selection

#########################################################################################################################
library(rstan) #https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
library(gtools)
library(gdata)
library(bayesplot)
library("R.utils")
source('../stan_models/stan_utility.R')
#Stan configuration
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load("results/data_NRTI.RData")
which_no_risk_bias_nrti= sr_nrti[,"risk_of_bias"]==0
###############################################################################################################################################################################
###############################################################################################################################################################################
#Sensitivity analysis: remove study with high risk of bias
#Data
k=k_nrti[which_no_risk_bias_nrti,-which(colnames(k_nrti)=="TAM")] #remove TAMs for now
nrti_names=colnames(k_nrti)
N=sum(which_no_risk_bias_nrti) #number of studies
M=M_nrti-1 #number of nrti mutations, -1 as we removed TAMs
n=n_nrti[which_no_risk_bias_nrti] #number of genotyped by study

#Time
time_bt = x_nrti_bt[which_no_risk_bias_nrti,"Time_to_fail"] #time backtransformed (real time)
time_bt[sr_nrti[which_no_risk_bias_nrti,"Publication"]=="Labhardt 2016"]= 4.1 *12
na_time_pos = which(is.na(time_bt)) #position of NA values
nna_time_pos = which(!is.na(time_bt)) #position of non NA values
print((sr_nrti[which_no_risk_bias_nrti,"Publication"])[na_time_pos])
time_nna = time_bt[nna_time_pos] #times with no NA
mean_t = mean(time_nna) #mean time
sd_t = sd(time_nna) #sd time
var_t = var(time_nna) #var time
time_tr=(time_bt-mean_t)/sd_t #transformed times
time_tr_nna=time_tr[nna_time_pos] #transformed times with no NA

#y: other variables: NRTI drugs
y_nrti_bt = x_nrti_bt[which_no_risk_bias_nrti,-1] #removing time
mean_y = apply(y_nrti_bt,2,mean) #mean by NRTI drugs
sd_y = apply(y_nrti_bt,2,sd) #sd by NRTI drugs
y_nrti = (y_nrti_bt - structure(rep(mean_y,each=dim(y_nrti_bt)[1]),dim=dim(y_nrti_bt)))/ #transformed NRTI drugs
  structure(rep(sd_y,each=dim(y_nrti_bt)[1]),dim=dim(y_nrti_bt))

#time_y: time binded with y
time_y_bt=cbind(time_bt,y_nrti_bt) #backtransformed time_y
time_y_bt[na_time_pos,1]=0 #Replace NA by 0 as Stan does not accept NA values
time_y=cbind(time_tr,y_nrti) #transformed time_y
time_y[na_time_pos,1]=0 #Replace NA by 0 as Stan does not accept NA values
x_mean=c(mean_t,mean_y) #Mean by variables
x_sd=c(sd_t,sd_y) #sd by variables

#-------------------------------------------------------------------------------------------------------------
#Model the prevalence of NRTI mutations
mod_NRTI<-stan_model(file="../stan_models/analysis/NRTI_stan2_no_risk_bias.stan") #calculate p_x_study

#Starting setting
load("results/fit_NRTI.RData")
B_tot=dim(time_y_bt)[2]
B_list_assumed=B_list_assumed_NRTI
B=sapply(B_list_assumed,function(x) length(x))

y_list<-lapply(B_list_assumed,function(x){
  return(structure(time_y[,x],dim=c(dim(time_y)[1],length(x))))
})
y_bt_list<-lapply(B_list_assumed,function(x){
  return(structure(time_y_bt[,x],dim=c(dim(time_y_bt)[1],length(x))))
})

for(i in 1:M){
  assign(paste("x_m",i,sep=""),y_list[[i]])
  assign(paste("x_study",i,sep=""),y_bt_list[[i]])
  if(B[i]>1) prov=c(36,as.numeric((B_list_assumed[[i]])[-1] %in% c(2,5))) else prov=36
  assign(paste("x_FTC_TDF",i,sep=""),prov)
  if(B[i]>1) prov=c(36,as.numeric((B_list_assumed[[i]])[-1] %in% c(2,6))) else prov=36
  assign(paste("x_FTC_ZDV",i,sep=""),prov)
}

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


fit_NRTI_no_risk_bias<- sampling(object=mod_NRTI,
                                 data = stan_NRTI,
                                 warmup = 1000,
                                 iter = 2500,
                                 chains = 4,
                                 cores = 4,
                                 thin = 1,
                                 control=list(adapt_delta=0.99,
                                              max_treedepth=15))


d<-summary(fit_NRTI_no_risk_bias)$summary
save(fit_NRTI_no_risk_bias,file="results/fit_NRTI_no_risk_bias.RData")
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
save(fit_TAM,file="results/fit_TAM_no_risk_bias.RData")

###############################################################################################################################################################################
###############################################################################################################################################################################
#Figure
load(file="results/fit_TAM_no_risk_bias.RData")
d4<-summary(fit_TAM)$summary
load("results/fit_NRTI_no_risk_bias.RData")
d<-summary(fit_NRTI_no_risk_bias)$summary

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
n=n_nrti[which_no_risk_bias_nrti] 
k=k_nrti[which_no_risk_bias_nrti,]
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
