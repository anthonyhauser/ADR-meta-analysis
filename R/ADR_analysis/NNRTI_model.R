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

###############################################################################################################################################################################
###############################################################################################################################################################################

#Data management
load(file="results/sr.RData")
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

#-------------------------------------------------------------------------------------------------------------
#NNRTI mutations data
#retrieve nrti drugs
nnrti_pos=which(colnames(sr)=="K101_EHPQ"):which(colnames(sr)=="G190_AERS")
nnrti_names<-colnames(sr)[nnrti_pos]
x_nnrti=cbind(sr$Time_to_fail,sr$EFV,sr$NVP)
colnames(x_nnrti)<-c("Time_to_fail","EFV","NVP")

#remove studies with NA drugs
na_drug_x=apply(x_nnrti[,-1],1,function(x) sum(is.na(x))==0)
x_nnrti=x_nnrti[na_drug_x,]
sr_nnrti=sr[na_drug_x,]
print(sr_nnrti$Publication)

#save backtransformed data (time in months, drug in proportion)
x_nnrti[,-grep("Time_to_fail",colnames(x_nnrti))]=x_nnrti[,-grep("Time_to_fail",colnames(x_nnrti))]/100
x_nnrti_bt=x_nnrti

#k, n, x
n_nnrti=sr_nnrti$Genotyped
k_nnrti=as.data.frame(sr_nnrti[,nnrti_pos])
k_nnrti[is.na(k_nnrti)]=0 #NA is considered as 0
N_nnrti=length(n_nnrti)
M_nnrti=dim(k_nnrti)[2]

save(n_nnrti,k_nnrti,N_nnrti,M_nnrti,x_nnrti_bt,sr_nnrti,file="results/data_NNRTI.RData")
#-------------------------------------------------------------------------------------------------------------
#Select data and variables
#Data
k=k_nnrti
nnrti_names=colnames(k_nnrti)
N=N_nnrti #number of studies
M=M_nnrti #number of nnrti mutations
n=n_nnrti #number of genotyped by study

#Time
time_bt = x_nnrti_bt[,"Time_to_fail"] #time backtransformed (real time)
time_bt[sr_nnrti$Publication=="Labhardt 2016"]= 4.1 *12
na_time_pos = which(is.na(time_bt)) #position of NA values
nna_time_pos = which(!is.na(time_bt)) #position of non NA values
sr_nnrti[na_time_pos,"Publication"]
time_nna = time_bt[!is.na(time_bt)] #times with no NA
mean_t = mean(time_nna) #mean time
sd_t = sd(time_nna) #sd time
var_t = var(time_nna) #var time
#position of non NA values
time_tr=(time_bt-mean_t)/sd_t #transformed times
time_tr_nna=time_tr[nna_time_pos] #transformed times with no NA

#y: other variables: NNRTI drugs
y_nnrti_bt = x_nnrti_bt[,-1] #removing time
mean_y = apply(y_nnrti_bt,2,mean) #mean by NRTI drugs
sd_y = apply(y_nnrti_bt,2,sd) #sd by NRTI drugs
y_nnrti = (y_nnrti_bt - structure(rep(mean_y,each=dim(y_nnrti_bt)[1]),dim=dim(y_nnrti_bt)))/ #transformed NRTI drugs
  structure(rep(sd_y,each=dim(y_nnrti_bt)[1]),dim=dim(y_nnrti_bt))

#time_y: time binded with y
time_y_bt=cbind(time_bt,y_nnrti_bt) #backtransformed time_y
time_y_bt[na_time_pos,1]=0 #Replace NA by 0 as Stan does not accept NA values
time_y=cbind(time_tr,y_nnrti) #transformed time_y
time_y[na_time_pos,1]=0 #Replace NA by 0 as Stan does not accept NA values
x_mean=c(mean_t,mean_y) #Mean by variables
x_sd=c(sd_t,sd_y) #sd by variables

###############################################################################################################################################################################
###############################################################################################################################################################################
#Model the prevalence of NRTI mutations
mod_NNRTI<-stan_model(file="../stan_models/analysis/NNRTI_stan.stan")


x_EFV=c(36,1)
x_NVP=c(36,0)
B_tot=2

stan_data=list(N=N,
               M=M,
               B_tot=B_tot,
               k=matrix(round(as.matrix(k)),nrow=N),
               n=n,
               y=structure(time_y[,2],dim=c(N,1)),
               x_mean=structure(x_mean[-3],dim=B_tot),
               x_sd=structure(x_sd[-3],dim=B_tot),
               x_EFV=structure(x_EFV,dim=B_tot),
               x_NVP=structure(x_NVP,dim=B_tot),
               time_nna=time_nna, #ART time in non-missing position
               mean_t=mean_t, #mean of ART time
               var_t =var_t, #variance of ART time
               na_time_pos=na_time_pos, #position of missing ART times
               nna_time_pos=nna_time_pos, #position of non-missing ART times
               init_res_1=rep(0,M),
               init_res_2=rep(4,M),
               sigma_1=rep(1,M),
               tau_1=1,
               beta_1=matrix(rep(0,B_tot*M),nrow=B_tot),
               beta_2=matrix(rep(2,B_tot*M),nrow=B_tot),
               inference=1)

fit_NNRTI <-sampling(object=mod_NNRTI,
             data = stan_data,
             warmup = 1000,
             iter = 2500,
             chains = 4,
             cores = 4,
             thin = 1,
             control=list(max_treedepth=15,adapt_delta=0.99))

save(fit_NNRTI,file="results/fit_NNRTI.RData")
d<-summary(fit_NNRTI)$summary

####################################################################################################################
#first look at the results
d[grep("beta",rownames(d)),]
d[grep("p_EFV",rownames(d)),]
d[grep("p_NVP",rownames(d)),]
d[grep("p_init",rownames(d)),]
d[grep("time_imput",rownames(d)),]
apply(k,2,function(x) sum(x)/sum(n))
(logit((d[grep("p_EFV",rownames(d)),"50%"])[1:7])-logit((d[grep("p_NVP",rownames(d)),"50%"])[1:7]))/36
d[grep("beta",rownames(d)),]
plot(time_y_bt[,2],k[,1]/n)

####################################################################################################################
#Figure
load("results/fit_NNRTI.RData")
d<-summary(fit_NNRTI)$summary
d1<-d[setdiff(grep("p_EFV",rownames(d)),grep("p_EFV_pred",rownames(d))),]
d1<-data.frame(Mutation=nnrti_names,
               Median=c(d1[,"50%"]),
               lower=c(d1[,"2.5%"]),
               upper=c(d1[,"97.5%"]))
d2<-d[setdiff(grep("p_NVP",rownames(d)),grep("p_NVP_pred",rownames(d))),]
d2<-data.frame(Mutation=nnrti_names,
               Median=c(d2[,"50%"]),
               lower=c(d2[,"2.5%"]),
               upper=c(d2[,"97.5%"]))
d<-rbind(d1,d2)

d$reg=factor(c(rep("After 3 years on EFV",7),
               rep("After 3 years on NVP",7)),
             levels=c("After 3 years on EFV","After 3 years on NVP"))
d$Mutation=factor(rep(c("K101","K103","V106","V108","Y181","Y188","G190"),2),
                  levels=c("K101","K103","V106","V108","Y181","Y188","G190"))

#study prevalences
n=n_nnrti
k=k_nnrti
p=k/matrix(rep(n,7),ncol=7)
study_est=data.frame(est=(as.vector(as.matrix(p))),
                     Mutation=factor(rep(c("K101","K103","V106","V108","Y181","Y188","G190"),each=dim(p)[1]),
                                     levels=c("K101","K103","V106","V108","Y181","Y188","G190")),
                     n_study=rep(n,7))

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
        scale_colour_manual(values=c("darkred","blue"),
                            breaks=c(levels(d$reg)[c(1,2)])))

print(ggplot(d, aes(x=Mutation, y=Median),group=interaction(Mutation,color)) +
        geom_point(aes(color=reg),shape=16,size=3, position = position_dodge(width = 0.6)) +
        geom_errorbar(width=.1, size=1,aes(ymin=lower, ymax=upper,color=reg),
                      position = position_dodge(width = 0.6)) +
        labs(x="Mutation",y= "Prevalence (%)",hjust=0)+
             labs(color="Prevalence")+
        scale_y_continuous(limits=c(0,1),labels = function(x) paste0(x*100, "%"))+
        theme_bw()+
        geom_point(data=study_est,aes(y=est,size=sqrt(n_study)),shape=1,alpha=0.5)+
        labs(size = "Study size")+
        scale_size(limits = c(3.5, 30),breaks = sqrt(c(20,100,500)),labels=c(20,100,500))+
        scale_colour_manual(values=c("darkred","blue"),
                            breaks=c(levels(d$reg)[c(1,2)])))
