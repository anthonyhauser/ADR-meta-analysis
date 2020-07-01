# Analysing the Systematic Review data

# remove the 'all authors', 'last author', 'cohort', 'title' and 'comments' columns before saving it as a csv!
# sr.dat <- read.table("C:/Users/ahauser/Documents/SystematicReview_Fardo/systematic_review_AH/Excel_data/ResHIV_dataExtract_AH.csv", header=TRUE, sep=";", dec=".", as.is = c(1,5,8,9,16,23))
#sr.dat <- read.table("../../data/ResHIV_dataExtract_withTAM_AH.csv", header=TRUE, sep=",", dec=".", as.is = c(1,5,8,9,16,23))
sr.dat <- read.csv("../../data/ResHIV_dataExtract_withTAM_AH.csv", header=TRUE, sep=",", dec=".", as.is = c(1,5,8,9,16,23))

sr<-sr.dat
sr<-sr[,-(which(colnames(sr)=="d4T"):which(colnames(sr)=="NFV"))]
colmutation <-seq(from = which(names(sr)=="M41_L"), to=which(names(sr)=="L98_IM"))

#keep only important mutations
sr<-sr[,!(names(sr) %in% setdiff(names(sr)[colmutation],
                                 c(c("M184_IV","D67_NGE","K65_R","K70_EGR","T215_FINSY","K219_DENQR","M41_L","L210_W"),
                                   c("K103_NRS","V106_AIM","G190_AERS","Y181_CIV","V108_I","Y188_CHL","K101_EHPQ"))))]
#########################################################################################################################
#########################################################################################################################
#Fardo and Anthony part: cleaning data
# Identify study year: either the end of enrollment, or if not known, end of followup.
surv.year <- sr$enroll_end
surv.year[is.na(surv.year)]<-sr$followup_end[is.na(surv.year)]

# Adding the publication column (and surv.year)
sr<-data.frame(Publication = paste(sr$First_author, sr$Year_published), sr, surv.year, stringsAsFactors=FALSE)

#remove El-Kathib 2010 (line 17), Dlamini, Steegen 2016, Doualla-Bell 2006, Levison, as it is second-line PI treatment
#remove El-Khatib (line 31) it may be duplicate
#remove Doualla Bell 2009, as second-line
sr<-sr[-c(3,8,17,18,24,26,31),]

# identifies publications with multiple entries (2 pops)
double_int<-which(sr$Publication %in% names(table(sr$Publication))[table(sr$Publication)>=2] )
sr$sub_n=0
double_pos<-lapply(as.list(unique(sr$Publication[double_int])),function(x) which(sr$Publication==x))
for(i in 1:length(double_pos)){
  sr$Publication[double_pos[[i]]]<-paste(sr$Publication[(double_pos[[i]])[1]],c("a","b","c","d")[1:length(double_pos[[i]])],sep="")
  sr$sub_n[double_pos[[i]]]<-i
}

# Data cleaning
sr$study_type[sr$study_type=="survey (prestudy)"]<-"survey" 
author<-sr$First_author
#########################################################################################################################
#########################################################################################################################
#Anthony part:proportion of drug
#Create art.table (initial treatment) and art.table2 (treatment at genotype)
art.table <- sr[,which(colnames(sr)=="ABC.ddI"):which(colnames(sr)=="PIdrug")]
art.table <- cbind(art.table,rep(NA,dim(art.table)[1]))
art.table2 <- sr[,which(colnames(sr)=="ABC.ddI.2"):which(colnames(sr)=="PIdrug.2")]
#------------------------------------------------------------------------------------------------------
#Column with nrti and nnrti treatment
c_nrti<-which(colnames(art.table2)=="ABC.ddI.2"):which(colnames(art.table2)=="d4T.ddI.2")
c_nnrti<-which(colnames(art.table2)=="EFV.2"):which(colnames(art.table2)=="PIdrug.2")
#------------------------------------------------------------------------------------------------------
#When row is all NA in table2 (treatment at genotype), replace by table (initial treatment)
#nrti
art.table2[apply(art.table2[,c_nrti],1,function(x) sum(!is.na(x))==0),c_nrti]<-
  art.table[apply(art.table2[,c_nrti],1,function(x) sum(!is.na(x))==0),c_nrti]
#nnrti
art.table2[apply(art.table2[,c_nnrti],1,function(x) sum(!is.na(x))==0),c_nnrti]<-
  art.table[apply(art.table2[,c_nnrti],1,function(x) sum(!is.na(x))==0),c_nnrti]
#------------------------------------------------------------------------------------------------------
#Replace NA by 0 when sum of nnrti or nrti proportion is 100 (or >99, due to rounding)
#Warning: do not do this for xxx studies as some people have 4 drugs and 3 nrtis
#nrti
apply(art.table2[,c_nrti],1,function(x) sum(x[!is.na(x)]))
apply(art.table2[,c_nnrti],1,function(x) sum(x[!is.na(x)]))
art.table2[,c_nrti]<-t(apply(art.table2[,c_nrti],1,function(x)
{r=x
if(sum(x[!is.na(x)])>99){
  r[is.na(x)]=0
}
return(r)}))
#nnrti
art.table2[,c_nnrti]<-t(apply(art.table2[,c_nnrti],1,function(x)
{r=x
if(sum(x[!is.na(x)])>99){
  r[is.na(x)]=0
}
return(r)}))

#number of studies with no art proportion, with complete nrti and complete nnrti proportion
sum(apply(art.table2,1,function(x) sum(!is.na(x))==0))
sum(apply(art.table2[,c_nrti],1,function(x) sum(is.na(x))==0))
sum(apply(art.table2[,c_nnrti],1,function(x) sum(is.na(x))==0))
#------------------------------------------------------------------------------------------------------
#Do a new table with single drug
drug<-unlist(strsplit(colnames(art.table2), "\\."))
drug<-unique(drug[-which(drug=="2")])
drug.table<-matrix(NA,dim(art.table2)[1],length(drug))
colnames(drug.table)=drug
for(i in 1:dim(drug.table)[1]){
  for(j in 1:dim(drug.table)[2]){
    drug_j=drug[j]
    drug.table[i,j]=sum(art.table2[i,grep(drug_j,colnames(art.table2))],na.rm=TRUE) #already consider here that NA=0
  }
}
#Small correction: Singh study has participants with 3 NRTIs (ddI or 3TC with AZT and d4T), ddI and 3TC has been counted twice for these patients, remove it
drug.table[which(author=="Sing"),"ddI"]=drug.table[which(author=="Sing"),"ddI"]-2.2
drug.table[which(author=="Sing"),"3TC"]=drug.table[which(author=="Sing"),"3TC"]-4.4

#Result: percentage of drug used by class
nrti_p <- apply(drug.table[, colnames(drug.table) %in% c("ddI","d4T","ZDV","TDF","3TC","FTC")],1,function(x) sum(x,na.rm=TRUE))
nnrti_p<- apply(drug.table[, colnames(drug.table) %in% c("EFV","NVP")],1,function(x) sum(x,na.rm=TRUE))

#run, but not useful anymore, as we already considered above that NA is 0
for(i in 1:dim(drug.table)[1]){
  x=drug.table[i,1:7]
  if(sum(x)<150){
    x[x==0]=NA
    drug.table[i,1:7]=x
  }
  x=drug.table[i,8:9]
  if(sum(x)<90){
    x[x==0]=NA
    drug.table[i,8:9]=x
  }}

#------------------------------------------------------------------------------------------------------
#combine sr and drug.table tables
sr<-sr[,c(1:(which(colnames(sr)=="ABC.ddI")-1),which(colnames(sr)=="risk_of_bias"))] #remove column of nnrti combination
sr<-cbind(sr,drug.table) #add nrti nnrti and pi drugs
sr$'FTC.3TC'=sr$`3TC`+sr$FTC
sr$`FTC.3TC`[sr$`FTC.3TC`>100]=100

save(sr,file="results/sr.RData")