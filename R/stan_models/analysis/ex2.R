N=20
M=2

n=seq(100,100+(N-1)*10,length.out = 20)
x=rnorm(N,2,2)/10
y=rnorm(N,4,2)/10
xy=cbind(x,y)
mu=-1
sigma=0.2
u= rnorm(N,0,1)
beta=c(2,2)
p=inv.logit(mu+sigma*u+xy%*%structure(beta,dim=c(M,1)))
k=rbinom(N,size=n,prob=p)


na_x_pos=c(2,5,8)
nna_x_pos=setdiff(1:N,na_x_pos)
x[na_x_pos]=0


list=list(N=N, #number of observations
          M_y=M-1, # number of y variables
          M_xy=M, # number of x and y variables
          x=structure(x,dim=N),
          y=structure(y,dim=c(N,1)),
          x_nna=x[nna_x_pos], #x in non-missing position
          na_x_pos=na_x_pos, #position of missing x
          nna_x_pos=nna_x_pos, #position of non-missing x
          k=k, #number of sucesses
          n=n, #number of trials
          inference=1)

model<-stan_model(file="C:/Users/ahauser/Desktop/systematic_review_AH2/R/stan/TAM/ex1.stan")
fit_ex<- sampling(object=model,
                    data = list,
                    warmup = 1000,
                    iter = 2500,
                    chains = 4,
                    cores = 4,
                    thin = 1,
                    control=list(adapt_delta=0.99,
                                 max_treedepth=10))
d = summary(fit_ex)$summary

