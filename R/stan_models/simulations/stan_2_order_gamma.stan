// Stan model for simple linear regression
functions {
  real min_2real(real x, real y) {
    vector[2] v;
    v[1]=x;
    v[2]=y;
    return min(v);
  }
  real normal_mixture_lpdf(real y, real s, real mu1, real mu2, real sigma1, real sigma2) {
    return log_sum_exp(log(s) + gamma_lpdf(y | mu1, sigma1),
                       log(1-s) + normal_lpdf(y | mu2, sigma2));
  }
}
data {
  int < lower = 1 > N; // number of studies
  int < lower = 1 > M; // number of mutations
  int < lower = 1 > B[M]; // number of covariates
  int < lower = 1 > B_tot;
  
  int k[N,M]; // Number of events by study
  int n[N]; // Number of people at risk by study
  
  matrix[N,B_tot] x; //matrix of dim N times B, B co-variates for each of the N (sub-)study
  
  row_vector[B_tot] x_mean;
  vector[B_tot] x_sd;
  
  vector[B_tot] x_pred;
  vector[B_tot] x_stat;
  
  real init_res_1[M];
  real init_res_2[M];
  real init_add[M];
  real < lower = 0 > sigma_1[M];
  real < lower = 0 > tau_1;
  
  real beta_p1;
  real beta_p2;
  vector[2] beta_time;
  vector[2] beta_3l;
  vector[2] beta_4l;
  
  real epsilon;
  real prob_bernoulli;
  
  int inference;
}

parameters {
  
  real <lower = 0> beta_11;
  real <lower = 0> beta_12;
  real <lower = 0> beta_13;
  real <lower = 0> beta_14;
  real <lower = 0> beta_15;
  real <lower = 0> beta_16;
  
  real <lower = 0> beta_21;
  real <lower = 0> beta_22;
  real <lower = 0> beta_26;
  real <lower = 0> beta_23;
  real <lower = 0> beta_24;
  real <lower = 0> beta_25;
  
  real <lower = 0> beta_31;
  real <lower = 0> beta_34;
  real <lower = 0> beta_32;
  real <lower = 0> beta_33; //<lower = beta_32>
  real <lower = 0> beta_35; //<lower = beta_32>
  real <lower = 0> beta_36;
    
    real <lower = 0>beta_41;
    real <lower = 0> beta_43;
    real <lower = beta_43> beta_42;
    real <lower = 0> beta_45;
    real <lower = 0> beta_44;
    real <lower = 0> beta_46;
    
    real <lower = 0> beta_51;
    real <lower = 0> beta_53;
    real <lower = beta_53> beta_52;
    real <lower = 0> beta_55;
    real <lower = 0> beta_54;
    real <lower = 0> beta_56;
    
    real <lower = 0> beta_61;
    real <lower = 0> beta_63;
    real <lower = beta_63> beta_62;
    real <lower = 0> beta_64;
    real <lower = 0, upper = beta_64> beta_65; //<lower = beta_62>
    real <lower = 0, upper = beta_64> beta_66;
      
    real <lower = 0> beta_71;
    real <lower = 0> beta_72;
    real <lower = 0> beta_73;
    real <lower = 0> beta_76;
    real <lower=0,upper = beta_76> beta_74;
    real <lower=0,upper = beta_76> beta_75;
        
    real <lower = 0> beta_81;
    real <lower = 0> beta_83;
    real <lower = beta_83> beta_82;
    real <lower = 0> beta_86;
    real <lower = beta_86> beta_84;
    real <lower = beta_86> beta_85;
        
    //real < lower = 0, upper = 1 > init_res[M];
        
    real u[N];
    real v[N,M];
    vector < lower = 0 > [M] sigma;
    real < lower = 0 > tau;
    real mu[M];
        
        
}

transformed parameters {
  //real mu[M];
  real alpha[N,M];
  //mu = logit(init_res);
  // for(i in 1:M){
  //   mu[i] = logit(init_res[i]) + init_add[i];
  // }
  for(i in 1:N){
    for(j in 1:M){
      alpha[i,j] = mu[j] + sigma[j] * v[i,j] + tau * u[i];
    }
  }
  
}

model {
  //Priors
  beta_11 ~ gamma(beta_p1,beta_p2);
  beta_12 ~ gamma(beta_p1,beta_p2);
  beta_13 ~ gamma(beta_p1,beta_p2);
  beta_14 ~ gamma(beta_p1,beta_p2);
  beta_15 ~ gamma(beta_p1,beta_p2);
  beta_16 ~ gamma(beta_p1,beta_p2);
  
  beta_21 ~ gamma(beta_p1,beta_p2);
  beta_22 ~ gamma(beta_p1,beta_p2);
  beta_23 ~ gamma(beta_p1,beta_p2);
  beta_24 ~ gamma(beta_p1,beta_p2);
  beta_25 ~ gamma(beta_p1,beta_p2);
  beta_26 ~ gamma(beta_p1,beta_p2);
  
  beta_31 ~ gamma(beta_p1,beta_p2);
  beta_32 ~ gamma(beta_p1,beta_p2);
  beta_33 ~ gamma(beta_p1,beta_p2);
  beta_34 ~ gamma(beta_p1,beta_p2);
  beta_35 ~ gamma(beta_p1,beta_p2);
  beta_36 ~ gamma(beta_p1,beta_p2);
  
  beta_41 ~ gamma(beta_p1,beta_p2);
  beta_42 ~ gamma(beta_p1,beta_p2);
  beta_43 ~ gamma(beta_p1,beta_p2);
  beta_44 ~ gamma(beta_p1,beta_p2);
  beta_45 ~ gamma(beta_p1,beta_p2);
  beta_46 ~ gamma(beta_p1,beta_p2);
  
  beta_51 ~ gamma(beta_p1,beta_p2);
  beta_52 ~ gamma(beta_p1,beta_p2);
  beta_53 ~ gamma(beta_p1,beta_p2);
  beta_54 ~ gamma(beta_p1,beta_p2);
  beta_55 ~ gamma(beta_p1,beta_p2);
  beta_56 ~ gamma(beta_p1,beta_p2);
  
  beta_61 ~ gamma(beta_p1,beta_p2);
  beta_62 ~ gamma(beta_p1,beta_p2);
  beta_63 ~ gamma(beta_p1,beta_p2);
  beta_64 ~ gamma(beta_p1,beta_p2);
  beta_65 ~ gamma(beta_p1,beta_p2);
  beta_66 ~ gamma(beta_p1,beta_p2);
  
  beta_71 ~ gamma(beta_p1,beta_p2);
  beta_72 ~ gamma(beta_p1,beta_p2);
  beta_73 ~ gamma(beta_p1,beta_p2);
  beta_74 ~ gamma(beta_p1,beta_p2);
  beta_75 ~ gamma(beta_p1,beta_p2);
  beta_76 ~ gamma(beta_p1,beta_p2);
  
  beta_81 ~ gamma(beta_p1,beta_p2);
  beta_82 ~ gamma(beta_p1,beta_p2);
  beta_83 ~ gamma(beta_p1,beta_p2);
  beta_84 ~ gamma(beta_p1,beta_p2);
  beta_85 ~ gamma(beta_p1,beta_p2);
  beta_86 ~ gamma(beta_p1,beta_p2);
  
  mu ~ normal(init_res_1,init_res_2);
  //init_res ~ beta(init_res_1,init_res_2);
  for(j in 1:M){
    sigma[j] ~ exponential(sigma_1[j]);
    for(i in 1:N){
      v[i,j] ~ normal(0,1);
    }
  }
  
  //Random intercept
  tau ~ exponential(tau_1);
  for(i in 1:N){
    u[i] ~ normal(0,1);
  }
  
  //Likelihood
  if(inference==1){
    for(i in 1:N){
      target += binomial_lpmf(k[i,1]|n[i],inv_logit(alpha[i,1]+beta_11*x[i,1]+beta_12*x[i,2]+beta_13*x[i,3]+beta_14*x[i,4]+beta_15*x[i,5]+beta_16*x[i,6]));
      target += binomial_lpmf(k[i,2]|n[i],inv_logit(alpha[i,2]+beta_21*x[i,1]+beta_22*x[i,2]+beta_23*x[i,3]+beta_24*x[i,4]+beta_25*x[i,5]+beta_26*x[i,6]));
      target += binomial_lpmf(k[i,3]|n[i],inv_logit(alpha[i,3]+beta_31*x[i,1]+beta_32*x[i,2]+beta_33*x[i,3]+beta_34*x[i,4]+beta_35*x[i,5]+beta_36*x[i,6]));
      target += binomial_lpmf(k[i,4]|n[i],inv_logit(alpha[i,4]+beta_41*x[i,1]+beta_42*x[i,2]+beta_43*x[i,3]+beta_44*x[i,4]+beta_45*x[i,5]+beta_46*x[i,6]));
      target += binomial_lpmf(k[i,5]|n[i],inv_logit(alpha[i,5]+beta_51*x[i,1]+beta_52*x[i,2]+beta_53*x[i,3]+beta_54*x[i,4]+beta_55*x[i,5]+beta_56*x[i,6]));
      target += binomial_lpmf(k[i,6]|n[i],inv_logit(alpha[i,6]+beta_61*x[i,1]+beta_62*x[i,2]+beta_63*x[i,3]+beta_64*x[i,4]+beta_65*x[i,5]+beta_66*x[i,6]));
      target += binomial_lpmf(k[i,7]|n[i],inv_logit(alpha[i,7]+beta_71*x[i,1]+beta_72*x[i,2]+beta_73*x[i,3]+beta_74*x[i,4]+beta_75*x[i,5]+beta_76*x[i,6]));
      target += binomial_lpmf(k[i,8]|n[i],inv_logit(alpha[i,8]+beta_81*x[i,1]+beta_82*x[i,2]+beta_83*x[i,3]+beta_84*x[i,4]+beta_85*x[i,5]+beta_86*x[i,6]));
    }
  }
}


generated quantities {
  matrix[B_tot,M] beta; //beta
  matrix[B_tot,M] beta_bt;
  
  vector[M] mu_bt;
  vector[M] p_init;
  vector[M] p_x;
  vector[M] p_xstat;
  
  beta[1,1]=beta_11;
  beta[1,2]=beta_21;
  beta[1,3]=beta_31;
  beta[1,4]=beta_41;
  beta[1,5]=beta_51;
  beta[1,6]=beta_61;
  beta[1,7]=beta_71;
  beta[1,8]=beta_81;
  
  beta[2,1]=beta_12;
  beta[2,2]=beta_22;
  beta[2,3]=beta_32;
  beta[2,4]=beta_42;
  beta[2,5]=beta_52;
  beta[2,6]=beta_62;
  beta[2,7]=beta_72;
  beta[2,8]=beta_82;
  
  beta[3,1]=beta_13;
  beta[3,2]=beta_23;
  beta[3,3]=beta_33;
  beta[3,4]=beta_43;
  beta[3,5]=beta_53;
  beta[3,6]=beta_63;
  beta[3,7]=beta_73;
  beta[3,8]=beta_83;
  
  beta[4,1]=beta_14;
  beta[4,2]=beta_24;
  beta[4,3]=beta_34;
  beta[4,4]=beta_44;
  beta[4,5]=beta_54;
  beta[4,6]=beta_64;
  beta[4,7]=beta_74;
  beta[4,8]=beta_84;
  
  beta[5,1]=beta_15;
  beta[5,2]=beta_25;
  beta[5,3]=beta_35;
  beta[5,4]=beta_45;
  beta[5,5]=beta_55;
  beta[5,6]=beta_65;
  beta[5,7]=beta_75;
  beta[5,8]=beta_85;
  
  beta[6,1]=beta_16;
  beta[6,2]=beta_26;
  beta[6,3]=beta_36;
  beta[6,4]=beta_46;
  beta[6,5]=beta_56;
  beta[6,6]=beta_66;
  beta[6,7]=beta_76;
  beta[6,8]=beta_86;
  
  
  for(i in 1:M){
    beta_bt[,i] = beta[,i] ./ x_sd;
  }
  
  
  for(i in 1:M){
    mu_bt[i] = mu[i] - x_mean * beta_bt[,i];
  }
  
  // credibility intervals
  p_init = inv_logit(mu_bt);
  p_x = inv_logit(mu_bt + beta_bt' *x_pred);
  p_xstat = inv_logit(mu_bt + beta_bt' *x_stat);
  
  
  
} // The prior or posterior predictive distribution
