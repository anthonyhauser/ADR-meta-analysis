// Stan model for simple linear regression
data {
  int < lower = 1 > N; // number of studies
  int < lower = 1 > M; // number of mutations
  int < lower = 1 > B_tot; //number of covariates
  
  int k[N,M]; // Number of events by study
  int n[N]; // Number of people at risk by study
  
  matrix[N,1] y; //matrix of dim N times B, B co-variates for each of the N (sub-)study
  
  row_vector[B_tot] x_mean; //mean of x=(time,y)
  vector[B_tot] x_sd; //sd of x=(time,y)
  
  vector[B_tot] x_EFV; //prediction for EFV
  vector[B_tot] x_NVP; //prediction for NVP
  
  vector[N-1] time_nna; //ART time in non-missing position
  real mean_t; //mean of ART time
  real var_t; //variance of ART time
  int na_time_pos; //position of missing ART times
  int nna_time_pos[N-1]; //position of non-missing ART times
  
  real init_res_1[M];
  real init_res_2[M];
  
  real < lower = 0 > sigma_1[M];
  real < lower = 0 > tau_1;
  
  matrix[B_tot,M] beta_1;
  matrix[B_tot,M] beta_2;
  
  int inference;
}

parameters {
  real < lower = 0 > time_imput_na;
  
  matrix[B_tot,M] beta;
  real mu[M];
  
  real u[N];
  real v[N,M];
  vector < lower = 0 > [M] sigma;
  real < lower = 0 > tau;
}

transformed parameters {
  real alpha[N,M];
  vector[N] time_imput;
  vector[N] time_tr_imput;
  //vector[2] time_tr_imput_na;
  
  for(i in 1:N){
    for(j in 1:M){
      alpha[i,j] = mu[j] + sigma[j] * v[i,j] + tau * u[i];
    }
  }
  
  //Time imput
  time_imput[nna_time_pos] = time_nna;
  time_imput[na_time_pos] = time_imput_na;
  time_tr_imput = (time_imput - mean_t) ./ sqrt(var_t);
  
}

model {
  //Priors
  mu ~ normal(init_res_1,init_res_2);
  
  for(j in 1:M){
    sigma[j] ~ exponential(sigma_1[j]);
    for(i in 1:N){
      v[i,j] ~ normal(0,1);
    }
  }
  for(i in 1:B_tot){
    for(j in 1:M){
      beta[i,j] ~ normal(beta_1[i,j],beta_2[i,j]);
    }
  }
  
  //Priors time imputation
  time_imput_na ~ gamma(mean_t^2/var_t,mean_t/var_t);
  
  //Random intercept
  tau ~ exponential(tau_1);
  for(i in 1:N){
    u[i] ~ normal(0,1);
  }
  
  //Likelihood
  if(inference==1){
    for(i in 1:N){
      for(j in 1:M){
        target += binomial_lpmf(k[i,j]|n[i],inv_logit(alpha[i,j]+beta[1,j]*time_tr_imput[i]+beta[2,j]*y[i,1]));
      }
    }
  }
}


generated quantities {
  matrix[B_tot,M] beta_bt; //beta back-transformed
  
  vector[M] mu_bt;
  vector[M] p_init;
  vector[M] p_EFV;
  vector[M] p_NVP;
  vector[M] p_init_pred;
  vector[M] p_EFV_pred;
  vector[M] p_NVP_pred;
  
  real randnorm;
  
  // back transformed beta and mu
  for(i in 1:M){
    beta_bt[,i] = beta[,i] ./ x_sd;
  }
  
  for(i in 1:M){
    mu_bt[i] = mu[i] - x_mean * beta_bt[,i];
  }
  
  // credibility intervals
  p_init = inv_logit(mu_bt);
  // p_EFV = inv_logit(mu_bt + beta_bt' * x_EFV);
  // p_NVP = inv_logit(mu_bt + beta_bt' * x_NVP);
  for(j in 1:M){
    p_EFV[j] = inv_logit(mu[j]+beta[1,j]*(36.0-x_mean[1])/x_sd[1]+beta[2,j]*(1.0-x_mean[2])/x_sd[2]);
    p_NVP[j] = inv_logit(mu[j]+beta[1,j]*(36.0-x_mean[1])/x_sd[1]+beta[2,j]*(0.0-x_mean[2])/x_sd[2]);
  }
  
  // // predictive intervals
  randnorm = normal_rng(0,1);
  p_init_pred = inv_logit(mu_bt + (tau + sigma) * randnorm);
  p_EFV_pred = inv_logit(mu_bt + beta_bt' * x_EFV + (tau + sigma) * randnorm);
  p_NVP_pred = inv_logit(mu_bt + beta_bt' * x_NVP + (tau + sigma) * randnorm);
} // The prior or posterior predictive distribution

