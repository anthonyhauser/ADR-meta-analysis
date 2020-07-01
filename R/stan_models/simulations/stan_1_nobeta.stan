// Stan model for simple linear regression
data {
  int < lower = 1 > N; // number of studies
  int < lower = 1 > M; // number of mutations
  int < lower = 1 > B[M]; // number of covariates
  int < lower = 1 > B_tot;
  
  int k[N,M]; // Number of events by study
  int n[N]; // Number of people at risk by study
  
  real < lower = 0 > init_res_1[M];
  real < lower = 0 > init_res_2[M];
  real < lower = 0 > sigma_1[M];
  real < lower = 0 > tau_1;
  
  int inference;
}

parameters {
  real < lower = 0, upper = 1 > init_res[M];
  
  real u[N];
  real v[N,M];
  real < lower = 0 > sigma[M];
  real < lower = 0 > tau;
  
}

transformed parameters {
  real mu[M];
  real alpha[N,M];
  mu = logit(init_res);
  for(i in 1:N){
    for(j in 1:M){
      alpha[i,j] = mu[j] + sigma[j] * v[i,j] + tau * u[i];
    }
  }
}

model {
  //Priors
  
  init_res ~ beta(init_res_1,init_res_2);
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
      target += binomial_lpmf(k[i,1]|n[i],inv_logit(alpha[i,1]));
      target += binomial_lpmf(k[i,2]|n[i],inv_logit(alpha[i,2]));
      target += binomial_lpmf(k[i,3]|n[i],inv_logit(alpha[i,3]));
      target += binomial_lpmf(k[i,4]|n[i],inv_logit(alpha[i,4]));
      target += binomial_lpmf(k[i,5]|n[i],inv_logit(alpha[i,5]));
      target += binomial_lpmf(k[i,6]|n[i],inv_logit(alpha[i,6]));
      target += binomial_lpmf(k[i,7]|n[i],inv_logit(alpha[i,7]));
      target += binomial_lpmf(k[i,8]|n[i],inv_logit(alpha[i,8]));
    }
  }
}


generated quantities {
  real p_init[M];
  vector[M] p_init_pred;
  real randnorm;
  
  p_init = inv_logit(mu);
  
  randnorm = normal_rng(0,1);
  p_init_pred[1] = inv_logit(mu[1] + (tau + sigma[1]) * randnorm);
  p_init_pred[2] = inv_logit(mu[2] + (tau + sigma[2]) * randnorm);
  p_init_pred[3] = inv_logit(mu[3] + (tau + sigma[3]) * randnorm);
  p_init_pred[4] = inv_logit(mu[4] + (tau + sigma[4]) * randnorm);
  p_init_pred[5] = inv_logit(mu[5] + (tau + sigma[5]) * randnorm);
  p_init_pred[6] = inv_logit(mu[6] + (tau + sigma[6]) * randnorm);
  p_init_pred[7] = inv_logit(mu[7] + (tau + sigma[7]) * randnorm);
  p_init_pred[8] = inv_logit(mu[8] + (tau + sigma[8]) * randnorm);


} // The prior or posterior predictive distribution

