// Stan model for simple linear regression
data {
  int < lower = 1 > N; // number of studies
  int < lower = 1 > B; // number of co-variates
  
  int k[N]; // Number of events by study
  int n[N]; // Number of people at risk by study
  
  matrix[N,B] x; //matrix of dim N times B, B co-variates for each of the N (sub-)study
  
  real < lower = 0 > init_res_1;
  real < lower = 0 > init_res_2;
  real < lower = 0 > sigma_1;
  real < lower = 0 > beta_1;
  //vector < lower = 0 > [B] beta_2;
  //real <lower = 0 > beta_2;
  int inference;
}

parameters {
  real < lower = 0, upper = 1 > init_res; 
  real < lower = 0 > sigma;
  real alpha_tilde [N];
  //vector < lower = 0 > [B] beta;
  vector < lower = 0 > [B] beta_2;
  vector [B] beta_tilde;
  //real <lower = 0 > beta_2;
}

transformed parameters {
  real mu = logit(init_res);
  real alpha[N];
  vector[B] beta;
  
  for(i in 1:N){
    alpha[i] = mu + sigma * alpha_tilde[i];
  }
  
  for(i in 1:B){
    beta[i] = beta_1 * beta_2[i] * beta_tilde[i];
  }
}

model {
  //Priors
  for(i in 1:B){
      //beta_2[i] ~ cauchy(0,1);
      beta_2[i] ~ cauchy(0,1);
  }
  init_res ~ beta(init_res_1,init_res_2);
  sigma ~ exponential(sigma_1);
  //b ~ exponential(1);
  //sigma ~ cauchy(0,sigma_1);
  //beta ~ normal(beta_1,beta_2);
 
  for(i in 1:B){
  beta_tilde[i] ~ normal(0,1);
  }
  
  
  //Random intercept
  for(i in 1:N){
    alpha_tilde[i] ~ normal(0,1);
  }
  
  //Likelihood
  if(inference==1){
    for(i in 1:N){
      target += binomial_lpmf(k[i]|n[i],inv_logit(alpha[i]+x[i,]*beta/100));
    }
  }
}

generated quantities {
  // real p_pred_min = inv_logit(normal_rng(mu,sigma));
  //real p_pred_max = inv_logit(normal_rng(mu+beta[1]*50/500+beta[2]*100/500,sigma));
} // The prior or posterior predictive distribution

