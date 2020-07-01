// Stan model for simple linear regression
data {
  int < lower = 1 > N; // number of studies
  int < lower = 1 > B; // number of covariates
  
  int k[N]; // Number of events by study
  int n[N]; // Number of people at risk by study
  
  matrix[N,B] x; //matrix of dim N times B, B co-variates for each of the N (sub-)study
  
  real < lower = 0 > init_res_1;
  real < lower = 0 > init_res_2;
  real < lower = 0 > sigma_1;
  
  int inference;
}

parameters {
  real < lower = 0, upper = 1 > init_res; 
  real < lower = 0 > sigma;
  real alpha_tilde[N];
}

transformed parameters {
  real mu = logit(init_res);
  real alpha[N];
  for(i in 1:N){
     alpha[i] = mu + sigma * alpha_tilde[i];
  }
}

model {
  //Priors
  init_res ~ beta(init_res_1,init_res_2);
  sigma ~ exponential(sigma_1);
  //sigma ~ uniform(1/sigma_1,sigma_1);
  //Random intercept
  for(i in 1:N){
    alpha_tilde[i] ~ normal(0,1);
  }
  
  //Likelihood
  if(inference==1){
    for(i in 1:N){
      target += binomial_lpmf(k[i]|n[i],inv_logit(alpha[i]));
    }
  }
}

generated quantities {
  real p_pred_min = inv_logit(normal_rng(mu,sigma));
  real p_pred_max = inv_logit(normal_rng(mu,sigma));
} // The prior or posterior predictive distribution

