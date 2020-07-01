// Stan model for simple linear regression
data {
  int < lower = 1 > N; // number of studies with TAMs
  int M; // Number of TAMs: M=6
  int k[N]; // Number of people with any TAMs by study
  int n[N]; // Number of people at risk by study
  
  vector[2] mu[N];
  matrix[2,2] cov [N];
  vector[2] mu_TDF;
  matrix[2,2] cov_TDF;
  vector[2] mu_ZDV;
  matrix[2,2] cov_ZDV;
  vector[2] mu_init;
  matrix[2,2] cov_init;
  
  int inference;
}

parameters {
  real < lower = 0, upper = 1 > alpha;
  ordered[2] inv_p[N];
  ordered[2] inv_p_TDF;
  ordered[2] inv_p_ZDV;
  ordered[2] inv_p_init;
}

transformed parameters {
  //P( any TAM ) = 1- P( no TAM )
  real p_TAM[N];
  
  for(i in 1:N){
    p_TAM[i] = 1 - ( inv_logit(inv_p[i,1]) + alpha * (inv_logit(inv_p[i,2]) - inv_logit(inv_p[i,1])));
  }
}

model {
  //Priors
  alpha  ~ uniform(0,1);
  for(i in 1:N){
    inv_p[i] ~ multi_normal(mu[i], cov[i]);
  }
  inv_p_TDF ~ multi_normal(mu_TDF, cov_TDF);
  inv_p_ZDV ~ multi_normal(mu_ZDV, cov_ZDV);
  inv_p_init ~ multi_normal(mu_init, cov_init);
  
  //Likelihood
  if(inference==1){
    for(i in 1:N){
      target += binomial_lpmf(k[i]|n[i], p_TAM[i]);
    }
  }
}

generated quantities {
  real p_FTC_ZDV;
  real p_FTC_TDF;
  real p_init;
  
  p_FTC_ZDV = 1 - ( inv_logit(inv_p_ZDV[1]) + alpha * (inv_logit(inv_p_ZDV[2]) -inv_logit(inv_p_ZDV[1])));
  p_FTC_TDF = 1 - ( inv_logit(inv_p_TDF[1]) + alpha * (inv_logit(inv_p_TDF[2]) - inv_logit(inv_p_TDF[1])));
  p_init = 1 - ( inv_logit(inv_p_init[1]) + alpha * (inv_logit(inv_p_init[2]) - inv_logit(inv_p_init[1])));
} // The prior or posterior predictive distribution
