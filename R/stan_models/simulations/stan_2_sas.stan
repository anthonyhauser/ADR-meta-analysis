// Stan model for simple linear regression
functions {
  real normal_mixture_lpdf(real y, real s, real mu1, real mu2, real sigma1, real sigma2) {
    return log_sum_exp(log(s) + normal_lpdf(y | mu1, sigma1),
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
  
  
  real c;
  real epsilon;
  real prob_bernoulli;
  matrix[B_tot,M] mean_c;
  
  int inference;
}

parameters {
  matrix[B_tot,M] beta;
  //real < lower = 0, upper = 1 > init_res[M];
  real mu[M];
  
  real u[N];
  real v[N,M];
  vector < lower = 0 > [M] sigma;
  real < lower = 0 > tau;
  
}

transformed parameters {
  //real mu[M];
  real alpha[N,M];
  matrix[N,M] bx;
  //mu = logit(init_res);
  // for(i in 1:M){
  //   mu[i] = logit(init_res[i]) + init_add[i];
  // }
  for(i in 1:N){
    for(j in 1:M){
      alpha[i,j] = mu[j] + sigma[j] * v[i,j] + tau * u[i];
    }
  }
  bx = x * beta;
}

model {
  //Priors
  mu ~ normal(init_res_1,init_res_2);
  //init_res ~ beta(init_res_1,init_res_2);
  for(j in 1:M){
    sigma[j] ~ exponential(sigma_1[j]);
    for(i in 1:N){
      v[i,j] ~ normal(0,1);
    }
  }
  for(i in 1:B_tot){
    for(j in 1:M){
      beta[i,j] ~ normal_mixture(prob_bernoulli, mean_c[i,j], 0, c, epsilon);
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
      for(j in 1:M){
        target += binomial_lpmf(k[i,j]|n[i],inv_logit(alpha[i,j]+bx[i,j]));
      }
    }
  }
}


generated quantities {
  matrix[B_tot,M] beta_bt; //beta back-transformed
  
  vector[M] mu_bt;
  vector[M] p_init;
  vector[M] p_x;
  vector[M] p_xstat;
  vector[M] p_init_pred;
  vector[M] p_x_pred;
  vector[M] p_xstat_pred;

  real randnorm;
  
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
  
  // predictive intervals
  randnorm = normal_rng(0,1);
  p_init_pred = inv_logit(mu_bt + (tau + sigma) * randnorm);
  p_x_pred = inv_logit(mu_bt + beta_bt' *x_pred + (tau + sigma) * randnorm);
  p_xstat_pred = inv_logit(mu_bt + beta_bt' *x_stat + (tau + sigma) * randnorm);
} // The prior or posterior predictive distribution

