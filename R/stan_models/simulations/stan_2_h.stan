// Stan model for simple linear regression
data {
  int < lower = 1 > N; // number of studies
  int < lower = 1 > M; // number of mutations
  int < lower = 1 > B[M]; // number of covariates
  int < lower = 1 > B_tot;
  
  int k[N,M]; // Number of events by study
  int n[N]; // Number of people at risk by study
  
  matrix[N,B_tot] x; //matrix of dim N times B, B co-variates for each of the N (sub-)study
  
  //for generated quantities
  row_vector[B_tot] x_mean;
  vector[B_tot] x_sd;
  vector[B_tot] x_pred;
  vector[B_tot] x_stat;
  
  //model hyperparameters for priors
  real init_res_1[M];
  real init_res_2[M];
  real init_add[M];
  real < lower = 0 > sigma_1[M];
  real < lower = 0 > tau_1;
  
  //horseshoe hyperparameters
  vector < lower =0 > [M] scale_global; // scale for the half -t prior for tau
  vector < lower =1 > [M] nu_global; // degrees of freedom for the half -t priors for tau
  vector < lower =1 > [M] nu_local; // degrees of freedom for the half - t priors for lambdas
  vector < lower =0 > [M] slab_scale; // slab scale for the regularized horseshoe
  vector < lower =0 > [M] slab_df; // slab degrees of freedom for the regularized horseshoe
  
  int inference;
}

parameters {
  //intercept (before logit)
  //real < lower = 0, upper = 1 > init_res[M];
  real mu[M];
  
  //heterogeneity
  real u[N];
  real v[N,M];
  vector < lower = 0 > [M] sigma;
  real < lower = 0 > tau;
  
  //horseshoe
  matrix [B_tot,M] z;
  vector < lower =0 > [M] tau_h; // global shrinkage parameter
  matrix < lower =0 > [B_tot,M] lambda; // local shrinkage parameter
  vector < lower =0 > [M] caux;
  
}

transformed parameters {
  //Intercept, intercept with heterogeneity
  //real mu[M];
  real alpha[N,M];
  matrix [B_tot,M] beta; 
  matrix[N,M] xb;
  vector < lower =0 > [M] c; // slab scale
  matrix < lower =0 >[B_tot,M] lambda_tilde ; // ’ truncated ’ local shrinkage parameter
  
  //mu = logit(init_res);
  // for(i in 1:M){
  //   mu[i] = logit(init_res[i]) + init_add[i];
  // }
  
  for(i in 1:N){
    for(j in 1:M){
      alpha[i,j] = mu[j] + sigma[j] * v[i,j] + tau * u[i];
    }
  }
  
  //beta horseshoe
  c = slab_scale .* sqrt (caux);
  for(i in 1:M){
    lambda_tilde[,i] = sqrt ( c[i] ^2 * square ( lambda[,i] ) ./ (c[i] ^2 + tau_h[i] ^2* square ( lambda[,i] )) );
    beta[,i] = z[,i] .* lambda_tilde[,i] * tau_h[i];
  }
  xb = x * beta;
}

model {
  //Priors intercept and heterogeneity
  //init_res ~ beta(init_res_1,init_res_2);
  mu ~ normal(init_res_1,init_res_2);
  for(j in 1:M){
    sigma[j] ~ exponential(sigma_1[j]);
    for(i in 1:N){
      v[i,j] ~ normal(0,1);
    }
  }
  tau ~ exponential(tau_1);
  for(i in 1:N){
    u[i] ~ normal(0,1);
  }
  
  
  //Priors beta
  for(i in 1:B_tot){
    for(j in 1:M){
      z[i,j] ~ normal(0,1);
      lambda[i,j] ~ student_t( nu_local[j], 0, 1);
    }
  }
  for(j in 1:M){
    tau_h[j] ~ student_t( nu_global[j], 0 , scale_global[j]);
    caux[j] ~ inv_gamma (0.5* slab_df[j] , 0.5* slab_df[j] );
  }
  
  
  //Likelihood
  if(inference==1){
    for(i in 1:N){
      for(j in 1:M){
        target += binomial_lpmf(k[i,j]|n[i],inv_logit(alpha[i,j]+xb[i,j]));
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

