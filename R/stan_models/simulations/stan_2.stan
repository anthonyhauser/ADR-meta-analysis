// Stan model for simple linear regression
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
  //real init_add[M];
  real < lower = 0 > sigma_1[M];
  real < lower = 0 > tau_1;
  
  matrix[B_tot,M] beta_1;
  matrix[B_tot,M] beta_2;
  
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
      beta[i,j] ~ normal(beta_1[i,j],beta_2[i,j]);
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

  // back transformed beta and mu
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

  // // predictive intervals
  randnorm = normal_rng(0,1);
  p_init_pred = inv_logit(mu_bt + (tau + sigma) * randnorm);
  p_x_pred = inv_logit(mu_bt + beta_bt' *x_pred + (tau + sigma) * randnorm);
  p_xstat_pred = inv_logit(mu_bt + beta_bt' *x_stat + (tau + sigma) * randnorm);
  
  // real random_v;
  // 
    // real p0_1;
  // 
    // real p1;
  // real p2;
  // real p3;
  // real p4;
  // real p5;
  // real p6;
  // real p7;
  // real p8;
  // real p1_pred;
  // real p2_pred;
  // real p3_pred;
  // real p4_pred;
  // real p5_pred;
  // real p6_pred;
  // real p7_pred;
  // real p8_pred;
  // 
    // 
    // random_v = normal_rng(0,1);
  // 
    // p0_1 = inv_logit(mu[1]-(x1_mean ./ x1_sd)*beta_m1);
  // p0_1 = inv_logit(mu[1]-(x1_mean ./ x1_sd)*beta_m1 + (x_pred1 ./ x1_sd) * beta_m1);
  // //p1 = inv_logit(mu[1]-(x_mean[B[1]] ./ x_sd[B[1]])*beta_m1);
  // 
    // p1 = inv_logit(mu[1]+x_pred1*beta_m1);
  // p2 = inv_logit(mu[2]+x_pred2*beta_m2);
  // p3 = inv_logit(mu[3]+x_pred3*beta_m3);
  // p4 = inv_logit(mu[4]+x_pred4*beta_m4);
  // p5 = inv_logit(mu[5]+x_pred5*beta_m5);
  // p6 = inv_logit(mu[6]+x_pred6*beta_m6);
  // p7 = inv_logit(mu[7]+x_pred7*beta_m7);
  // p8 = inv_logit(mu[8]+x_pred8*beta_m8);
  // 
    // 
    // 
    // p1_pred = inv_logit(mu[1]+x_pred1*beta_m1+(tau+sigma[1])*random_v);
  // p2_pred = inv_logit(mu[2]+x_pred2*beta_m2+(tau+sigma[2])*random_v);
  // p3_pred = inv_logit(mu[3]+x_pred3*beta_m3+(tau+sigma[3])*random_v);
  // p4_pred = inv_logit(mu[4]+x_pred4*beta_m4+(tau+sigma[4])*random_v);
  // p5_pred = inv_logit(mu[5]+x_pred5*beta_m5+(tau+sigma[5])*random_v);
  // p6_pred = inv_logit(mu[6]+x_pred6*beta_m6+(tau+sigma[6])*random_v);
  // p7_pred = inv_logit(mu[7]+x_pred7*beta_m7+(tau+sigma[7])*random_v);
  // p8_pred = inv_logit(mu[8]+x_pred8*beta_m8+(tau+sigma[8])*random_v);
} // The prior or posterior predictive distribution

