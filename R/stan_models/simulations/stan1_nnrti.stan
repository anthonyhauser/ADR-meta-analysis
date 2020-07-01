// Stan model for simple linear regression
data {
  int < lower = 1 > N; // number of studies
  int < lower = 1 > M; // number of mutations
  int < lower = 1 > B[M]; // number of covariates
  int < lower = 1 > B_tot;
  
  int k[N,M]; // Number of events by study
  int n[N]; // Number of people at risk by study
  
  matrix[N,B[1]] x_m1; //matrix of dim N times B, B co-variates for each of the N (sub-)study
  matrix[N,B[2]] x_m2;
  matrix[N,B[3]] x_m3;
  matrix[N,B[4]] x_m4;
  matrix[N,B[5]] x_m5;
  matrix[N,B[6]] x_m6;
  matrix[N,B[7]] x_m7;

  
  row_vector[B[1]] x1_mean;
  row_vector[B[2]] x2_mean;
  row_vector[B[3]] x3_mean;
  row_vector[B[4]] x4_mean;
  row_vector[B[5]] x5_mean;
  row_vector[B[6]] x6_mean;
  row_vector[B[7]] x7_mean;
  
  vector[B[1]] x1_sd;
  vector[B[2]] x2_sd;
  vector[B[3]] x3_sd;
  vector[B[4]] x4_sd;
  vector[B[5]] x5_sd;
  vector[B[6]] x6_sd;
  vector[B[7]] x7_sd;
  
  row_vector[B[1]] x_pred1;
  row_vector[B[2]] x_pred2;
  row_vector[B[3]] x_pred3;
  row_vector[B[4]] x_pred4;
  row_vector[B[5]] x_pred5;
  row_vector[B[6]] x_pred6;
  row_vector[B[7]] x_pred7;
  
  row_vector[B[1]] x_stat1;
  row_vector[B[2]] x_stat2;
  row_vector[B[3]] x_stat3;
  row_vector[B[4]] x_stat4;
  row_vector[B[5]] x_stat5;
  row_vector[B[6]] x_stat6;
  row_vector[B[7]] x_stat7;
  
  real < lower = 0 > init_res_1[M];
  real < lower = 0 > init_res_2[M];
  real init_add[M];
  real < lower = 0 > sigma_1[M];
  real < lower = 0 > tau_1;
  
  vector[B[1]] beta_m1_1;
  vector[B[2]] beta_m2_1;
  vector[B[3]] beta_m3_1;
  vector[B[4]] beta_m4_1;
  vector[B[5]] beta_m5_1;
  vector[B[6]] beta_m6_1;
  vector[B[7]] beta_m7_1;
  
  vector < lower = 0 > [B[1]] beta_m1_2;
  vector < lower = 0 > [B[2]] beta_m2_2;
  vector < lower = 0 > [B[3]] beta_m3_2;
  vector < lower = 0 > [B[4]] beta_m4_2;
  vector < lower = 0 > [B[5]] beta_m5_2;
  vector < lower = 0 > [B[6]] beta_m6_2;
  vector < lower = 0 > [B[7]] beta_m7_2;
  
  int inference;
}

parameters {
  vector[B[1]] beta_m1;
  vector[B[2]] beta_m2;
  vector[B[3]] beta_m3;
  vector[B[4]] beta_m4;
  vector[B[5]] beta_m5;
  vector[B[6]] beta_m6;
  vector[B[7]] beta_m7;
  
  real < lower = 0, upper = 1 > init_res[M];
  
  real u[N];
  real v[N,M];
  real < lower = 0 > sigma[M];
  real < lower = 0 > tau;
  
}

transformed parameters {
  real mu[M];
  real alpha[N,M];
  for(i in 1:M){
    mu[i] = logit(init_res[i]) + init_add[i];
  }
  //mu = logit(init_res) + init_add;
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
  beta_m1 ~ normal(beta_m1_1,beta_m1_2);
  beta_m2 ~ normal(beta_m2_1,beta_m2_2);
  beta_m3 ~ normal(beta_m3_1,beta_m3_2);
  beta_m4 ~ normal(beta_m4_1,beta_m4_2);
  beta_m5 ~ normal(beta_m5_1,beta_m5_2);
  beta_m6 ~ normal(beta_m6_1,beta_m6_2);
  beta_m7 ~ normal(beta_m7_1,beta_m7_2);
  
  
  //Random intercept
  tau ~ exponential(tau_1);
  for(i in 1:N){
    u[i] ~ normal(0,1);
  }
  
  //Likelihood
  if(inference==1){
    for(i in 1:N){
      target += binomial_lpmf(k[i,1]|n[i],inv_logit(alpha[i,1]+x_m1[i,]*beta_m1));
      target += binomial_lpmf(k[i,2]|n[i],inv_logit(alpha[i,2]+x_m2[i,]*beta_m2));
      target += binomial_lpmf(k[i,3]|n[i],inv_logit(alpha[i,3]+x_m3[i,]*beta_m3));
      target += binomial_lpmf(k[i,4]|n[i],inv_logit(alpha[i,4]+x_m4[i,]*beta_m4));
      target += binomial_lpmf(k[i,5]|n[i],inv_logit(alpha[i,5]+x_m5[i,]*beta_m5));
      target += binomial_lpmf(k[i,6]|n[i],inv_logit(alpha[i,6]+x_m6[i,]*beta_m6));
      target += binomial_lpmf(k[i,7]|n[i],inv_logit(alpha[i,7]+x_m7[i,]*beta_m7));
    }
  }
}


generated quantities {
  // matrix[B_tot,M] beta_bt; //matrix
  
  vector[B[1]] beta_bt_m1; //beta back-transformed
  vector[B[2]] beta_bt_m2;
  vector[B[3]] beta_bt_m3;
  vector[B[4]] beta_bt_m4;
  vector[B[5]] beta_bt_m5;
  vector[B[6]] beta_bt_m6;
  vector[B[7]] beta_bt_m7;
  
  vector[M] mu_bt;
  vector[M] p_init;
  vector[M] p_x;
  vector[M] p_xstat;
  vector[M] p_init_pred;
  vector[M] p_x_pred;
  vector[M] p_xstat_pred;
  
  real randnorm;
  
  // back transformed beta and mu
  beta_bt_m1 = beta_m1 ./ x1_sd;
  beta_bt_m2 = beta_m2 ./ x2_sd;
  beta_bt_m3 = beta_m3 ./ x3_sd;
  beta_bt_m4 = beta_m4 ./ x4_sd;
  beta_bt_m5 = beta_m5 ./ x5_sd;
  beta_bt_m6 = beta_m6 ./ x6_sd;
  beta_bt_m7 = beta_m7 ./ x7_sd;
  
  mu_bt[1] = mu[1] - x1_mean * beta_bt_m1;
  mu_bt[2] = mu[2] - x2_mean * beta_bt_m2;
  mu_bt[3] = mu[3] - x3_mean * beta_bt_m3;
  mu_bt[4] = mu[4] - x4_mean * beta_bt_m4;
  mu_bt[5] = mu[5] - x5_mean * beta_bt_m5;
  mu_bt[6] = mu[6] - x6_mean * beta_bt_m6;
  mu_bt[7] = mu[7] - x7_mean * beta_bt_m7;

  // credibility intervals
  p_init = inv_logit(mu_bt);
  
  p_x[1] = inv_logit(mu_bt[1] + x_pred1 * beta_bt_m1);
  p_x[2] = inv_logit(mu_bt[2] + x_pred2 * beta_bt_m2);
  p_x[3] = inv_logit(mu_bt[3] + x_pred3 * beta_bt_m3);
  p_x[4] = inv_logit(mu_bt[4] + x_pred4 * beta_bt_m4);
  p_x[5] = inv_logit(mu_bt[5] + x_pred5 * beta_bt_m5);
  p_x[6] = inv_logit(mu_bt[6] + x_pred6 * beta_bt_m6);
  p_x[7] = inv_logit(mu_bt[7] + x_pred7 * beta_bt_m7);
  
  p_xstat[1] = inv_logit(mu_bt[1] + x_stat1 * beta_bt_m1);
  p_xstat[2] = inv_logit(mu_bt[2] + x_stat2 * beta_bt_m2);
  p_xstat[3] = inv_logit(mu_bt[3] + x_stat3 * beta_bt_m3);
  p_xstat[4] = inv_logit(mu_bt[4] + x_stat4 * beta_bt_m4);
  p_xstat[5] = inv_logit(mu_bt[5] + x_stat5 * beta_bt_m5);
  p_xstat[6] = inv_logit(mu_bt[6] + x_stat6 * beta_bt_m6);
  p_xstat[7] = inv_logit(mu_bt[7] + x_stat7 * beta_bt_m7);
  
  // predictive intervals
  randnorm = normal_rng(0,1);
  p_init_pred[1] = inv_logit(mu_bt[1] + (tau + sigma[1]) * randnorm);
  p_init_pred[2] = inv_logit(mu_bt[2] + (tau + sigma[2]) * randnorm);
  p_init_pred[3] = inv_logit(mu_bt[3] + (tau + sigma[3]) * randnorm);
  p_init_pred[4] = inv_logit(mu_bt[4] + (tau + sigma[4]) * randnorm);
  p_init_pred[5] = inv_logit(mu_bt[5] + (tau + sigma[5]) * randnorm);
  p_init_pred[6] = inv_logit(mu_bt[6] + (tau + sigma[6]) * randnorm);
  p_init_pred[7] = inv_logit(mu_bt[7] + (tau + sigma[7]) * randnorm);
  
  p_x_pred[1] = inv_logit(mu_bt[1] + x_pred1 * beta_bt_m1 + (tau + sigma[1]) * randnorm);
  p_x_pred[2] = inv_logit(mu_bt[2] + x_pred2 * beta_bt_m2 + (tau + sigma[2]) * randnorm);
  p_x_pred[3] = inv_logit(mu_bt[3] + x_pred3 * beta_bt_m3 + (tau + sigma[3]) * randnorm);
  p_x_pred[4] = inv_logit(mu_bt[4] + x_pred4 * beta_bt_m4 + (tau + sigma[4]) * randnorm);
  p_x_pred[5] = inv_logit(mu_bt[5] + x_pred5 * beta_bt_m5 + (tau + sigma[5]) * randnorm);
  p_x_pred[6] = inv_logit(mu_bt[6] + x_pred6 * beta_bt_m6 + (tau + sigma[6]) * randnorm);
  p_x_pred[7] = inv_logit(mu_bt[7] + x_pred7 * beta_bt_m7 + (tau + sigma[7]) * randnorm);
  
  p_xstat_pred[1] = inv_logit(mu_bt[1] + x_stat1 * beta_bt_m1 + (tau + sigma[1]) * randnorm);
  p_xstat_pred[2] = inv_logit(mu_bt[2] + x_stat2 * beta_bt_m2 + (tau + sigma[2]) * randnorm);
  p_xstat_pred[3] = inv_logit(mu_bt[3] + x_stat3 * beta_bt_m3 + (tau + sigma[3]) * randnorm);
  p_xstat_pred[4] = inv_logit(mu_bt[4] + x_stat4 * beta_bt_m4 + (tau + sigma[4]) * randnorm);
  p_xstat_pred[5] = inv_logit(mu_bt[5] + x_stat5 * beta_bt_m5 + (tau + sigma[5]) * randnorm);
  p_xstat_pred[6] = inv_logit(mu_bt[6] + x_stat6 * beta_bt_m6 + (tau + sigma[6]) * randnorm);
  p_xstat_pred[7] = inv_logit(mu_bt[7] + x_stat7 * beta_bt_m7 + (tau + sigma[7]) * randnorm);
  
} // The prior or posterior predictive distribution

