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
  
  matrix[N,B[1]] x_m1; //matrix of dim N times B, B co-variates for each of the N (sub-)study
  matrix[N,B[2]] x_m2;
  matrix[N,B[3]] x_m3;
  matrix[N,B[4]] x_m4;
  matrix[N,B[5]] x_m5;
  matrix[N,B[6]] x_m6;
  matrix[N,B[7]] x_m7;
  matrix[N,B[8]] x_m8;
  
  row_vector[B[1]] x1_mean;
  row_vector[B[2]] x2_mean;
  row_vector[B[3]] x3_mean;
  row_vector[B[4]] x4_mean;
  row_vector[B[5]] x5_mean;
  row_vector[B[6]] x6_mean;
  row_vector[B[7]] x7_mean;
  row_vector[B[8]] x8_mean;
  
  vector[B[1]] x1_sd;
  vector[B[2]] x2_sd;
  vector[B[3]] x3_sd;
  vector[B[4]] x4_sd;
  vector[B[5]] x5_sd;
  vector[B[6]] x6_sd;
  vector[B[7]] x7_sd;
  vector[B[8]] x8_sd;
  
  row_vector[B[1]] x_pred1;
  row_vector[B[2]] x_pred2;
  row_vector[B[3]] x_pred3;
  row_vector[B[4]] x_pred4;
  row_vector[B[5]] x_pred5;
  row_vector[B[6]] x_pred6;
  row_vector[B[7]] x_pred7;
  row_vector[B[8]] x_pred8;
  
  row_vector[B[1]] x_stat1;
  row_vector[B[2]] x_stat2;
  row_vector[B[3]] x_stat3;
  row_vector[B[4]] x_stat4;
  row_vector[B[5]] x_stat5;
  row_vector[B[6]] x_stat6;
  row_vector[B[7]] x_stat7;
  row_vector[B[8]] x_stat8;
  
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
  vector[B[1]] beta_m1;
  vector[B[2]] beta_m2;
  vector[B[3]] beta_m3;
  vector[B[4]] beta_m4;
  vector[B[5]] beta_m5;
  vector[B[6]] beta_m6;
  vector[B[7]] beta_m7;
  vector[B[8]] beta_m8;
  
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
  mu ~ normal(init_res_1,init_res_2);
  //init_res ~ beta(init_res_1,init_res_2);
  for(j in 1:M){
    sigma[j] ~ exponential(sigma_1[j]);
    for(i in 1:N){
      v[i,j] ~ normal(0,1);
    }
  }

  for(i in 1:B[1]){
    beta_m1[i] ~ normal_mixture(prob_bernoulli, mean_c[i,1], 0, c, epsilon);
  }
  for(i in 1:B[2]){
    beta_m2[i] ~ normal_mixture(prob_bernoulli, mean_c[i,2], 0, c, epsilon);
  }
  for(i in 1:B[3]){
    beta_m3[i] ~ normal_mixture(prob_bernoulli, mean_c[i,3], 0, c, epsilon);
  }
  for(i in 1:B[4]){
    beta_m4[i] ~ normal_mixture(prob_bernoulli, mean_c[i,4], 0, c, epsilon);
  }
  for(i in 1:B[5]){
    beta_m5[i] ~ normal_mixture(prob_bernoulli, mean_c[i,5], 0, c, epsilon);
  }
  for(i in 1:B[6]){
    beta_m6[i] ~ normal_mixture(prob_bernoulli, mean_c[i,6], 0, c, epsilon);
  }
  for(i in 1:B[7]){
    beta_m7[i] ~ normal_mixture(prob_bernoulli, mean_c[i,7], 0, c, epsilon);
  }
  for(i in 1:B[8]){
    beta_m8[i] ~ normal_mixture(prob_bernoulli, mean_c[i,8], 0, c, epsilon);
  }
  

  
  
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
      target += binomial_lpmf(k[i,8]|n[i],inv_logit(alpha[i,8]+x_m8[i,]*beta_m8));
    }
  }
}


generated quantities {
  vector[B[1]] beta_bt_m1; //beta back-transformed
  vector[B[2]] beta_bt_m2;
  vector[B[3]] beta_bt_m3;
  vector[B[4]] beta_bt_m4;
  vector[B[5]] beta_bt_m5;
  vector[B[6]] beta_bt_m6;
  vector[B[7]] beta_bt_m7;
  vector[B[8]] beta_bt_m8;
  
  
  vector[M] mu_bt;
  vector[M] p_init;
  vector[M] p_x;
  vector[M] p_xstat;

  
  // back transformed beta and mu
  beta_bt_m1 = beta_m1 ./ x1_sd;
  beta_bt_m2 = beta_m2 ./ x2_sd;
  beta_bt_m3 = beta_m3 ./ x3_sd;
  beta_bt_m4 = beta_m4 ./ x4_sd;
  beta_bt_m5 = beta_m5 ./ x5_sd;
  beta_bt_m6 = beta_m6 ./ x6_sd;
  beta_bt_m7 = beta_m7 ./ x7_sd;
  beta_bt_m8 = beta_m8 ./ x8_sd;
  
  // for(i in 1:M){
  //   for(j in 1:B_tot){
  //      beta_bt[j,i] = 0;
  //   }
  // }
  //  beta_bt[B1,1] = beta_bt_m1;
  //  beta_bt[B2,2] = beta_bt_m2;
  //  beta_bt[B3,3] = beta_bt_m3;
  //  beta_bt[B4,4] = beta_bt_m4;
  //  beta_bt[B5,5] = beta_bt_m5;
  //  beta_bt[B6,6] = beta_bt_m6;
  //  beta_bt[B7,7] = beta_bt_m7;
  //  beta_bt[B8,8] = beta_bt_m8;

     
  
  mu_bt[1] = mu[1] - x1_mean * beta_bt_m1;
  mu_bt[2] = mu[2] - x2_mean * beta_bt_m2;
  mu_bt[3] = mu[3] - x3_mean * beta_bt_m3;
  mu_bt[4] = mu[4] - x4_mean * beta_bt_m4;
  mu_bt[5] = mu[5] - x5_mean * beta_bt_m5;
  mu_bt[6] = mu[6] - x6_mean * beta_bt_m6;
  mu_bt[7] = mu[7] - x7_mean * beta_bt_m7;
  mu_bt[8] = mu[8] - x8_mean * beta_bt_m8;
  
  // credibility intervals
  p_init = inv_logit(mu_bt);
  
  p_x[1] = inv_logit(mu_bt[1] + x_pred1 * beta_bt_m1);
  p_x[2] = inv_logit(mu_bt[2] + x_pred2 * beta_bt_m2);
  p_x[3] = inv_logit(mu_bt[3] + x_pred3 * beta_bt_m3);
  p_x[4] = inv_logit(mu_bt[4] + x_pred4 * beta_bt_m4);
  p_x[5] = inv_logit(mu_bt[5] + x_pred5 * beta_bt_m5);
  p_x[6] = inv_logit(mu_bt[6] + x_pred6 * beta_bt_m6);
  p_x[7] = inv_logit(mu_bt[7] + x_pred7 * beta_bt_m7);
  p_x[8] = inv_logit(mu_bt[8] + x_pred8 * beta_bt_m8);
  
  p_xstat[1] = inv_logit(mu_bt[1] + x_stat1 * beta_bt_m1);
  p_xstat[2] = inv_logit(mu_bt[2] + x_stat2 * beta_bt_m2);
  p_xstat[3] = inv_logit(mu_bt[3] + x_stat3 * beta_bt_m3);
  p_xstat[4] = inv_logit(mu_bt[4] + x_stat4 * beta_bt_m4);
  p_xstat[5] = inv_logit(mu_bt[5] + x_stat5 * beta_bt_m5);
  p_xstat[6] = inv_logit(mu_bt[6] + x_stat6 * beta_bt_m6);
  p_xstat[7] = inv_logit(mu_bt[7] + x_stat7 * beta_bt_m7);
  p_xstat[8] = inv_logit(mu_bt[8] + x_stat8 * beta_bt_m8);
  
  
} // The prior or posterior predictive distribution

