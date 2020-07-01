// Stan model for simple linear regression
data {
  int < lower = 1 > N; // number of studies
  int < lower = 1 > M; // number of mutations
  int < lower = 1 > B[M]; // number of covariates
  int < lower = 1 > B_tot;
  int tam_pos[6]; //position of the 6 TAMs
  
  int k[N,M]; // Number of events by study
  int n[N]; // Number of people at risk by study
  
  matrix[N,B[1]] x_m1;
  matrix[N,B[2]] x_m2;
  matrix[N,B[3]] x_m3;
  matrix[N,B[4]] x_m4;
  matrix[N,B[5]] x_m5;
  matrix[N,B[6]] x_m6;
  matrix[N,B[7]] x_m7;
  matrix[N,B[8]] x_m8;
  
  vector[N-3] time_tr_nna; //ART time in non-missing position
  real mean_t; //mean of ART time
  real var_t; //variance of ART time
  int na_time_pos[3]; //position of missing ART times
  int nna_time_pos[N-3]; //position of non-missing ART times
  
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
  
  
  row_vector[B[1]] x_FTC_TDF1;
  row_vector[B[2]] x_FTC_TDF2;
  row_vector[B[3]] x_FTC_TDF3;
  row_vector[B[4]] x_FTC_TDF4;
  row_vector[B[5]] x_FTC_TDF5;
  row_vector[B[6]] x_FTC_TDF6;
  row_vector[B[7]] x_FTC_TDF7;
  row_vector[B[8]] x_FTC_TDF8;
  
  
  row_vector[B[1]] x_FTC_ZDV1;
  row_vector[B[2]] x_FTC_ZDV2;
  row_vector[B[3]] x_FTC_ZDV3;
  row_vector[B[4]] x_FTC_ZDV4;
  row_vector[B[5]] x_FTC_ZDV5;
  row_vector[B[6]] x_FTC_ZDV6;
  row_vector[B[7]] x_FTC_ZDV7;
  row_vector[B[8]] x_FTC_ZDV8;
  
  
  real init_res_1[M];
  real init_res_2[M];
  //real init_add[M];
  real < lower = 0 > sigma_1[M];
  real < lower = 0 > tau_1;
  
  vector[B[1]] beta_m1_1;
  vector[B[2]] beta_m2_1;
  vector[B[3]] beta_m3_1;
  vector[B[4]] beta_m4_1;
  vector[B[5]] beta_m5_1;
  vector[B[6]] beta_m6_1;
  vector[B[7]] beta_m7_1;
  vector[B[8]] beta_m8_1;
  
  vector < lower = 0 > [B[1]] beta_m1_2;
  vector < lower = 0 > [B[2]] beta_m2_2;
  vector < lower = 0 > [B[3]] beta_m3_2;
  vector < lower = 0 > [B[4]] beta_m4_2;
  vector < lower = 0 > [B[5]] beta_m5_2;
  vector < lower = 0 > [B[6]] beta_m6_2;
  vector < lower = 0 > [B[7]] beta_m7_2;
  vector < lower = 0 > [B[8]] beta_m8_2;
  
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
  real < lower = 0 > sigma[M];
  real < lower = 0 > tau;
  vector < lower = 0 > [3] time_imput;
}

transformed parameters {
  real alpha[N,M];
  real < lower = 0, upper = 1 > init_res[M];
  
  matrix[N,B[1]] x_imput1;
  matrix[N,B[2]] x_imput2;
  matrix[N,B[3]] x_imput3;
  matrix[N,B[4]] x_imput4;
  matrix[N,B[5]] x_imput5;
  matrix[N,B[6]] x_imput6;
  matrix[N,B[7]] x_imput7;
  matrix[N,B[8]] x_imput8;
  
  vector[3] time_imput_tr;
  time_imput_tr = (time_imput - mean_t) ./ sqrt(var_t);
  
  
  init_res=inv_logit(mu);
  
  for(i in 1:N){
    for(j in 1:M){
      alpha[i,j] = mu[j] + sigma[j] * v[i,j] + tau * u[i];
    }
  }
  
  
  x_imput1[nna_time_pos,1] = time_tr_nna;
  x_imput1[na_time_pos,1] = time_imput_tr;
  if(B[1]>1){
    for(i in 2:B[1]){
      x_imput1[,i] = x_m1[,i];
    }
  }
  x_imput2[nna_time_pos,1] = time_tr_nna;
  x_imput2[na_time_pos,1] = time_imput_tr;
  if(B[2]>1){
    for(i in 2:B[2]){
      x_imput2[,i] = x_m2[,i];
    }
  }
  x_imput3[nna_time_pos,1] = time_tr_nna;
  x_imput3[na_time_pos,1] = time_imput_tr;
  if(B[3]>1){
    for(i in 2:B[3]){
      x_imput3[,i] = x_m3[,i];
    }
  }
  x_imput4[nna_time_pos,1] = time_tr_nna;
  x_imput4[na_time_pos,1] = time_imput_tr;
  if(B[4]>1){
    for(i in 2:B[4]){
      x_imput4[,i] = x_m4[,i];
    }
  }
  x_imput5[nna_time_pos,1] = time_tr_nna;
  x_imput5[na_time_pos,1] = time_imput_tr;
  if(B[5]>1){
    for(i in 2:B[5]){
      x_imput5[,i] = x_m5[,i];
    }
  }
  x_imput6[nna_time_pos,1] = time_tr_nna;
  x_imput6[na_time_pos,1] = time_imput_tr;
  if(B[6]>1){
    for(i in 2:B[6]){
      x_imput6[,i] = x_m6[,i];
    }
  }
  x_imput7[nna_time_pos,1] = time_tr_nna;
  x_imput7[na_time_pos,1] = time_imput_tr;
  if(B[7]>1){
    for(i in 2:B[7]){
      x_imput7[,i] = x_m7[,i];
    }
  }
  x_imput8[nna_time_pos,1] = time_tr_nna;
  x_imput8[na_time_pos,1] = time_imput_tr;
  if(B[8]>1){
    for(i in 2:B[8]){
      x_imput8[,i] = x_m8[,i];
    }
  }
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
  beta_m1 ~ normal(beta_m1_1,beta_m1_2);
  beta_m2 ~ normal(beta_m2_1,beta_m2_2);
  beta_m3 ~ normal(beta_m3_1,beta_m3_2);
  beta_m4 ~ normal(beta_m4_1,beta_m4_2);
  beta_m5 ~ normal(beta_m5_1,beta_m5_2);
  beta_m6 ~ normal(beta_m6_1,beta_m6_2);
  beta_m7 ~ normal(beta_m7_1,beta_m7_2);
  beta_m8 ~ normal(beta_m8_1,beta_m8_2);
  
  
  //Priors time imputation
  time_imput[1] ~ gamma(37.7^2/var_t,mean_t/var_t);
  time_imput[2:3] ~ gamma(mean_t^2/var_t,mean_t/var_t);
  
  //Random intercept
  tau ~ exponential(tau_1);
  for(i in 1:N){
    u[i] ~ normal(0,1);
  }
  
  //Likelihood
  if(inference==1){
    for(i in 1:N){
      target += binomial_lpmf(k[i,1]|n[i],inv_logit(alpha[i,1]+x_imput1[i,]*beta_m1));
      target += binomial_lpmf(k[i,2]|n[i],inv_logit(alpha[i,2]+x_imput2[i,]*beta_m2));
      target += binomial_lpmf(k[i,3]|n[i],inv_logit(alpha[i,3]+x_imput3[i,]*beta_m3));
      target += binomial_lpmf(k[i,4]|n[i],inv_logit(alpha[i,4]+x_imput4[i,]*beta_m4));
      target += binomial_lpmf(k[i,5]|n[i],inv_logit(alpha[i,5]+x_imput5[i,]*beta_m5));
      target += binomial_lpmf(k[i,6]|n[i],inv_logit(alpha[i,6]+x_imput6[i,]*beta_m6));
      target += binomial_lpmf(k[i,7]|n[i],inv_logit(alpha[i,7]+x_imput7[i,]*beta_m7));
      target += binomial_lpmf(k[i,8]|n[i],inv_logit(alpha[i,8]+x_imput8[i,]*beta_m8));
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
  vector[B[8]] beta_bt_m8;
  
  
  vector[M] mu_bt;
  vector[M] p_init;
  matrix[N,M] p_x_study;
  vector[M] p_x_FTC_TDF;
  vector[M] p_x_FTC_ZDV;
  
  vector[N] p_noTAM_min;
  vector[N] p_noTAM_max;
  real p_noTAM_TDF_min;
  real p_noTAM_TDF_max;
  real p_noTAM_ZDV_min;
  real p_noTAM_ZDV_max;
  real p_noTAM_init_min;
  real p_noTAM_init_max;
  
  
  
  // back transformed beta and mu
  beta_bt_m1 = beta_m1 ./ x1_sd;
  beta_bt_m2 = beta_m2 ./ x2_sd;
  beta_bt_m3 = beta_m3 ./ x3_sd;
  beta_bt_m4 = beta_m4 ./ x4_sd;
  beta_bt_m5 = beta_m5 ./ x5_sd;
  beta_bt_m6 = beta_m6 ./ x6_sd;
  beta_bt_m7 = beta_m7 ./ x7_sd;
  beta_bt_m8 = beta_m8 ./ x8_sd;
  
  
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
  p_noTAM_init_min = (1-p_init[1]) * (1-p_init[3]) * (1-p_init[4]) * (1-p_init[6]) * (1-p_init[7]) * (1-p_init[8]);
  p_noTAM_init_max = 1 -max(p_init[tam_pos]);
  
  for(i in 1:N){
    p_x_study[i,1] = inv_logit(mu[1] + tau * u[i] + x_imput1[i,] * beta_m1);
    p_x_study[i,2] = inv_logit(mu[2] + tau * u[i] + x_imput2[i,] * beta_m2);
    p_x_study[i,3] = inv_logit(mu[3] + tau * u[i] + x_imput3[i,] * beta_m3);
    p_x_study[i,4] = inv_logit(mu[4] + tau * u[i] + x_imput4[i,] * beta_m4);
    p_x_study[i,5] = inv_logit(mu[5] + tau * u[i] + x_imput5[i,] * beta_m5);
    p_x_study[i,6] = inv_logit(mu[6] + tau * u[i] + x_imput6[i,] * beta_m6);
    p_x_study[i,7] = inv_logit(mu[7] + tau * u[i] + x_imput7[i,] * beta_m7);
    p_x_study[i,8] = inv_logit(mu[8] + tau * u[i] + x_imput8[i,] * beta_m8);
    
    p_noTAM_min[i] = (1-p_x_study[i,1]) * (1-p_x_study[i,3]) * (1-p_x_study[i,4]) * (1-p_x_study[i,6]) * (1-p_x_study[i,7]) * (1-p_x_study[i,8]);
    p_noTAM_max[i] = 1 - max(p_x_study[i,tam_pos]);
  }
  
  
  p_x_FTC_TDF[1] = inv_logit(mu_bt[1] + x_FTC_TDF1 * beta_bt_m1);
  p_x_FTC_TDF[2] = inv_logit(mu_bt[2] + x_FTC_TDF2 * beta_bt_m2);
  p_x_FTC_TDF[3] = inv_logit(mu_bt[3] + x_FTC_TDF3 * beta_bt_m3);
  p_x_FTC_TDF[4] = inv_logit(mu_bt[4] + x_FTC_TDF4 * beta_bt_m4);
  p_x_FTC_TDF[5] = inv_logit(mu_bt[5] + x_FTC_TDF5 * beta_bt_m5);
  p_x_FTC_TDF[6] = inv_logit(mu_bt[6] + x_FTC_TDF6 * beta_bt_m6);
  p_x_FTC_TDF[7] = inv_logit(mu_bt[7] + x_FTC_TDF7 * beta_bt_m7);
  p_x_FTC_TDF[8] = inv_logit(mu_bt[8] + x_FTC_TDF8 * beta_bt_m8);
  
  p_noTAM_TDF_min = (1-p_x_FTC_TDF[1]) * (1-p_x_FTC_TDF[3]) * (1-p_x_FTC_TDF[4]) * (1-p_x_FTC_TDF[6]) * (1-p_x_FTC_TDF[7]) * (1-p_x_FTC_TDF[8]);
  p_noTAM_TDF_max = 1 -max(p_x_FTC_TDF[tam_pos]);
  
  p_x_FTC_ZDV[1] = inv_logit(mu_bt[1] + x_FTC_ZDV1 * beta_bt_m1);
  p_x_FTC_ZDV[2] = inv_logit(mu_bt[2] + x_FTC_ZDV2 * beta_bt_m2);
  p_x_FTC_ZDV[3] = inv_logit(mu_bt[3] + x_FTC_ZDV3 * beta_bt_m3);
  p_x_FTC_ZDV[4] = inv_logit(mu_bt[4] + x_FTC_ZDV4 * beta_bt_m4);
  p_x_FTC_ZDV[5] = inv_logit(mu_bt[5] + x_FTC_ZDV5 * beta_bt_m5);
  p_x_FTC_ZDV[6] = inv_logit(mu_bt[6] + x_FTC_ZDV6 * beta_bt_m6);
  p_x_FTC_ZDV[7] = inv_logit(mu_bt[7] + x_FTC_ZDV7 * beta_bt_m7);
  p_x_FTC_ZDV[8] = inv_logit(mu_bt[8] + x_FTC_ZDV8 * beta_bt_m8);
  
  p_noTAM_ZDV_min = (1-p_x_FTC_ZDV[1]) * (1-p_x_FTC_ZDV[3]) * (1-p_x_FTC_ZDV[4]) * (1-p_x_FTC_ZDV[6]) * (1-p_x_FTC_ZDV[7]) * (1-p_x_FTC_ZDV[8]);
  p_noTAM_ZDV_max = 1 -max(p_x_FTC_ZDV[tam_pos]);
} // The prior or posterior predictive distribution
