// Stan model for simple linear regression
data {
  int < lower = 1 > N; // number of studies
  int < lower = 1 > M; // number of mutations
  int < lower = 1 > B[M]; // number of covariates
  
  int k[N,M]; // Number of events by study
  int n[N]; // Number of people at risk by study
  
  matrix[N,B[1]] x_m1; //matrix of dim N times B, B co-variates for each of the N (sub-)study
  matrix[N,B[2]] x_m2;
  matrix[N,B[3]] x_m3;
  matrix[N,B[4]] x_m4;
  matrix[N,B[5]] x_m5;
  matrix[N,B[6]] x_m6;
  matrix[N,B[7]] x_m7;
  
  row_vector[B[1]] x_pred1;
  row_vector[B[2]] x_pred2;
  row_vector[B[3]] x_pred3;
  row_vector[B[4]] x_pred4;
  row_vector[B[5]] x_pred5;
  row_vector[B[6]] x_pred6;
  row_vector[B[7]] x_pred7;
  
  // real x_pred3[B[3]];
  // real x_pred4[B[4]];
  // real x_pred5[B[5]];
  // real x_pred6[B[6]];
  // real x_pred7[B[7]];
  // real x_pred8[B[8]];
  
  
  real < lower = 0 > init_res_1[M];
  real < lower = 0 > init_res_2[M];
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
  //real a;
  
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
  
  //a ~ normal(0,1);
  
  // print("---------------------------------")
  // print(x_m7[1,]*beta_m7)
  // print("+++++++++++++++++++++++++++++++++++")
  // print(inv_logit(alpha[1,7]+x_m7[1,]*beta_m7))
  // print("x_pred2")
  // print(x_pred2)
  // print("beta_m2")
  // print(beta_m2)
  
  //Likelihood
  if(inference==1){
    for(i in 1:N){
      target += binomial_lpmf(k[i,1]|n[i],inv_logit(alpha[i,1]+x_m1[i,]*beta_m1/1));
      target += binomial_lpmf(k[i,2]|n[i],inv_logit(alpha[i,2]+x_m2[i,]*beta_m2/1));
      target += binomial_lpmf(k[i,3]|n[i],inv_logit(alpha[i,3]+x_m3[i,]*beta_m3/1));
      target += binomial_lpmf(k[i,4]|n[i],inv_logit(alpha[i,4]+x_m4[i,]*beta_m4/1));
      target += binomial_lpmf(k[i,5]|n[i],inv_logit(alpha[i,5]+x_m5[i,]*beta_m5/1));
      target += binomial_lpmf(k[i,6]|n[i],inv_logit(alpha[i,6]+x_m6[i,]*beta_m6/1));
      target += binomial_lpmf(k[i,7]|n[i],inv_logit(alpha[i,7]+x_m7[i,]*beta_m7/1));
    }
  }
}


generated quantities {
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
  // p1_pred = inv_logit(mu[1]+x_pred1*beta_m1+(tau+sigma[1])*v[1,1]);
  // p2_pred = inv_logit(mu[2]+x_pred2*beta_m2+(tau+sigma[2])*v[1,2]);
  // p3_pred = inv_logit(mu[3]+x_pred3*beta_m3+(tau+sigma[3])*v[1,3]);
  // p4_pred = inv_logit(mu[4]+x_pred4*beta_m4+(tau+sigma[4])*v[1,4]);
  // p5_pred = inv_logit(mu[5]+x_pred5*beta_m5+(tau+sigma[5])*v[1,5]);
  // p6_pred = inv_logit(mu[6]+x_pred6*beta_m6+(tau+sigma[6])*v[1,6]);
  // p7_pred = inv_logit(mu[7]+x_pred7*beta_m7+(tau+sigma[7])*v[1,7]);
  // p8_pred = inv_logit(mu[8]+x_pred8*beta_m8+(tau+sigma[8])*v[1,8]);
  //p_min[1] = inv_logit(mu+x_m1[1,]*beta_m1/100);
  
  // real p_pred_min = inv_logit(normal_rng(mu,sigma));
  // real p_pred_max = inv_logit(normal_rng(mu+beta[1]*50/500+beta[2]*100/500,sigma));
} // The prior or posterior predictive distribution

