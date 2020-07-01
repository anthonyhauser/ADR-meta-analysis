// Stan model for simple linear regression
functions {
  real min_abs_2real(real x, real y) {
    vector[2] v;
    v[1]=fabs(x);
    v[2]=fabs(y);
    return min(v);
  }
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
  
  real < lower = 0 > init_res_1[M];
  real < lower = 0 > init_res_2[M];
  real init_add[M];
  real < lower = 0 > sigma_1[M];
  real < lower = 0 > tau_1;
  
  real beta_p1;
  real beta_p2;
  vector[2] beta_time;
  vector[2] beta_3l;
  vector[2] beta_4l;
  
  real epsilon;
  real prob_bernoulli;
  
  int inference;
}

parameters {
  
  real beta_11;
  real beta_12;
  real beta_13;
  real beta_14;
  real beta_15;
  
  real beta_21;
  real beta_22;
  // real beta_25;
  // real <lower = beta_25> beta_23;
  // real <lower = beta_25> beta_24;
  real beta_23;
  real beta_24;
  real <lower = -min_abs_2real(beta_23,beta_24),upper = min_abs_2real(beta_23,beta_24)> beta_25;
  
  real beta_31;
  real beta_34;
  real beta_32;
  real beta_33; //<lower = beta_32>
  real beta_35; //<lower = beta_32>
      
  real beta_41;
  real beta_42;
  real beta_45;
  real beta_43; //<lower = beta_45>
  real beta_44; //<lower = beta_45>
        
  real beta_51;
  real beta_52;
  real < lower = -fabs(beta_52), upper = fabs(beta_52)> beta_53;
  real < lower = -fabs(beta_52), upper = fabs(beta_52)> beta_54;
  real < lower = -fabs(beta_52), upper = fabs(beta_52)> beta_55;
      
  real beta_61;
  real beta_64;
  real beta_62;
  real beta_63; //<lower = beta_62>
  real beta_65; //<lower = beta_62>
          
  real beta_71;
  real beta_72;
  // real beta_74;
  // real <lower = beta_72> beta_73;
  // real <lower = beta_72> beta_75;
  
  real beta_73;
  real beta_75;
  real <lower = -min_abs_2real(beta_73,beta_75),upper = min_abs_2real(beta_73,beta_75)> beta_74;
        
  // real beta_71;
  // real beta_72;
  // real <lower = beta_72> beta_74;
  // real <lower = beta_74> beta_73;
  // real <lower = beta_74> beta_75;
        
  real beta_81;
  real beta_82;
  real beta_83;
  real beta_84;
  real beta_85;
        
        
  real < lower = 0, upper = 1 > init_res[M];
        
  real u[N];
  real v[N,M];
  vector < lower = 0 > [M] sigma;
  real < lower = 0 > tau;
        
        
}

transformed parameters {
  real mu[M];
  real alpha[N,M];
  //mu = logit(init_res);
  for(i in 1:M){
    mu[i] = logit(init_res[i]) + init_add[i];
  }
  for(i in 1:N){
    for(j in 1:M){
      alpha[i,j] = mu[j] + sigma[j] * v[i,j] + tau * u[i];
    }
  }
  
}

model {
  //Priors
  beta_11 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_12 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_13 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_14 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_15 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  beta_21 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_22 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_23 ~ normal_mixture(prob_bernoulli,beta_3l[1],0,beta_3l[2], epsilon);
  beta_24 ~ normal_mixture(prob_bernoulli,beta_3l[1],0,beta_3l[2], epsilon);
  beta_25 ~ normal_mixture(prob_bernoulli,beta_3l[1],0,beta_3l[2], epsilon);
  
  beta_31 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_32 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_33 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_34 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_35 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  beta_41 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_42 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_43 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_44 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_45 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  beta_51 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_52 ~ normal_mixture(prob_bernoulli,beta_4l[1],0,beta_4l[2], epsilon);
  beta_53 ~ normal_mixture(prob_bernoulli,beta_4l[1],0,beta_4l[2], epsilon);
  beta_54 ~ normal_mixture(prob_bernoulli,beta_4l[1],0,beta_4l[2], epsilon);
  beta_55 ~ normal_mixture(prob_bernoulli,beta_4l[1],0,beta_4l[2], epsilon);
  
  beta_61 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_62 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_63 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_64 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_65 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  beta_71 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_72 ~ normal_mixture(prob_bernoulli,beta_3l[1],0,beta_3l[2], epsilon);
  beta_73 ~ normal_mixture(prob_bernoulli,beta_3l[1],0,beta_3l[2], epsilon);
  beta_74 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_75 ~ normal_mixture(prob_bernoulli,beta_3l[1],0,beta_3l[2], epsilon);
  
  beta_81 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_82 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_83 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_84 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_85 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  
  init_res ~ beta(init_res_1,init_res_2);
  for(j in 1:M){
    sigma[j] ~ exponential(sigma_1[j]);
    for(i in 1:N){
      v[i,j] ~ normal(0,1);
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
      target += binomial_lpmf(k[i,1]|n[i],inv_logit(alpha[i,1]+beta_11*x[i,1]+beta_12*x[i,2]+beta_13*x[i,3]+beta_14*x[i,4]+beta_15*x[i,5]));
      target += binomial_lpmf(k[i,2]|n[i],inv_logit(alpha[i,2]+beta_21*x[i,1]+beta_22*x[i,2]+beta_23*x[i,3]+beta_24*x[i,4]+beta_25*x[i,5]));
      target += binomial_lpmf(k[i,3]|n[i],inv_logit(alpha[i,3]+beta_31*x[i,1]+beta_32*x[i,2]+beta_33*x[i,3]+beta_34*x[i,4]+beta_35*x[i,5]));
      target += binomial_lpmf(k[i,4]|n[i],inv_logit(alpha[i,4]+beta_41*x[i,1]+beta_42*x[i,2]+beta_43*x[i,3]+beta_44*x[i,4]+beta_45*x[i,5]));
      target += binomial_lpmf(k[i,5]|n[i],inv_logit(alpha[i,5]+beta_51*x[i,1]+beta_52*x[i,2]+beta_53*x[i,3]+beta_54*x[i,4]+beta_55*x[i,5]));
      target += binomial_lpmf(k[i,6]|n[i],inv_logit(alpha[i,6]+beta_61*x[i,1]+beta_62*x[i,2]+beta_63*x[i,3]+beta_64*x[i,4]+beta_65*x[i,5]));
      target += binomial_lpmf(k[i,7]|n[i],inv_logit(alpha[i,7]+beta_71*x[i,1]+beta_72*x[i,2]+beta_73*x[i,3]+beta_74*x[i,4]+beta_75*x[i,5]));
      target += binomial_lpmf(k[i,8]|n[i],inv_logit(alpha[i,8]+beta_81*x[i,1]+beta_82*x[i,2]+beta_83*x[i,3]+beta_84*x[i,4]+beta_85*x[i,5]));
    }
  }
}


generated quantities {
  // matrix[B_tot,M] beta_bt; //beta back-transformed
  // 
    // 
    // 
    // vector[M] mu_bt;
  // vector[M] p_init;
  // vector[M] p_x;
  // vector[M] p_xstat;
  // vector[M] p_init_pred;
  // vector[M] p_x_pred;
  // vector[M] p_xstat_pred;
  // 
    // real randnorm;
  // 
    // // back transformed beta and mu
  // for(i in 1:M){
    //   beta_bt[,i] = beta[,i] ./ x_sd;
    // }
  // 
    // for(i in 1:M){
      //   mu_bt[i] = mu[i] - x_mean * beta_bt[,i];
      // }
  // 
    // // credibility intervals
  // p_init = inv_logit(mu_bt);
  // p_x = inv_logit(mu_bt + beta_bt' *x_pred);
  // p_xstat = inv_logit(mu_bt + beta_bt' *x_stat);
  // 
    // // predictive intervals
  // randnorm = normal_rng(0,1);
  // p_init_pred = inv_logit(mu_bt + (tau + sigma) * randnorm);
  // p_x_pred = inv_logit(mu_bt + beta_bt' *x_pred + (tau + sigma) * randnorm);
  // p_xstat_pred = inv_logit(mu_bt + beta_bt' *x_stat + (tau + sigma) * randnorm);
  
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

