// Stan model for simple linear regression
functions {
  real min_2real(real x, real y) {
    vector[2] v;
    v[1]=x;
    v[2]=y;
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
  real beta_16;
  
  real beta_21;
  real beta_22;
  real beta_26;
  real beta_23;
  real beta_24;
  real beta_25;
  
  real beta_31;
  real beta_34;
  real beta_32;
  real beta_33; //<lower = beta_32>
  real beta_35; //<lower = beta_32>
  real beta_36;
    
  real beta_41;
  real beta_43;
  real <lower = beta_43> beta_42;
  real beta_45;
  real beta_44;
  real beta_46;
      
  real beta_51;
  real beta_53;
  real <lower = beta_53> beta_52;
  real beta_55;
  real beta_54;
  real beta_56;
      
  real beta_61;
  real beta_63;
  real <lower = beta_63> beta_62;
  real beta_64;
  real <upper = beta_64> beta_65; //<lower = beta_62>
  real <upper = beta_64> beta_66;
        
  real beta_71;
  real beta_72;
  real beta_73;
  real beta_76;
  real <upper = beta_76> beta_74; //<lower = beta_62>
  real <upper = beta_76> beta_75;
        
  real beta_81;
  real beta_83;
  real <lower = beta_83> beta_82;
  real beta_86;
  real <lower = beta_86> beta_84;
  real <lower = beta_86> beta_85;
        
        
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
  beta_16 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  beta_21 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_22 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_23 ~ normal_mixture(prob_bernoulli,beta_3l[1],0,beta_3l[2], epsilon);
  beta_24 ~ normal_mixture(prob_bernoulli,beta_3l[1],0,beta_3l[2], epsilon);
  beta_25 ~ normal_mixture(prob_bernoulli,beta_3l[1],0,beta_3l[2], epsilon);
  beta_26 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  beta_31 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_32 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_33 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_34 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_35 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_36 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  beta_41 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_42 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_43 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_44 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_45 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_46 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  beta_51 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_52 ~ normal_mixture(prob_bernoulli,beta_4l[1],0,beta_4l[2], epsilon);
  beta_53 ~ normal_mixture(prob_bernoulli,beta_4l[1],0,beta_4l[2], epsilon);
  beta_54 ~ normal_mixture(prob_bernoulli,beta_4l[1],0,beta_4l[2], epsilon);
  beta_55 ~ normal_mixture(prob_bernoulli,beta_4l[1],0,beta_4l[2], epsilon);
  beta_56 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  beta_61 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_62 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_63 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_64 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_65 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_66 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  beta_71 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_72 ~ normal_mixture(prob_bernoulli,beta_3l[1],0,beta_3l[2], epsilon);
  beta_73 ~ normal_mixture(prob_bernoulli,beta_3l[1],0,beta_3l[2], epsilon);
  beta_74 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_75 ~ normal_mixture(prob_bernoulli,beta_3l[1],0,beta_3l[2], epsilon);
  beta_76 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  beta_81 ~ normal_mixture(prob_bernoulli,beta_time[1],0,beta_time[2], epsilon);
  beta_82 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_83 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_84 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_85 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  beta_86 ~ normal_mixture(prob_bernoulli,beta_p1,0,beta_p2, epsilon);
  
  
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
} // The prior or posterior predictive distribution

