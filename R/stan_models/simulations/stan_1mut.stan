// Stan model for simple linear regression
data {
      int < lower = 1 > N; // number of studies
      int < lower = 1 > B; // number of co-variates
  
      int k[N]; // Number of events by study
      int n[N]; // Number of people at risk by study
      
      matrix[N,B] x; //matrix of dim N times B, B co-variates for each of the N (sub-)study

      real < lower = 0 > init_res_1;
      real < lower = 0 > init_res_2;
      real < lower = 0 > sigma_1;
      vector[B] beta_1;
      vector < lower = 0 > [B] beta_2;

      int inference;
}
      
parameters {
      real < lower = 0, upper = 1 > init_res; 
      real < lower = 0 > sigma;
      vector[B] beta_tilde;
      real alpha_tilde;
}

transformed parameters {
      real mu = logit(init_res);
      real alpha = mu + sigma * alpha_tilde;
      vector[B] beta;
      for(j in 1:B){
        beta[j] = beta_1[j] + beta_2[j] * beta_tilde[j];
      }
}
      
model {
      //Priors
      init_res ~ beta(init_res_1,init_res_2);
      sigma ~ exponential(sigma_1);
      beta_tilde ~ normal(0,1);
      
      //Random intercept
      alpha_tilde ~ normal(0,1);
      
      //Likelihood
      if(inference==1){
        for(i in 1:N){
          target += binomial_lpmf(k[i]|n[i],inv_logit(alpha+x[i,]*beta));
        }
      }
}
      
generated quantities {
      real p_pred_min = inv_logit(normal_rng(mu,sigma));
      real p_pred_max = inv_logit(normal_rng(mu+beta[1]*50/500+beta[2]*100/500,sigma));
} // The prior or posterior predictive distribution

