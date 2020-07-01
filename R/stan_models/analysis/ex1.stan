// Stan model for simple linear regression
data {
  int < lower = 0 > N; //number of observations
  int < lower = 1 > M_y; // number of y variables
  int < lower = 1 > M_xy; // number of x and y variables
  
  vector[N] x; //variable with missing entry
  matrix[N,M_y] y; //other variables (possibly of length 0)
  
  vector[N-3] x_nna; //x in non-missing position
  int na_x_pos[3]; //position of missing x
  int nna_x_pos[N-3]; //position of non-missing x
  
  int k [N]; //number of sucesses
  int n [N]; //number of trials
  
  int inference;
}

parameters {
  vector[M_xy] beta_xy; //regression coefficient of x and y
  vector[3] x_imput; //3 missing values
  
  real < lower = 0 > sigma;
  vector[N] u;
  
  real mu;
}

transformed parameters {
  //vector[N] x;
  matrix[N,M_xy] xy;
  vector[N] alpha;
  
  xy[nna_x_pos,1] = x_nna;
  xy[na_x_pos,1] = x_imput;
  if(M_xy>1){
    for(i in 2:M_xy){
      xy[,i] = y[,i-1];
    }
  }
  alpha = inv_logit(mu + sigma * u + xy * beta_xy);
}

model {
  //Priors
  mu ~ normal(0,5);
  beta_xy ~ normal(0,5);
  sigma ~ exponential(1);
  //Random effect
  for(i in 1:N){
    u[i] ~ normal(0,1);
  }
  //Imputation
  x_imput[1:3] ~ normal(2,2);
  
  //Likelihood
  if(inference==1){
    for(i in 1:N){
      target += binomial_lpmf(k[i]|n[i],alpha[i]);  
    }
  }
}


generated quantities {
  
} // The prior or posterior predictive distribution

