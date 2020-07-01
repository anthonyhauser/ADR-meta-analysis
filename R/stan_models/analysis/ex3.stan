// Stan model for simple linear regression
data {
  int < lower = 0 > N; //number of observations
  int < lower = 1 > M_y; // number of y variables
  int < lower = 1 > M_xy; // number of x and y variables
  
  //vector[N] x; //variable with missing entry
  matrix[N,M_y] y; //other variables (possibly of length 0)
  //matrix[N,M_xy] xy; //all variables
  
  real x_nna [N-3]; //x in non-missing position
  int na_time_pos[3]; //position of missing x
  int nna_time_pos[N-3]; //position of non-missing x
  
  vector[N] k; //number of sucesses
  vector[N] n; //number of trials
  
  int inference;
}

parameters {
  real beta_x; //regression coefficient of x
  vector[M_y] beta_y; //regression coefficient of y
  vector[M_xy] beta_xy; //regression coefficient of x and y
  vector[3] x_imput; //3 missing values
  
  real < lower = 0 > sigma;
  vector[N] u;
  
  real mu;
}

transformed parameters {
  //vector[N] x;
  matrix[N,M_xy] xy;
  
  xy[nna_x_pos,1] = x_nna;
  xy[na_x_pos,1] = x_imput;
  if(M_xy>1){
    for(i in 2:M_xy){
      xy[,i] = y[,i-1]
    }
  }
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
  target += binomial_lpmf(k|n,inv_logit(mu + sigma * u + xy * beta_xy));
}


generated quantities {
  
} // The prior or posterior predictive distribution

