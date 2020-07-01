// Stan model for simple linear regression
data {
  int < lower = 1 > N; // number of studies with TAMs
  int M; // Number of TAMs: M=6
  int k[N]; // Number of people with any TAMs by study
  int n[N]; // Number of people at risk by study
  
  real p_noTAM_min[N]; //minimum prevalence of no TAM, i.e. given TAMs are independant
  real p_noTAM_max[N]; //maximum prevalence of no TAM, i.e. given TAMs are highly correlated
  real p_FTC_ZDV_min; //min prevalence of no TAM, product of (1-estimated prevalence with 100% FTC and ZDV, after 3 years)
  real p_FTC_TDF_min; //min prevalence of no TAM, product of (1-estimated prevalence with 100% FTC and TDF, after 3 years)
  real p_FTC_ZDV_max; //max prevalence of no TAM, min of (1 - estimated prevalence with 100% FTC and ZDV, after 3 years)
  real p_FTC_TDF_max; //max prevalence of no TAM, min of (1 - estimated prevalence with 100% FTC and TDF, after 3 years)
  
  int inference;
}

parameters {
  real < lower = 0, upper = 1 > alpha;
}

transformed parameters {
  //P( any TAM ) = 1- P( no TAM )
  real p_TAM[N];

  for(i in 1:N){
    p_TAM[i] = 1 - ( p_noTAM_min[i] + alpha * (p_noTAM_max[i] - p_noTAM_min[i]));
  }
}

model {
  //Priors
  alpha  ~ uniform(0,1);
  
  //Likelihood
  if(inference==1){
    for(i in 1:N){
      target += binomial_lpmf(k[i]|n[i], p_TAM[i]);
    }
  }
}

generated quantities {
  real p_FTC_ZDV;
  real p_FTC_TDF;
  
  p_FTC_ZDV = 1 - ( p_FTC_ZDV_min + alpha * (p_FTC_ZDV_max - p_FTC_ZDV_min));
  p_FTC_TDF = 1 - ( p_FTC_TDF_min + alpha * (p_FTC_TDF_max - p_FTC_TDF_min));
} // The prior or posterior predictive distribution
