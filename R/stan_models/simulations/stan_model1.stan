// Stan model for simple linear regression

      data {
      int < lower = 1 > N; // Sample size
      vector[N] x; // Predictor
      vector[N] y; // Outcome
      }
      
      parameters {
      real alpha; // Intercept
      real beta; // Slope (regression coefficients)
      real < lower = 0 > sigma; // Error SD
      }
      
      model {
      y ~ normal(alpha + x * beta , sigma);
      }
      
      generated quantities {
      } // The posterior predictive distribution
