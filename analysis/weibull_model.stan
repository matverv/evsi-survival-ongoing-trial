// Weibull survival model
// Modified from the R package "survHE" by Gianluca Baio. https://CRAN.R-project.org/package=survHE

functions {
  // Defines the log hazard
  vector log_h (vector t, real shape, real scale) {
    vector[num_elements(t)] v_log_h;
    v_log_h = log(shape)+(shape-1)*log(t ./ scale)-log(scale);
    return v_log_h;
  }
  
  // Defines the log survival
  vector log_S (vector t, real shape, real scale) {
    vector[num_elements(t)] v_log_S;
    for (i in 1:num_elements(t)) {
      v_log_S[i] = -pow((t[i]/scale),shape);
    }
    return v_log_S;
  }
  
  // Defines the sampling distribution
  real surv_weibull_lpdf (vector t, vector d, vector LL, real shape, real scale) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t, shape, scale) + log_S(t, shape, scale) - log_S(LL, shape, scale);
    prob = sum(log_lik);
    return prob;
  }
}

data {
  int n;                  // number of observations
  vector[n] t;            // observed times
  vector[n] d;            // censoring indicator (1=observed, 0=censored)
  int<lower=1> K;         // number of model parameters
  vector[K] theta0;       // prior means 
	cov_matrix[K] sigma0;   // prior variance-covariance matrix
	vector<lower=0>[n] LL;  // lower truncation limit

}

parameters {
	vector[K] theta;                                                                                  
}

transformed parameters {
  real shape = exp(theta[1]);
  real scale = exp(theta[2]);
}

model {
  // Likelihood
  target += surv_weibull_lpdf(t | d, LL, shape, scale);

  // Priors
	target += multi_normal_lpdf(theta | theta0, sigma0);
}
