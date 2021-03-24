// log-Normal survival model
// Modified from the R package "survHE" by Gianluca Baio. https://CRAN.R-project.org/package=survHE

functions {
  // Defines the log survival
  vector log_S (vector t, real mean, real sd) {
    vector[num_elements(t)] v_log_S;
    for (i in 1:num_elements(t)) {
      v_log_S[i] = log(1-Phi((log(t[i])-mean)/sd));
    }
    return v_log_S;
  }
  
  // Defines the log hazard
  vector log_h (vector t, real mean, real sd) {
    vector[num_elements(t)] v_log_h;
    vector[num_elements(t)] ls;
    ls = log_S(t,mean,sd);
    for (i in 1:num_elements(t)) {
      v_log_h[i] = lognormal_lpdf(t[i]|mean, sd) - ls[i];
    }
    return v_log_h;
  }

  // Defines the sampling distribution
  real surv_lognormal_lpdf (vector t, vector d, vector LL, real mean, real sd) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t, mean, sd) + log_S(t, mean, sd) - log_S(LL, mean, sd);
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
  real meanlog = theta[1];
  real sdlog = exp(theta[2]); //exp(log_sigma);
}

model {
  // Likelihood
  target += surv_lognormal_lpdf(t | d, LL, meanlog, sdlog);

  // Priors
	target += multi_normal_lpdf(theta | theta0, sigma0);
}
