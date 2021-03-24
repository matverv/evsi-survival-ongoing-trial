// Exponential survival model
// Modified from the R package "survHE" by Gianluca Baio. https://CRAN.R-project.org/package=survHE

functions {
  // Defines the log hazard
  vector log_h (vector t, real rate) {
    vector[num_elements(t)] v_log_h;
    for (i in 1:num_elements(t)) {
      v_log_h[i] = log(rate);
    }
    return v_log_h;
  }
  
  // Defines the log survival
  vector log_S (vector t, real rate) {
    vector[num_elements(t)] v_log_S;
    for (i in 1:num_elements(t)) {
      v_log_S[i] = -rate * t[i];
    }
    return v_log_S;
  }
  
  // Defines the sampling distribution
  real surv_exponential_lpdf (vector t, vector d, vector LL, real rate) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t, rate) + log_S(t, rate) - log_S(LL, rate);
    prob = sum(log_lik);
    return prob;
  }
}

data {
  int n;                  // number of observations
  vector[n] t;            // observed times
  vector[n] d;            // censoring indicator (1=observed, 0=censored)
  vector<lower=0>[n] LL;  // lower truncation limit
  int<lower=0> events0;   // prior number of events
  int<lower=0> trisk0;    // prior time at risk
}

parameters {
  real<lower=0> rate;     // rate
}

model {
  // Likelihood
  target += surv_exponential_lpdf(t | d, LL, rate);

	// Priors
	target += gamma_lpdf(rate | events0, trisk0); 
}

