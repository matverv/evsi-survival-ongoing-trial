// Gamma survival model

data {
  int n;                          // number of observations
  int<lower=0> nobs;              // number of observed events
  int<lower=0> ncens;             // number of censored events
  vector<lower=0>[nobs] yobs;     // observed event times
  vector<lower=0>[ncens] ycens;   // censoring times
  int<lower=1> K;                 // number of model parameters
  vector[K] theta0;               // prior means 
  cov_matrix[K] sigma0;           // prior variance-covariance matrix
  vector<lower=0>[n] LL;          // lower truncation limit
}

parameters {
  vector[K] theta;                                                                                  
}

transformed parameters {
  real shape = exp(theta[1]);
  real rate = exp(theta[2]);
}

model {
  // Likelihood
  target += gamma_lpdf(yobs | shape, rate);
  target += gamma_lccdf(ycens | shape, rate);
  target += -gamma_lccdf(LL | shape, rate); // subtract the log complementary cdf given the truncation times from the log likelihood. See: https://discourse.mc-stan.org/t/reparameterizing-a-truncated-normal/3522/13
  
  // Priors
  target += multi_normal_lpdf(theta | theta0, sigma0);
}

