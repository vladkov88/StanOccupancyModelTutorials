data {
  int<lower=0> nSites;
  int<lower=0> nSurveys;
  int<lower=0,upper=1> y[nSites,nSurveys];
}
parameters {
  real muPsi;
  real muP;
}
model {

  real log_muPsi;
  real log1m_muPsi;

  log_muPsi   = log_inv_logit(muPsi);
  log1m_muPsi = log1m_inv_logit(muPsi);

  muP ~ normal(0, 2);
  muPsi ~ normal(0, 2);
  
  // likelihood
  for (r in 1:nSites) {
    if (sum(y[r]) > 0)
      target +=
	log_muPsi + bernoulli_logit_lpmf(y[r] | muP);
    else
      target +=
	log_sum_exp(log_muPsi +
		    bernoulli_logit_lpmf(y[r] | muP), log1m_muPsi);
  }
}
generated quantities{
  real<lower = 0, upper = 1> p;
  real<lower = 0, upper = 1> psi;

  p   = inv_logit(muP);
  psi = inv_logit(muPsi);
}
