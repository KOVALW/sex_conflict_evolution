
// The input data
data {
  int<lower=0> N;
  int<lower=0> Ktheta;
  
  int<lower = 0> b_coefs;
  int<lower=0> eta_coefs;
  
  int offspring[N];
  int eta_id[N];
  
  int theta_id[N];
  int<lower=0,upper=1>mortality[N];
  int<lower=0,upper=1>zero_bool[N];
  
  matrix[N, b_coefs] tmt_spec;
  
}

transformed data {
  int eta_coefs1 = eta_coefs + 1;
  int K1 = Ktheta + 1;
}

// The parameters accepted by the model. 
parameters {
  
  real<lower=3,upper=5> base_reproduction;
  real<lower=0,upper=1> base_theta;
  real<lower=0,upper=1> mort_theta;
  
  vector<lower=0,upper=1>[Ktheta] tmt_theta;
  
  vector[b_coefs] beta_vals;
  vector[eta_coefs] etas;
  
}

transformed parameters {
  vector[eta_coefs1] alpha = rep_vector(base_reproduction, eta_coefs1);
  vector<lower=0, upper=1>[K1] theta_spec = rep_vector(base_theta, K1);
  vector<lower=0, upper=1>[N] theta;
  vector<lower=0, upper=1>[N] mort_fx = rep_vector(0, N);
  vector<upper=0>[N] log_theta;
  vector<upper=0>[N] log1m_theta;
  
  for (i in 2:eta_coefs1) {
    alpha[i] += etas[i-1];
  }
  
  for (i in 2:K1) {
    theta_spec[i] = tmt_theta[i-1];
  }
  
  for ( i in 1:N){
    mort_fx[i] += mort_theta * mortality[i];
    theta[i] = mort_fx[i] + (1 - mort_fx[i]) * theta_spec[theta_id[i]];
  }
  
  log_theta = log(theta);
  log1m_theta = log1m(theta);
  
}

// The model to be estimated. 
model {
  
  base_reproduction ~ normal(4, 0.25);
  base_theta ~ beta(2, 38);
  etas ~ normal(0, 5);
  beta_vals ~ normal(0, 5);
  theta ~ uniform(0, 1);
  
  target += poisson_log_glm_lpmf(offspring | tmt_spec, alpha[eta_id], beta_vals) + log1m_theta;
  
}

generated quantities {
  vector[N] log_lik;
  
  log_lik = poisson_log_glm_lpmf(offspring | tmt_spec, alpha[eta_id], beta_vals) + log1m_theta;
  
}
