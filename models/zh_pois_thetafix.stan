
// The input data
data {
  int<lower=0> N;
  
  int<lower = 0> b_coefs;
  int<lower=0> eta_coefs;
  
  int offspring[N];
  int eta_id[N];
  
  int<lower=0,upper=1>mortality[N];
  int<lower=0,upper=1>zero_bool[N];
  
  matrix[N, b_coefs] tmt_spec;
  
}

transformed data {
  int eta_coefs1 = eta_coefs + 1;
  int zeroes[N] = rep_array(0, N);
}

// The parameters accepted by the model. 
parameters {
  
  real<lower=3,upper=5> base_reproduction;
  real<lower=0,upper=1> base_theta;
  real<lower=0,upper=1> mort_theta;
  
  vector[b_coefs] beta_vals;
  vector[eta_coefs] etas;
  
}

transformed parameters {
  vector[eta_coefs1] alpha = rep_vector(base_reproduction, eta_coefs1);
  vector<lower=0, upper=1>[N] theta;
  real<lower = 0, upper = 1> mort_fx;
  //vector<lower=0, upper=1>[N] mort_fx = rep_vector(0, N);
  
  for (i in 2:eta_coefs1) {
    alpha[i] += etas[i-1];
  }
  
  //theta = mort_theta * mortality + (1 - mortality) * base_theta;
  
  for ( i in 1:N){
    mort_fx = mort_theta * mortality[i];
    theta[i] = mort_fx + (1 - mort_fx) * base_theta;
  }
  
}

// The model to be estimated. 
model {
  
  base_reproduction ~ normal(3.9, 0.3);
  base_theta ~ beta(2, 38);
  mort_theta ~ beta(2, 38);
  etas ~ normal(0, 5);
  beta_vals ~ normal(0, 5);
  
  target += poisson_log_glm_lpmf(offspring | tmt_spec, alpha[eta_id], beta_vals) - 
            log1m_exp(poisson_log_glm_lpmf(zeroes | tmt_spec, alpha[eta_id], beta_vals));
  
  target += bernoulli_lpmf(zero_bool | theta);
  
}

generated quantities {
  vector[N] log_lik;
  real beta_fx_row;
  
  for (i in 1:N)
  {
    beta_fx_row = alpha[eta_id[i]] + dot_product( beta_vals,  tmt_spec[i,]);
    
    log_lik[i] = bernoulli_lpmf(zero_bool[i] | theta[i]) +
    poisson_log_lpmf(offspring[i] | beta_fx_row) - 
    log1m_exp(poisson_log_lpmf(0 | beta_fx_row));
  }
}
