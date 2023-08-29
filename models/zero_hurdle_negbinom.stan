
// The input data
data {
  int<lower=0> N;
  int<lower=0> NZ;
  int<lower = 0> b_coefs;
  int<lower=0> eta_coefs;
  
  int offspring[N];
  int eta_id[N];
  
  matrix[N, b_coefs] tmt_spec;
  
}

transformed data {
  int eta_coefs1 = eta_coefs + 1;
  int zeroes[N] = rep_array(0, N);
}

// The parameters accepted by the model. 
parameters {
  
  real<lower=3,upper=5> base_reproduction;
  real<lower=0, upper=1> theta;
  real<lower=0, upper=100> phi;
  
  vector[b_coefs] beta_vals;
  vector[eta_coefs] etas;
  
}

transformed parameters {
  vector[eta_coefs1] alpha = rep_vector(base_reproduction, eta_coefs1);
  
  for (i in 2:eta_coefs1) {
    alpha[i] += etas[i-1];
  }
  
}

// The model to be estimated. 
model {
  
  base_reproduction ~ normal(4, 0.25);
  theta ~ uniform(0,1);
  etas ~ normal(0, 5);
  beta_vals ~ normal(0, 5);
  phi ~ uniform(0,100);
  
  target += neg_binomial_2_log_glm_lpmf(offspring | tmt_spec, alpha[eta_id], beta_vals, phi) - 
            log1m_exp(neg_binomial_2_log_glm_lpmf(zeroes | tmt_spec, alpha[eta_id], beta_vals, phi));
  
  target += NZ * bernoulli_lpmf(1 | theta) + N * bernoulli_lpmf(0 | theta);
  
}
