
data {
  int<lower=1> N;
  vector<lower=0>[N] traps;
  int<lower=0> complaints[N];
}
parameters {
  real alpha;
  real beta;
  real<lower=0> inv_phi;
}
transformed parameters {
  real phi = inv(inv_phi);
}
model {
  alpha ~ normal(log(4), 1);
  beta ~ normal(0.25, 1);
  inv_phi ~ normal(0, 1);

  complaints ~ neg_binomial_2_log(alpha + beta * traps , phi);
}

/*
generated quantities {
  int y_rep[N];
  for (n in 1:N)
    y_rep[n] = neg_binomial_2_log_rng(alpha + beta * traps[n], phi);

}
*/

