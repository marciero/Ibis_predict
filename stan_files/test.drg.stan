// test file to check recaptured parameter values
data {
    int<lower=1> N;
    real alpha;
    // real phi;
}
generated quantities {
  //  real age[N];
  //  real age_adj[N];
    int counts[N];
   // real alpha = normal_rng(log(4), 0.5);
  //  real beta = normal_rng(log(1.75), 0.2);
  //  real inv_phi = fabs(normal_rng(1, .5));

    for (n in 1:N) {
   // age[n] = normal_rng(70, 8);
  //  age_adj[n] = (age[n] - 70)/8;
    counts[n] =
             poisson_log_rng(alpha);
    }
}
