//
// truncated poisson
functions {
  /*
  * Alternative to poisson_log_rng() that
  * avoids potential numerical problems during warmup
  */
  int poisson_log_safe_rng(real eta) {
    real pois_rate = exp(eta);
    if (pois_rate >= exp(20.79))
      return -9;
    return poisson_rng(pois_rate);
  }
}

 data {
   int N;
   int visits[N];
   vector[N] hcc_conds;
   }
  // transformed data {
   //   real<lower = 0> lambda = mean(visits);
  // }
   parameters {
   real<lower = 0> alpha;
   real beta;
   }
   model {
    alpha ~ normal(log(3.5), 1);
    beta ~ normal(0.5, 1);
    visits ~ poisson_log(alpha + beta*hcc_conds);
   }

 generated quantities {
       int y_rep[N];
       for (i in 1:N)
    y_rep[i] = poisson_log_safe_rng(alpha + beta*hcc_conds[i]);
    }
