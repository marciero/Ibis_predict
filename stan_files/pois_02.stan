//  See visits.R
// visits ~ Pois( alpha_j + beta_j* hcc)


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
   int ibis_idx[N];
   }
  // transformed data {
   //   real<lower = 0> lambda = mean(visits);
  // }
   parameters {
   vector[2] alpha;
   vector[2] beta;
 // real<lower = 0> sigma_alpha;
  //  real<lower = 0> sigma_beta;

   }
   model {
       vector[N] lambda;
   //    sigma_alpha ~ exponential(1);
   //    sigma_beta ~ exponential(1);
    alpha ~ normal(log(3.5), 1);
     beta ~ normal(0, 1);
    for (i in 1:N) {
        lambda[i] = alpha[ibis_idx[i]] + beta[ibis_idx[i]] * hcc_conds[i];
    }
    visits ~ poisson_log(lambda);
   }
/*
 generated quantities {
       int y_rep[N];
       for (i in 1:N)
    y_rep[i] =
    poisson_log_safe_rng(alpha[ibis_idx[i]] + beta[ibis_idx[i]] * visits[i]);
    }
*/

