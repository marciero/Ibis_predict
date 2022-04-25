// Stan shell file to test code to be used inline within Rmd file ibis_01, to replace the rstanarm
   functions {
  /*
  * Alternative to neg_binomial_2_log_rng() that
  * avoids potential numerical problems during warmup
  */
  int neg_binomial_2_log_safe_rng(real eta, real phi) {
    real gamma_rate = gamma_rng(phi, phi / exp(eta));
    if (gamma_rate >= exp(20.79))
      return -9;

    return poisson_rng(gamma_rate);
  }
}
 data {
   int N;
   int visits[N];
   vector[N]  hcc_count;
   int  ibis_idx[N];
   }
   parameters {
   real<lower = 0> alpha;
   vector[2] beta;
   real<lower = 0> inv_phi;
   }

  transformed parameters {
     real<lower = 0> phi = inv(inv_phi);
   }
   model {
    alpha ~ normal(log(4.6), 1);
    beta ~ normal(1, 1);
    inv_phi ~ exponential(1);
    visits ~ neg_binomial_2_log(alpha + beta[ibis_idx] .* hcc_count, phi);
   }
      generated quantities {
       int ibis_no;
       int ibis_yes;
       int y_rep[N];
       for (i in 1:N){
    y_rep[i] = neg_binomial_2_log_safe_rng(alpha +  beta[ibis_idx[i]] * hcc_count[i], phi);
     }
       ibis_no = neg_binomial_2_log_safe_rng(alpha +  5*beta[1], phi);
       ibis_yes = neg_binomial_2_log_safe_rng(alpha + 5* beta[2], phi);
      }
