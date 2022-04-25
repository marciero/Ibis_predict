data {
   int N;
   int hcc_count[N];
   }
   parameters {
   real<lower = 0> alpha;
   }
   model {
    alpha ~ normal(log(3), 1);
    hcc_count ~ poisson_log(alpha);
      }
