 data {
   int N;
   real cost[N];
   vector[N] visits;
   vector[N] hcc_count;
   vector[N]  char_score;
   vector[N] age_norm;
  /*
      real elix_score;
   real gender_idx;
   real fin_idx;
    */
   }
/*transformed data {
    vector[N] age_norm;
    for (i in 1:N)
    age_norm[i] = (age[i] - mean(age))/sd(age);
}
*/
   parameters {
   real<lower = 0> alpha;
   real beta_vs;
   real beta_hcc;
   real beta_char;
   real beta_age;
   real<lower = 0> sigma;
   }
 transformed parameters {
     vector[N] mu = alpha + beta_vs * visits +
      beta_hcc * hcc_count +
              beta_char * char_score  + beta_age * age_norm;
   }
   model {
    alpha ~ normal(10, 5);
    beta_vs ~ normal(1, 4);
    beta_hcc ~ normal(1, 4);
    beta_char ~ normal(1, 4);
    beta_age ~ normal(1, 4);  // for every sd in age...
   // visits ~ normal(1, 1);  // was mistake but ran anyway-was not used
    sigma ~ exponential(1);
    cost ~ normal(mu, sigma);
   }
      generated quantities {
          real y_rep[N];
       for (i in 1:N){
    y_rep[i] = normal_rng(mu[i], sigma);
     }
      }
