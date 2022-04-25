 // cost_01 model with priors only

   generated quantities {
          int visits[200];
          int hcc_count[200];
          int char_score[200];
          real age_norm[200];
         real cost[200];
         real alpha  = normal_rng(15, 4);
          real beta_vs = normal_rng(1, 1);
      real beta_hcc = normal_rng(2, 2);
         real beta_char  = normal_rng(3, 1);
          real beta_age = normal_rng(2, 1);
        real<lower = 0> sigma = exponential_rng(1);


       for (i in 1:200){
       visits[i] = poisson_rng(2.5);
       hcc_count[i] = neg_binomial_2_rng(3.5, 3.2);
       char_score[i] = binomial_rng(8, 0.25);
        age_norm[i] = normal_rng(0, 1);
        cost[i] = normal_rng(alpha + beta_vs * visits[i] +
      beta_hcc * hcc_count[i] +
              beta_char * char_score[i]  + beta_age * age_norm[i], sigma);
         }
     }
