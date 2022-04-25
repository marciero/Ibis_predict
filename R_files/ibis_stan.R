vis_01 <- "
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
   vector[N]  hcc_count;
   }
   parameters {
   real<lower = 0> alpha;
   real beta;
 //  real<lower = 0> sigma_alpha;
  // real<lower = 0> sigma_beta;
   }
   model {
    alpha ~ normal(log(4.6), 1);
    beta ~ normal(1, 1);
    visits ~ poisson_log(alpha + beta * hcc_count);
   }
      generated quantities {
       int y_rep[N];
       for (i in 1:N)
    y_rep[i] = poisson_log_safe_rng(alpha +  beta * hcc_count[i]);
    }
"
vis_01_sim <- stan(model_code = vis_01, data = list(N = nrow(augusta_unique_hcc),
                                                    visits = augusta_unique_hcc$visits,
                                                              hcc_count = augusta_unique_hcc$hcc_conds),
                     chains = 4, iter = 1000, seed = 84732
)

tidy(vis_01_sim) %>% head(2)
vis_01_samps <- extract(vis_01_sim)

samps_df <- as.data.frame(vis_01_sim)
pred_new <- samps_df %>% mutate(vis_5cond = rpois(2000, exp(alpha + 5*beta))) %>%
         mutate(vis_10cond = rpois(2000, exp(alpha + 10*beta))) %>% select(c(vis_10cond, vis_5cond))

library(rethinking)
precis(samps_df %>% select(c(alpha, beta)))
precis(pred_new)

mcmc_hist(pred_new)

mean(pred_new$vis_10cond > 4)/mean(pred_new$vis_5cond > 4)

y_rep <- as.matrix(vis_01_sim, pars = "y_rep")
ppc_dens_overlay(y = augusta_unique_hcc$visits, yrep = y_rep[1:200, ]) +
    labs(x = "visits") + xlim(c(0, 20))


############
############
############  Above with negative binomial (not included in ibis_01)

vis_02 <- "
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
   }
   parameters {
   real<lower = 0> alpha;
   real beta;
   real<lower = 0> inv_phi;
   }

  transformed parameters {
     real<lower = 0> phi = inv(inv_phi);
   }
   model {
    alpha ~ normal(log(4.6), 1);
    beta ~ normal(1, 1);
    inv_phi ~ exponential(1);
    visits ~ neg_binomial_2_log(alpha + beta * hcc_count, phi);
   }
      generated quantities {
       int y_rep[N];
       for (i in 1:N)
    y_rep[i] = neg_binomial_2_log_safe_rng(alpha +  beta * hcc_count[i], phi);
    }
"

vis_02_sim <- stan(model_code = vis_02, data = list(N = nrow(augusta_unique_hcc),
                                                    visits = augusta_unique_hcc$visits,
                                                    hcc_count = augusta_unique_hcc$hcc_conds),
                   chains = 4, iter = 1000, seed = 84732
)

#tidy(vis_02_sim) %>% head(2)
vis_02_samps <- extract(vis_02_sim)

samps_nb_df <- as.data.frame(vis_02_sim)

## problem here. need to resolve phi vs probability
#pred_new_nb <- samps_nb_df %>% mutate(vis_5cond = rpois(2000, exp(alpha + 5*beta))) %>%
    mutate(vis_10cond = rpois(2000, exp(alpha + 10*beta))) %>% select(c(vis_10cond, vis_5cond))

library(rethinking)
precis(samps_nb_df %>% select(c(alpha, beta, phi)))

mcmc_dens(samps_nb_df %>% select(c(alpha, beta, phi)))
mcmc_hist(pred_new_nb)

precis(pred_new_nb)

## problem here
mean(pred_new_nb$vis_10cond > 4)/mean(pred_new_nb$vis_5cond > 4)

y_rep <- as.matrix(vis_02_sim, pars = "y_rep")
ppc_dens_overlay(y = augusta_unique_hcc$visits, yrep = y_rep[1:200, ]) +
    labs(x = "visits") + xlim(c(0, 20))

### various other ways access parameters, y_reps
samps_nb_df$`y_rep[1]`
vis_02_samps <- extract(vis_02_sim)
as.data.frame(vis_02_sim, pars = "y_rep") %>% View()

##############
##############  Now with NB using ibis.

vis_03 <- "
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
       for (i in 1:N) {
    y_rep[i] = neg_binomial_2_log_safe_rng(alpha +  beta[ibis_idx[i]] * hcc_count[i], phi);
       }
      ibis_no = neg_binomial_2_log_safe_rng(alpha +  5*beta[1], phi);
       ibis_yes = neg_binomial_2_log_safe_rng(alpha + 5* beta[2], phi);
      }
"
df <- augusta_unique_hcc
ibis_fake <- df %>% mutate(wt = case_when(visits < 5 ~ 0.65,
                                          visits >= 5 ~ 0.45)) %>%
    group_by(1:n()) %>%
    mutate(ibis = sample(c(0, 1), size = 1, replace = TRUE, prob =  c(1-wt, wt))) %>%
    mutate(ibis_idx = ibis +1) %>% ungroup()

vis_03_sim <- stan(model_code = vis_03, data = list(N = nrow(ibis_fake),
                                                    visits = ibis_fake$visits,
                                                    hcc_count = ibis_fake$hcc_conds,
                                                  ibis_idx = ibis_fake$ibis_idx),
                   chains = 4, iter = 1000, seed = 84732
)

#tidy(vis_02_sim) %>% head()
vis_03_samps <- extract(vis_03_sim)

#library(rethinking)
samps_ibis_df <- as.data.frame(vis_03_sim)
precis(samps_ibis_df %>% select(c(alpha, `beta[1]`, `beta[2]`, phi)), depth = 2)

mcmc_dens(samps_ibis_df %>% select(c(alpha, `beta[1]`, `beta[2]`, phi)))

y_rep <- as.matrix(vis_03_sim, pars = "y_rep")
ppc_dens_overlay(y = augusta_unique_hcc$visits, yrep = y_rep[1:200, ]) +
    labs(x = "visits") + xlim(c(0, 20))

#### Need to check mean vs probability calling rnbinom from R vs stan
pred_new_ibis <- samps_ibis_df %>% mutate(ibis_no = rnbinom(2000, phi, 1/(1 + exp(alpha + 5*`beta[1]`))),
                                          ibis_yes = rnbinom(2000, phi, 1/(1 + exp(alpha + 5*`beta[2]`)))
                                              ) %>%
     select(c(ibis_yes, ibis_no))

mean(pred_new_ibis$ibis_no > 3)/mean(pred_new_ibis$ibis_yes > 3)
## Should check these against what we get generating them in the stan file:
mean(vis_03_samps$ibis_no > 3)/mean(vis_03_samps$ibis_yes > 3)

### very few zero visits, but y_reps has a lot of them...
as.data.frame(vis_03_samps$y_rep) %>%  map_dbl(~mean(.x == 0)) %>% mean()
## again these are for five hcc's.
## Should check against what we get generating them in the stan file


precis(pred_new_ibis)
mcmc_hist(pred_new_ibis)
mean(pred_new_ibis$vis_10cond > 4)/mean(pred_new_nb$vis_5cond > 4)
