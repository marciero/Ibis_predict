##// test file to check recaptured parameter values
dgp_cp <- stan_model("test.drg.stan")
dgp_sim <- sampling(dgp_cp, data = list(N = 100, alpha = log(3)), chains = 1, iter = 1,
                    algorithm = 'Fixed_param')

samps <- rstan::extract(dgp_sim)
tidy(dgp_sim)
samps_df <- as.data.frame(dgp_sim)
samps_df %>% ggplot(aes(`counts[1]`)) + geom_histogram()
summary(dgp_sim)

test_model <- "
   data {
   int N;
   int counts[N];
   }
   parameters {
   real<lower = 0> alpha;
   }
   model {
    alpha ~ normal(log(5), 1);
    counts ~ poisson_log(alpha);
   }
      generated quantities {
       int y_rep[N];
       for (i in 1:N)
    y_rep[i] = poisson_log_rng(alpha);
    }
"
dat <- list(N = length(samps$counts), counts = samps$counts[1, ])
test_sim <- stan(model_code = test_model, data = dat,
                     chains = 4, iter = 1000
)

samps_mod <- rstan::extract(test_sim)
tidy(test_sim)
samps_df <- as.data.frame(test_sim)
samps_df %>% ggplot(aes(`counts[1]`)) + geom_histogram()
summary(test_sim)

y_rep <- as.matrix(test_sim, pars = "y_rep")
ppc_dens_overlay(y = samps$counts[1, ], yrep = y_rep[1:20, ])


hcc_pois_sim %>% as.data.frame() %>% View()
posterior <- extract(hcc_pois_sim)

## To use ppc_ we need to generate N y_reps in the model. Above we had 100 for simplicity.


prop <- function(x) mean(x == 2)
ppc_stat(y = augusta_unique_hcc$hcc_conds, yrep = y_rep, stat = "prop")




