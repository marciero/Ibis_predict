library(comorbidity)

## See update on this package https://ellessenne.github.io/comorbidity/articles/C-changes.html

### uses hcc_augusta.R
df <- unnest(augusta_unique_hcc, data)

    elix <- comorbidity(df, id = "PATIENT_NBR", code = "code", map = "elixhauser_icd10_quan", assign0 = FALSE)
    char <- comorbidity(df, id = "PATIENT_NBR", code = "code", map = "charlson_icd10_quan", assign0 = FALSE)
comorbidity(df, id = "PATIENT_NBR", code = "code", map = "charlson_icd10_quan", assign0 = FALSE)

morb_idx <- cbind(elix, elix_score = score(elix, weights = "swiss", assign0 = T)) %>% full_join(
cbind(char, char_score = score(char, assign0 = T)), by = "PATIENT_NBR") %>%
    select(PATIENT_NBR, char_score, elix_score)

cor(morb_idx$char_score, morb_idx$elix_score)

df <- augusta_unique_hcc %>% full_join(morb_idx, by = "PATIENT_NBR")

df %>% filter(is.na(char_score))

is.na(df$PATIENT_NBR)



cor(df$visits, df$char_score)

df %>% ggplot(aes(char_score)) + geom_histogram()

df %>% ggplot() + geom_boxplot(aes(x = as.factor(char_score), y = visits))


####
#### Now use to simulate fake cost data. Then use model in cost_01.stan

## Try without age first, try with smaller data set
cost_df <- df %>% mutate(cost = rnorm(nrow(df), 15 + visits + 2*hcc_conds + 3*char_score, 10))
glimpse(cost_df)

cost_df_spl <- slice_sample(cost_df, n = 200, replace = FALSE)


cost_fake = list(cost = cost_df_spl$cost, N = nrow(cost_df_spl), visits = cost_df_spl$visits,
                 hcc_count = cost_df_spl$hcc_conds, char_score = cost_df_spl$char_score)

 ## cost_cp <- stan_model("cost_01.stan")

cost_sim <- sampling(cost_cp, data = cost_fake, chains = 4, iter = 1000)

cost_sim_df <- as.data.frame(cost_sim)
precis(cost_sim_df %>%  select(c(alpha, beta_vs, beta_hcc, beta_char)))

mcmc_dens(cost_sim, pars = c("alpha", "beta_vs", "beta_hcc", "beta_char"))


############
############ Now add age
###########
###



cost_df <- df %>% mutate(cost = rnorm(nrow(df), 15 + visits + 2*hcc_conds +
                                       3*char_score + 2*(age - 70)/8 , 10))
glimpse(cost_df)

cost_df_spl <- slice_sample(cost_df, n = 200, replace = FALSE)

### Use recipes to normalize age
 library(tidymodels)
rec <- recipe(cost ~  visits + hcc_conds + char_score + age, data = cost_df_spl) %>%
    step_normalize(age)


### weird to bake  a subset but whatever
cost_prep <- rec %>% prep(cost_df) %>% bake(cost_df_spl)

cost_fake = list(cost = cost_prep$cost, N = nrow(cost_prep), visits = cost_prep$visits,
                 hcc_count = cost_prep$hcc_conds, char_score = cost_prep$char_score,
                 age_norm = cost_prep$age)

## Use same stan file
cost_cp <- stan_model("cost_01.stan")
library(rethinking)
cost_sim <- sampling(cost_cp, data = cost_fake, chains = 4, iter = 1000)

cost_sim_df <- as.data.frame(cost_sim)
precis(cost_sim_df %>%  select(c(alpha, beta_vs, beta_hcc, beta_char, beta_age)))

mcmc_dens(cost_sim, pars = c("alpha", "beta_vs", "beta_hcc", "beta_char", "beta_age"))

### See how well params are recovered- fake data generated with "true" params vs model fit to fake data
posterior_params <- as.data.frame(cost_sim, pars = c("alpha", "beta_vs", "beta_hcc", "beta_char", "beta_age"))
true_params <- c(15, 1, 2, 3, 2)
mcmc_recover_hist(posterior_params, true_params)

y_rep <- as.matrix(cost_sim, pars = "y_rep")
ppc_dens_overlay(y = cost_df_spl$cost, yrep = y_rep[1:100, ])

### Do stratified, etc for the above
###
###
###
###  FAKE, fake dpg data- generating cost using just the priors and generating function

priors_cp <- stan_model("cost_01_priors.stan")
priors_sim <- sampling(priors_cp, chains = 1, iter = 1, algorithm = 'Fixed_param')

samps <- rstan::extract(priors_sim)

### fit the model to this fake generated data as in Gabry Helsinki

cost_fake_fake <- list(N = 200, visits = samps$visits[1, ],
                            hcc_count = samps$hcc_count[1, ], char_score = samps$char_score[1, ],
                            age_norm = samps$age[1, ], cost = samps$cost[1, ])

samps$cost[1,] %>% class()

cost_sim_fake <- sampling(cost_cp, data = cost_fake_fake, chains = 4, iter = 1000)

cost_sim_fake_df <- as.data.frame(cost_sim_fake)
precis(cost_sim_fake_df %>%  select(c(alpha, beta_vs, beta_hcc, beta_char, beta_age)))

mcmc_dens(cost_sim_fake, pars = c("alpha", "beta_vs", "beta_hcc", "beta_char", "beta_age"))

### See how well params are recovered- fake data generated with "true" params vs model fit to fake data
posterior_params <- as.data.frame(cost_sim_fake, pars = c("alpha", "beta_vs", "beta_hcc", "beta_char", "beta_age"))
true_params <- c(samps$alpha, samps$beta_vs, samps$beta_hcc, samps$beta_char, samps$beta_age)
mcmc_recover_hist(posterior_params, true_params)

## Check density overlay- fake data vs model fit to fake data
y_rep <- as.matrix(cost_sim_fake, pars = "y_rep")
ppc_dens_overlay(y = cost_fake_fake$cost, yrep = y_rep[1:200, ])

## rootogram() only for count data
ppc_rootogram()






