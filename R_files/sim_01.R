### Use data generated in hcc_augusta.R
### generate fake data
###
## create columns for each HCC
##
augusta_unique_hcc %>% unnest(c(data))
augusta_wide <- augusta_unique_hcc %>%  unnest(HCC_unique) %>% mutate(HCC = 1) %>%
    pivot_wider(names_from = HCC_unique, names_glue = "{.value}_{HCC_unique}",
                values_from = HCC, values_fill = 0)

is.factor(augusta_wide$HCC_111)

library(caret)
colnames(augusta_wide[, nearZeroVar(augusta_wide %>% select(-data))])
nzv<- nearZeroVar(augusta_wide %>% select(-data), saveMetrics = TRUE)
colnames(augusta_wide[, nearZeroVar(augusta_wide %>% select(-data), freqCut = 99/1)])

table(augusta_wide$HCC_83)

augusta_unique_hcc %>% unnest(HCC_unique) %>%  group_by(HCC_unique) %>%
    summarize(count = n()) %>% arrange(desc(count))

##  See loose notes. sim_01.stan is No model, only generated data.

mod_01 <- stan_model("sim_01.stan")
mod_01_sim <- sampling(mod_01, data = list(N = 100), chains = 1, iter = 1, algorithm = 'Fixed_param')

tidy(mod_01_sim)   %>% filter(term == "inv_phi") ### only generates a single value of alpha with iter = 1
 ## mod_01_sim %>% as.data.frame() %>% View()
posterior <- extract(mod_01_sim)
posterior$hcc_conds[1, ]  %>% as.data.frame() %>%  ggplot(aes(.)) + geom_bar()

#########3
#########  Now feed fake data into model


## To use ppc_ we need to generate N y_reps in the model. Above we had 100 for simplicity.
y_rep <- as.matrix(hcc_pois_sim, pars = "y_rep")
ppc_dens_overlay(y = augusta_unique_hcc$hcc_conds, yrep = y_rep[1:200, ])

prop <- function(x) mean(x == 2)
ppc_stat(y = augusta_unique_hcc$hcc_conds, yrep = y_rep, stat = "prop")




####  Same thing with cmndstan
library(cmdstanr)
## See stan_exercises, or vignette
## again only generated data
mod_01_sim_cmd_stan <- cmdstan_model("sim_01.stan")

fit <- mod_01_sim_cmd_stan$sample(data = list(N=nrow(augusta_unique_hcc)), seed = 123,
                      chains = 1,
                      parallel_chains = 1,
                      fixed_param = TRUE,    ## compare "algorithm = 'Fixed_param'" in rstan version
                     iter_sampling  = 1000,
                      refresh = 1000)   ## is just the frequency progress is displayed on screen

fit$summary()
draws_arr <- fit$draws()
str(draws_arr)
draws <- fit$draws(format = "df")
str(draws)
draws %>% select(starts_with("hcc_conds")) %>% mcmc_hist()  ## need df format to use select


y_rep <- as.matrix(draws %>% select(starts_with("hcc")))
ppc_dens_overlay(y = augusta_unique_hcc$hcc_conds, yrep = y_rep)



########
########  Now try hcc_conds ~ NB(alpha +beta*age). sim_02.stan
########
ggplot() + geom_function(fun = ~.x^2)

ggplot() + geom_function(fun = ~dpois(.x, 3))

mod_02 <- stan_model("sim_02.stan")

idx <- sample(1:nrow(augusta_unique_hcc), 100)
sbst <- augusta_unique_hcc[idx, ]

mod_02_sim <- sampling(mod_02, data = list(N = length(idx), age = sbst$age,
                                           hcc_conds = sbst$hcc_conds),
                       chains = 4, iter = 1000)

##############


print(mod_02_sim)
mcmc_intervals(mod_02_sim, pars = c("alpha", "recip_phi"))
#mcmc_parcoord(posterior)

tidy(mod_02_sim, conf.int = TRUE, conf.level = 0.95)
tidy(mod_02_sim)
mcmc_trace(mod_02_sim, pars = c("alpha", "inv_phi"))
mcmc_hist(mod_02_sim, pars = c("alpha", "inv_phi"))
summary(mod_02_sim)
mcmc_dens_overlay(mod_02_sim, pars = c("alpha", "inv_phi"))
mcmc_acf(mod_02_sim, pars = c("alpha", "inv_phi"))
rhat(mod_02_sim, pars = c("alpha", "inv_phi"))
neff_ratio(mod_02_sim, pars = c("alpha", "inv_phi"))

mcmc_dens(mod_02_sim, pars = c("alpha", "inv_phi"))  ## Is all the chains combined
mcmc_areas(mod_02_sim, pars = c("alpha", "inv_phi"), prob = 0.95)

## To use ppc_ we need to generate N y_reps in the model. Above we had 100 for simplicity.
y_rep <- as.matrix(mod_02_sim, pars = "y_rep")
ppc_dens_overlay(y = sbst$hcc_conds, yrep = y_rep[1:200, ]) + xlim(0,20)

prop <- function(x) mean(x == 5)
ppc_stat(y = sbst$hcc_conds, yrep = y_rep, stat = "prop")

#### do ppc intervals- see hcc_augusta nb intercept only model

mean_y_rep <- colMeans(y_rep)
std_resid <- (augusta_unique_hcc$hcc_conds - mean_y_rep) / sqrt(mean_y_rep)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)

ppc_rootogram(sbst$hcc_conds, yrep = y_rep)
ppc_bars(sbst$hcc_conds, yrep = y_rep)

###
###  Full data set
  mod_02_sim_full <- sampling(mod_02, data = list(N = nrow(augusta_unique_hcc), age = augusta_unique_hcc$age,
  hcc_conds = augusta_unique_hcc$hcc_conds),
  chains = 4, iter = 1000)

  post_mod2 <- extract(mod_02_sim_full)
 tidy(mod_02_sim_full) %>% head()
 y_rep <- as.matrix(mod_02_sim_full, pars = "y_rep")
 ppc_dens_overlay(y = augusta_unique_hcc$hcc_conds, yrep = y_rep[1:200, ]) + xlim(0,20)


 prop <- function(x) mean(x == 3)
 ppc_stat(y = augusta_unique_hcc$hcc_conds, yrep = y_rep, stat = "prop")




######## T
######## troubleshoot with sim_02_trouble.stan. Try pest example as template
########
mod_tr <- stan_model("sim_02_trouble.stan")

fitted_model_dgp_NB <-
    sampling(
        comp_dgp_multiple_NB,
        data = list(N =4000),                      ## N = nrow(pest_data)),
        chains = 1,
        cores = 1,
        iter = 1,
        algorithm = 'Fixed_param',
        seed = 123
    )
samps_dgp_NB <- rstan::extract(fitted_model_dgp_NB)

stan_dat_fake_NB <- list(
    N = 4000,
    traps = augusta_unique_hcc$age[1:4000],
    complaints = augusta_unique_hcc[1:4000]
)


mod_tr_sim <- sampling(mod_tr, data = stan_dat_fake_NB,
                       chains = 4, iter = 200)

mod_tr_sim <- sampling(mod_tr, data = list(N = nrow(augusta_unique_hcc), age = augusta_unique_hcc$age,
                                           hcc_conds = augusta_unique_hcc$hcc_conds),
                       chains = 4, iter = 200)

mod_tr <- stan_model("sim_02_trouble.stan")

tidy(mod_tr_sim)
mcmc_trace(mod_tr_sim)











## try rstanarm
hcc_rstn <- stan_glm(
    hcc_conds ~ age,
    data = augusta_unique_hcc, family = neg_binomial_2,
    prior_intercept = normal(0, 2.5, autoscale = TRUE),
    prior = normal(1, 1, autoscale = TRUE),
    prior_aux = normal(1,1, autoscale = TRUE),  ##
    chains = 4, iter = 2000, seed = 84735)


########
########

pp_check(hcc_rstn) +
    xlab("hcc_count")   ### same as mcmc_dens_overlay
tidy(hcc_rstn)

mcmc_trace(hcc_rstn)
mcmc_dens_overlay(hcc_rstn)
mcmc_acf(hcc_rstn)

pp_check(hcc_rstn, plotfun = "hist", nreps = 5)+
    xlab("hcc_count")

hcc_rstn_df <- as.data.frame(hcc_rstn_sim)

predictions_no <- posterior_predict(hcc_rstn)                 ##newdata = artist_means...

# Plot the posterior predictive intervals for some
idx <- sample(1:nrow(augusta_unique_hcc), 50)
ppc_intervals(augusta_unique_hcc$hcc_conds[idx], yrep = predictions_no[, idx],
              prob_outer = 0.80) +
    ggplot2::scale_x_continuous()











