### Augusta data clean and match hcc from CMS icd mapping table.
### Then determine distributions of hccs, model hcc as poisson, negative binomial.

library(tidyverse)
library(rstan)
library(rstanarm)
library(bayesplot)
library(bayesrules)
library(broom)
library(broom.mixed)
## library(rethinking)
library(janitor)

augusta_raw <- read_csv("data/augusta_2020.csv")

## get rid of spaces and also remove decimals in ICD codes, as the lookup table has none.
augusta <- augusta_raw %>% rename_with(~ str_replace_all(., " ", "_"), "PATIENT NBR": "CPT 15") %>%
    mutate(across(DIAG_1:DIAG_15, ~str_remove(., "\\.")))

nrow(augusta)
n_distinct(augusta$PATIENT_NBR)

augusta %>% group_by(PATIENT_NBR) %>% summarize(count = n())  %>%
    ggplot(aes(count)) +
    geom_histogram(aes(y = ..density..), binwidth = 1, color = "white")

### create column with visit number
visit_df <- augusta %>% group_by(PATIENT_NBR) %>% summarize(count = n()) %>%
    mutate(visit = map(count, function(i){c(1:i)})) %>%
    unnest(visit)

augusta_vis <- cbind(visit_df[, 3], augusta)

## The HCC lookup table
hcc <- read_csv("data/HCC_mappings_2020.csv")
colnames(hcc)[3] <- "HCC"
hcc <- hcc %>% filter(!is.na(HCC))
hcc$HCC <- as.factor(hcc$HCC)
hcc$Code <- as.factor(hcc$Code)

augusta_long <- augusta_vis %>% pivot_longer(DIAG_1:DIAG_15, names_to = "diagnosis_no", values_to = "code") %>%
         select(PATIENT_NBR, visit, FIN_CLASS, DRG, diagnosis_no, code) %>% select(-DRG) %>%
    filter(!is.na(code))
augusta_long$code <- as.factor(augusta_long$code)

## some codes mapped to more that one HCC.
augusta_long_hcc <- augusta_long %>% left_join(hcc, by = c("code" = "Code")) # %>%
    ##  filter(!is.na(HCC))  ## previously we dropped records with no HCC.

 n_distinct(augusta_long_hcc$HCC)

 n_distinct(augusta_long_hcc$HCC)

 ### Pick out the distinct HCC for each patient, number of distinct HCCs, number of visits, generate age based on
 ### number of hcc categories, generate financial/income, gender.

 N= nrow(augusta_long_hcc %>% group_nest(PATIENT_NBR))
 set.seed(123)
 augusta_unique_hcc <- augusta_long_hcc %>% group_nest(PATIENT_NBR) %>%
     mutate(HCC_unique = map(data, ~unique(.x$HCC))) %>%
     mutate(visits = map_int(data, ~max(.x$visit))) %>%
   mutate(hcc_conds = map_int(HCC_unique, ~length(.) - sum(is.na(.x)))) %>%  ## subtracts NA, which were counted distinct HCC
     mutate(age = 65 + 2*hcc_conds + rnorm(N, 0 , sd = 5))  %>%   ## fake data
     mutate(across(age, ~floor(.x))) %>%
     mutate(gender = as.factor(sample(c("M", "F", "NB"), N, prob = c(0.49, 0.49, 0.02),  replace = TRUE)))%>%  ## more fake data
     mutate(financial = ordered(sample(c(1:3), N, replace = TRUE)))

 is.factor(augusta_unique_hcc$financial)


  ## distinct patients
 n_distinct(augusta_unique_hcc$PATIENT_NBR)

 augusta_unique_hcc %>%    ggplot() + geom_bar(aes(age))

 ## dist of hcc frequencies
 augusta_unique_hcc %>% unnest(HCC_unique) %>% na.omit() %>%
     ggplot(aes(HCC_unique)) + geom_bar()

 augusta_unique_hcc %>% unnest(HCC_unique)

 ## dist of  patient hcc counts
   library(janitor)
   tabyl(augusta_unique_hcc$hcc_conds == 0)
   tabyl(augusta_unique_hcc$hcc_conds)

   augusta_unique_hcc %>%
       ggplot() + geom_bar(aes(hcc_conds))

   # visits
   augusta_unique_hcc %>%
       ggplot() + geom_bar(aes(visits)) +
       xlim(0, 20)

  augusta_unique_hcc %>% summarize(mean = mean(hcc_conds), sd = sd(hcc_conds))

## model hcc counts by patient as neg binomial
  data.frame(x = c(0:20)) %>% mutate(dens = dnbinom(x, mu = , size = 6.3)) %>%
      ggplot(aes(x, dens)) + geom_bar(stat = "identity")

hcc_pois_model <- "
   data {
   int N;
   int hcc_count[N];
   }
   parameters {
   real<lower = 0> alpha;
   }
   model {
    alpha ~ normal(log(4.6), 1);
    hcc_count ~ poisson_log(alpha);
   }
      generated quantities {
       int y_rep[N];
       for (i in 1:N)
    y_rep[i] = poisson_log_rng(alpha);
    }
"
hcc_pois_sim <- stan(model_code = hcc_pois_model, data = list(N = nrow(augusta_unique_hcc),
                                                     hcc_count = augusta_unique_hcc$hcc_conds),
               chains = 4, iter = 1000
                )
hcc_pois_sim %>% as.data.frame() %>% View()
mcmc_hist(hcc_pois_sim, pars = "y_rep[1]")
posterior <- extract(hcc_pois_sim)
mcmc_hist(hcc_pois_sim, pars = "alpha")


## To use ppc_ we need to generate N y_reps in the model. Above we had 100 for simplicity.
y_rep <- as.matrix(hcc_pois_sim, pars = "y_rep")
ppc_dens_overlay(y = augusta_unique_hcc$hcc_conds, yrep = y_rep[1:200, ])

prop <- function(x) mean(x == 2)
ppc_stat(y = augusta_unique_hcc$hcc_conds, yrep = y_rep, stat = "prop")

### See {bayeplot}, and Gabry, Bettancourt vignettes for more ppc_  post predictive

posterior$lp
## the posterior alphas are just the ave hcc with ;very small sd... makes sense.
## log(3)
   ###  hcc_pois_sim, pars = "lp__", include=FALSE) # include passed to rstan::extract
mcmc_intervals(hcc_pois_sim, pars = "lp__")
print(hcc_pois_sim)


#mcmc_parcoord(posterior)
#library(broom.mixed)
tidy(hcc_pois_sim, conf.int = TRUE, conf.level = 0.95)
tidy(hcc_pois_sim)
mcmc_trace(hcc_pois_sim, pars = "alpha")
mcmc_hist(hcc_pois_sim, pars = "alpha")

mcmc_trace(hcc_pois_sim, pars = "alpha", size = 0.5) +
    xlab("iteration")
mcmc_dens_overlay(hcc_pois_sim, pars = "alpha")
mcmc_acf(hcc_pois_sim, pars = "alpha")
rhat(hcc_pois_sim, pars = "alpha")
neff_ratio(hcc_pois_sim, pars = "alpha")

mcmc_dens(hcc_pois_sim, pars = "alpha")  ## Is all the chains combined
mcmc_areas(hcc_pois_sim, pars = "alpha", prob = 0.95)

## Try gamma poisson, which we have analytic form for
plot_gamma(2,.5)
sum(augusta_unique_hcc$hcc_conds)
plot_gamma_poisson(2, .5, sum_y = sum(augusta_unique_hcc$hcc_conds), n = nrow(augusta_unique_hcc))

## Note that sd of y_reps is less than the actual sd. This is because Poisson has mean = variance and cannot
## adjust separately.
sum(augusta_unique_hcc$hcc_conds)/nrow(augusta_unique_hcc)
sd(augusta_unique_hcc$hcc_conds)

#####

mean_y_rep <- colMeans(y_rep)
std_resid <- (augusta_unique_hcc$hcc_conds - mean_y_rep) / sqrt(mean_y_rep)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)

ppc_rootogram(augusta_unique_hcc$hcc_conds, yrep = y_rep)
ppc_bars(augusta_unique_hcc$hcc_conds, yrep = y_rep)

## also have ppc_intervals but would want x = predictor in linear model. See Gabry Helsinki Pest example
ppc_intervals(y = augusta_unique_hcc$hcc_conds,yrep = y_rep) # wont make sense

### Now try with nb model. Note three parametrizations for nb model in stan. We dont need
### to use log forms; eg neg_binomial_2_log,  because it is only alpha that would get exp. If we
### had a full linear model with predictors we would.
###

hcc_nb_model <- "
   data {
   int N;
   int hcc_count[N];
   }
   parameters {
   real<lower = 0> alpha;
   real<lower = 0> inv_phi; // variance is lambda + lambda^2/phi
   }
   transformed parameters {
     real phi = inv(inv_phi);
   }
   model {
    alpha ~ normal(log(3.5), 1);
    inv_phi ~ normal(0,1);
    hcc_count ~ neg_binomial_2_log(alpha, phi);
   }
   generated quantities {
       int y_rep[N];
       for (i in 1:N)
    y_rep[i] = neg_binomial_2_log_rng(alpha, phi);
    }
"
hcc_nb_sim <- stan(model_code = hcc_nb_model, data = list(N = nrow(augusta_unique_hcc),
                                                         hcc_count = augusta_unique_hcc$hcc_conds),
                chains = 4, iter = 2000
)


print(hcc_nb_sim)
mcmc_intervals(hcc_nb_sim, pars = c("alpha", "inv_phi"))
#mcmc_parcoord(posterior)
ext <- rstan::extract(hcc_nb_sim)
tidy(hcc_nb_sim, conf.int = TRUE, conf.level = 0.95)
tidy(hcc_nb_sim) %>% head()


mcmc_trace(hcc_nb_sim, pars = c("alpha", "inv_phi"))
mcmc_hist(hcc_nb_sim, pars = c("alpha", "inv_phi"))
summary(hcc_nb_sim)
mcmc_dens_overlay(hcc_nb_sim, pars = c("alpha", "inv_phi"))
mcmc_acf(hcc_nb_sim, pars = c("alpha", "inv_phi"))
rhat(hcc_nb_sim, pars = c("alpha", "inv_phi"))
neff_ratio(hcc_nb_sim, pars = c("alpha", "inv_phi"))

mcmc_dens(hcc_nb_sim, pars = c("alpha", "inv_phi"))  ## Is all the chains combined
mcmc_areas(hcc_nb_sim, pars = c("alpha", "inv_phi"), prob = 0.95)

## To use ppc_ we need to generate N y_reps in the model. Above we had 100 for simplicity.
y_rep <- as.matrix(hcc_nb_sim, pars = "y_rep")
ppc_dens_overlay(y = augusta_unique_hcc$hcc_conds, yrep = y_rep[1:200, ])

prop <- function(x) mean(x == 0)
ppc_stat(y = augusta_unique_hcc$hcc_conds, yrep = y_rep, stat = "prop") +  labs(x = "prop 0")


### random sample of data, check against prediction interval
 ext <- rstan::extract(hcc_nb_sim)
idx <- sample(1:nrow(augusta_unique_hcc), 25)
predictions_pois <- ext$y_rep[, idx]
y_rep_subs <- as.matrix(hcc_nb_sim, pars = "y_rep")
ppc_intervals(augusta_unique_hcc$hcc_conds[idx], yrep = predictions_pois[, idx],
              prob_outer = 0.90) +
    ggplot2::scale_x_continuous()






mean_y_rep <- colMeans(y_rep)
std_resid <- (augusta_unique_hcc$hcc_conds - mean_y_rep) / sqrt(mean_y_rep)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)

ppc_rootogram(augusta_unique_hcc$hcc_conds, yrep = y_rep)
ppc_bars(augusta_unique_hcc$hcc_conds, yrep = y_rep)

#####
##### plot age vs hcc

augusta_unique_hcc %>% ggplot(aes(age, hcc_conds, color = financial)) + geom_point(position = "jitter")

########
######## Try using rstanarm and ulam
######## Using BayesRules as guide.

hcc_rstn_sim <- stan_glm(
    hcc_conds ~ 1,
    data = augusta_unique_hcc, family = neg_binomial_2,
    prior_intercept = normal(0, 2.5, autoscale = TRUE),
    # prior = normal(0, 2.5, autoscale = TRUE),
    prior_aux = exponential(1, autoscale = TRUE),  ##
    chains = 4, iter = 2000, seed = 84735)

## Note: hcc_rstn_pois used below is the above with poisson model
# Check out the priors
prior_summary(hcc_rstn_sim)
pp_check(hcc_rstn_sim) +
    xlab("hcc_count")   ### same as mcmc_dens_overlay
tidy(hcc_rstn_sim)

mcmc_trace(hcc_rstn_sim)
mcmc_dens_overlay(hcc_rstn_sim)
mcmc_acf(hcc_rstn_sim)

set.seed(1)
pp_check(hcc_rstn_sim, plotfun = "hist", nreps = 5)+
    xlab("hcc_count")

## From {bayesrules}
cv_procedure <- prediction_summary_cv(
    model = hcc_rstn_sim, data = augusta_unique_hcc, k = 3)

cv_procedure$folds
cv_procedure$cv

cv_procedure_pois <- prediction_summary_cv(
    model = hcc_rstn_pois, data = augusta_unique_hcc, k = 3)

cv_procedure_pois$cv

loo_nb <- loo(hcc_rstn_sim)  ## would fit 4k+ models, but only 20 seconds?
loo_pois <- loo(hcc_rstn_pois)
loo_compare(loo_pois, loo_nb)

##### quick look at visits, but wont do full modelling
augusta_unique_hcc %>% ggplot(aes(visits)) + geom_bar()

## now generate simulated patient data with these distributions of hcc.
