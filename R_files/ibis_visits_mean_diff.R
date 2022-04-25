

## We can compare or contrast the difference in the number of visits, for `ibis` versus non users.
## One way to do this is to simulate data from populations differing only in ibis use, or simply change the
## `ibis` variable in the patient data we have to indicate use; that is, as if all existing patients had
##  used Ibis, then use the model to predict visits. Now do the same with `ibis` changed to
##  indicate non use. We can then get a distribution for the expected number of visits in each case,
##  and hence a distribution for the difference in the expected number of visits.
##  There are problems with doing it this way, however because this is the data we used to fit the model
##  in the first place, and also because ibis use was simulated based on this data.
##
##  Also check and compare lm() vs Bayes on cost model- should be consistent.


ibis_no <- ibis_fake %>% mutate(ibis = 0)
ibis_yes <- ibis_fake %>% mutate(ibis = 1)
ibis_no_pred <- posterior_predict(p.2.rstn, newdata = ibis_no) %>% as.data.frame()
ibis_yes_pred <- posterior_predict(p.2.rstn, newdata = ibis_yes) %>% as.data.frame()


ibis_means <- data.frame(yes = ibis_yes_pred %>%
                             summarize(across(.cols = everything(), mean)) %>% t(),
                         no = ibis_no_pred %>%
                             summarize(across(.cols = everything(), mean)) %>% t()) %>%
    mutate(diff_means = no - yes)



ibis_means %>%  pivot_longer(1:2, names_to = "ibis_use", values_to = "mean_visits") %>%
    ggplot(aes(mean_visits, fill = ibis_use))  +
    geom_histogram(alpha = 0.75, position = 'identity') +
    lims(x = c(0, 10))


ibis_means %>% ggplot(aes(x = diff_means))  +
    geom_histogram(aes(y = ..density..), colour="black", fill="white")  +
    geom_density(alpha=.2) +    ###  fill="#FF6666"
    lims(x = c(0, 1.5))

mean(ibis_means$diff_means)
### Why only 0.52 mean difference? How does this square with our distribution for beta_ibis? It has
### mean
precis(ibis_means)

#######
#######

## Bayes vs lm()
###  from hcc_comorb.Rmd
###  cost_df <- df %>% mutate(cost = rnorm(nrow(df), 15 + visits + 2*hcc_conds +
###  3*char_score + 2*(age - 71.6)/7.3 , 10))
###

lm_fit <- lm(cost ~ ., cost_prep)
tidy(lm_fit)

## Compare with bayes model-they are worse! (This seems to be fixed below)
precis(cost_sim_df %>% select(c(alpha, beta_vs, beta_hcc, beta_char, beta_age)))

## Try from scratch and  compare using rstanarm. Use df from morbidity.R line18
##


### take sample 200 for ease computation

df_spl <- slice_sample(df, n = 200, replace = FALSE)

rec <- recipe(df_spl) %>%
    step_normalize(age)

cost_prep <- rec %>% prep(df_spl) %>% bake(df_spl)

cost_data <- cost_prep %>% mutate(cost = rnorm(nrow(cost_prep), 15 + visits + 2*hcc_conds +
                                          3*char_score + 2*age, 10))

cost_data %>% ggplot(aes(age, cost)) + geom_point() + geom_smooth(method = "lm", se = F)
cost_data %>% ggplot(aes(hcc_conds, cost)) + geom_point() + geom_smooth(method = "lm", se = F)

lm(cost ~ visits + hcc_conds + char_score + age, cost_data) %>% summary()

stn_mod <- stan_glm(cost ~ visits + hcc_conds + char_score + age, cost_data,
         family = gaussian,
         prior_intercept = normal(0, 1, autoscale = T),
         prior = normal(1, 1, autoscale = T),
         prior_aux = exponential(1, autoscale = TRUE),
         chains = 4, iter = 2000)


prior_summary(stn_mod, digits = 2)
### the means agree  with lm() model
stn_mod %>% tidy()
as.data.frame(stn_mod) %>% precis()

##Now try stan model
cost_fake = list(cost = cost_data$cost, N = nrow(cost_data), visits = cost_data$visits,
                 hcc_count = cost_data$hcc_conds, char_score = cost_data$char_score,
                 age_norm = cost_data$age)

## Use same stan file
cost_cp <- stan_model("cost_01.stan")

cost_sim <- sampling(cost_cp, data = cost_fake, chains = 4, iter = 1000)

cost_sim_df <- as.data.frame(cost_sim)
precis(cost_sim_df %>%  select(c(alpha, beta_vs, beta_hcc, beta_char, beta_age, sigma)))



