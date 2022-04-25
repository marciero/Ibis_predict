library(tidymodels)
library(embed)


## Using embed() vs rstanarm as in Ch_18.R.
## Change success to factor rather than TRUE/FALSE


climbers2 <- climbers %>%
    mutate(success_fact = as.factor(success))

climbers2 <- climbers %>%
    mutate(success = case_when(success == TRUE ~ "successful",
                               success == FALSE ~ "fail"))


climb2 <-
    recipe(success ~ age + oxygen_used + expedition_id, data = climbers2) %>%
    step_lencode_bayes(
        expedition_id,
        outcome = vars(success),
        options = opts  ## defined in Ch_18.R
    ) %>%
    prep(training = climbers2)

bayes_param <- tidy(climb2, number = 1)

## Now execute/define climb_glmer on line 65 Ch_18.R using this data
## Still does not agree but at least getting what appears to be log odds.
## Try the grant data both ways.

grants_glmer <-
    recipe(class ~ ., data = grants_other) %>%
    step_lencode_bayes(
        sponsor_code,
        outcome = vars(class),
        options = opts
    ) %>%
    prep(training = grants_other)

tidy(grants_glmer, number = 1)

grant_model <- stan_glmer(
    class ~ . + (1|sponsor_code),
    data = grants_other, family = binomial,
    prior_intercept = normal(0, 2.5, autoscale = TRUE),
    prior = normal(0, 2.5, autoscale = TRUE),
    prior_covariance = decov(reg = 1, conc = 1, shape = 1, scale = 1),
    chains = 4, iter = 1000, seed = 84735
)

tidy(grant_model, effects = "ran_vals")

######
######  Try non-bayes, lencode_glm

climb2_glm <-  recipe(success ~ age + oxygen_used + expedition_id, data = climbers2) %>%
    step_lencode_glm(
        expedition_id,
        outcome = vars(success),
    ) %>%
    prep(training = climbers2)

glm_params <- tidy(climb2_glm, number = 1)


log_odds <- climbers %>% group_by(expedition_id) %>% summarize(prop = mean(success)) %>%
    mutate(log_odds = log(prop/(1-prop)))

###   note that the last line below requires the model that follows it- line 87

log_odds %>% inner_join(glm_params, by = c( "expedition_id" = "level" )) %>%
    transmute(level = expedition_id, log_odds, lencode_glm = value ) %>%
    inner_join(bayes_param, by = "level") %>%
    transmute(level, log_odds, lencode_glm, lencode_bayes = value) %>%
    inner_join(rstan_params, by = "level") %>%
    transmute(level, log_odds, lencode_glm, lencode_bayes, rstan_est = estimate) %>%
    mutate(rstan_adj = rstan_est - 1.78)

##
### So these still do not agree. BUT...
### edit: rstanarm is estimating the beta coeffs in the (generalized) linear model
###  while these embed recipe steps are estimating the log odds

climb_model <- stan_glmer(
    success ~ age + oxygen_used + (1 | expedition_id),
    data = climbers, family = binomial,
    prior_intercept = normal(0, 2.5, autoscale = TRUE),
    prior = normal(0, 2.5, autoscale = TRUE),
    prior_covariance = decov(reg = 1, conc = 1, shape = 1, scale = 1),
    chains = 4, iter = 1000, seed = 84735
)

rstan_params <- tidy(climb_model, effects = "ran_vals")

climb_model %>% as.data.frame() %>% summarize(median(`b[(Intercept) expedition_id:AMAD03107]`))
climb_model %>% as.data.frame() %>% summarize(median(`(Intercept)`))

### NMedian for this expedition should be approx sum of the above.


select(c("level", "value", "log_odds")) %>% mutate(lencode_glm = value) %>%
    select(-value) %>% inner_join()

bayes_param %>% transmute(level, lencode_bayes = value)
tidy(climb2, number = 1)

tidy(climb_model, effects = "ran_vals")


tidy(climb_model, effects = "ran_vals")

########
########  Below uses step_lencode_mixed. Is a no-intercept model See documentation
climb3 <-
    recipe(success ~ age + oxygen_used + expedition_id, data = climbers2) %>%
    step_lencode_mixed(
        expedition_id,
        outcome = vars(success),
        ##  options = opts  ## defined in Ch_18.R
    ) %>%
    prep(training = climbers2)

tidy(climb3, number = 1)

###

