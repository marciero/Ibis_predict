### model visits as poisson or nb
### Is basis for ibis_01.Rmd
library(rethinking)
df <- augusta_unique_hcc
augusta_unique_hcc %>% ggplot(aes(visits)) + geom_bar()

## Do not seem correlated
augusta_unique_hcc %>% ggplot(aes(hcc_conds, visits)) + geom_point() +
    geom_smooth(method = lm, se = F)
augusta_unique_hcc %>% ggplot(aes(hcc_conds, visits)) + geom_point(position = "jitter")

mean(df$visits)
sd(df$visits)

cor(df$visits, df$hcc_conds)

## Model as lm. First non bayes
vm.0 <- lm(visits ~ hcc_conds, df)
 summary(lm.1)
 tidy(lm.1)

 predict(p.0.1, data.frame(hcc_conds = c(2,5,10)))
 # Bayes  Note that below is GAUSSIAN so will be bad
vm.1 <- stan_glm(visits ~ hcc_conds, df,
                 family = gaussian,
                  prior_intercept = normal(0, 1, autoscale = T),
                 prior = normal(1, 1, autoscale = T),
                 prior_aux = exponential(1, autoscale = TRUE),
                 chains = 4, iter = 4000
                 )

prior_summary(vm.1)

tidy(vm.1)

posterior_predict(vm.1, newdata = data.frame(hcc_conds = c(2,5,10)))

 mcmc_intervals(vm.1)
mcmc_hist(vm.1)

pp_check(vm.1) +
    xlab("visits") + xlim(-10, 20)

pp <- posterior_predict(vm.1, newdata = data.frame(hcc_conds = c(5,5))) %>%
    as.data.frame()

# 50 simulated model lines
df  %>% select(visits, hcc_conds) %>%
    add_epred_draws(vm.1, ndraws = 50) %>%
    ggplot(aes(x = hcc_conds, y = visits)) +
    geom_line(aes(y = .epred, group = .draw), alpha = 0.15) +
    geom_point(data = df, position = "jitter", size = 0.05)

######3
###### Now model as poisson or nb. Visits may be exponential dist
######  non bayes

p.0 <- glm(visits ~ 1, df, family = poisson)
summary(p.0)
tidy(p.0)

p.0.1 <- glm(visits ~ hcc_conds, df, family = poisson)
summary(p.0.1)
tidy(p.0.1)

##Bayes

p.1 <- stan_model("pois_01.stan")
p1sim <- sampling(p.1, data = list(N = nrow(df), visits = df$visits, hcc_conds = df$hcc_conds),
                  chains = 4, iter = 2000)
mcmc_hist(p1sim, pars = c("alpha", "beta"))
tidy(p1sim)

### Note precis on model vs the extracted- we get nice little histogram thumbnail in latter.
ext <- extract(p1sim)
precis(p1sim, pars = c("alpha", "beta"))
precis(ext, pars = c("alpha"))   ## ignores alpha
ext_reth <- extract.samples(p1sim)
precis(ext_reth$alpha)
precis(ext_reth, pars = c("alpha"))

## Or use rstanarm
vm.2 <- stan_glm(visits ~ hcc_conds, df,
                 family = poisson,
                 prior_intercept = normal(0, 1, autoscale = T),
                 prior = normal(1, 1, autoscale = T),
                 prior_aux = exponential(1, autoscale = TRUE),
                 chains = 4, iter = 2000
)

prior_summary(vm.2)

tidy(vm.2)

mcmc_intervals(vm.2)
mcmc_hist(vm.2)

pp_check(vm.2) +
    xlab("visits") + xlim(-5, 15)
##
## intercept only
vm.0 <- stan_glm(visits ~ 1, df,
                 family = poisson,
                 prior_intercept = normal(0, 1, autoscale = T),
                # prior = normal(1, 1, autoscale = T),
                 prior_aux = exponential(1, autoscale = TRUE),
                 chains = 4, iter = 4000
)

prior_summary(vm.0)

tidy(vm.0)

mcmc_intervals(vm.0)
mcmc_hist(vm.0)

pp_check(vm.0) +
    xlab("visits") + xlim(-5, 15)


#########
#########  Truncated poisson. Have to use stan. Not able to generate y_reps yet

tpois <- stan_model("trunc_pois.stan")
tpsim <- sampling(tpois, data = list(N = nrow(df), visits = df$visits, hcc_conds = df$hcc_conds),
                      chains = 4, iter = 2000)
 mcmc_hist(tpois_sim)
 tidy(tpois_sim)
 post <- extract(tpois_sim)  ## rstan::extract(tpois_sim)



### Fake ibis data-includes indicator for ibis use

ibis_fake <- df %>% mutate(wt = case_when(visits < 5 ~ 0.65,
                                          visits >= 5 ~ 0.45)) %>%
    group_by(1:n()) %>%
    mutate(ibis = sample(c(0, 1), size = 1, replace = TRUE, prob =  c(1-wt, wt))) %>%
    mutate(ibis_idx = ibis +1) %>% ungroup()
  ## %>%  mutate(across(ibis, ~as.factor(.)))

# library(janitor)
ibis_fake %>% tabyl(visits, ibis)

ibis_fake %>% ggplot(aes(age, visits)) + geom_point(aes(color = as.factor(ibis)))

ibis_fake %>% filter(visits <= 10) %>% group_by(visits, ibis) %>%
    summarize(count = n())

ibis_fake %>% filter(visits <= 10) %>%
    ggplot(aes(visits)) +
    geom_bar(aes(fill = as.factor(ibis)), position = "fill") +
    scale_x_discrete(limits = factor(c(1:10)))

ibis_fake %>% ggplot(aes(visits)) + geom_bar(aes(color = as.factor(ibis)))

## can see interaction, ibis:hcc_conds
ibis_fake %>% ggplot(aes(hcc_conds, visits)) + geom_point(aes(color = as.factor(ibis))) +
    geom_smooth(aes(color = as.factor(ibis)),method = lm, se = F)


#########3
#########  poisson model with  hcc, ibis, then add interaction of ibis and hcc counts
#########  To see if this makes sense, try  (of course it make sense becaues ibis modeled as hcc counts )

bin.1 <- glm(ibis ~ visits + hcc_conds, ibis_fake, family = binomial)
tidy(bin.1, exponentiate = T)

p.2.0 <- glm(visits ~ ibis + hcc_conds, ibis_fake,
             family = poisson)
tidy(p.2.0, exponentiate = F)

#####


p.2 <- stan_model("pois_02.stan")
p2sim <- sampling(p.2, data = list(N = nrow(df), visits = ibis_fake$visits, hcc_conds = ibis_fake$hcc_conds,
                                   ibis_idx = ibis_fake$ibis_idx),
                  chains = 4, iter = 2000)

ext2 <- extract.samples(p2sim)
precis(ext2, depth = 2, pars = c("alpha", "beta"))   ### Note these are on the log scale

##### troubleshoot above. modify model to reproduce p.2.0
p.20stan <- stan_model("pois_02_0.stan")
p20stnsim <- sampling(p.20stan, data = list(N = nrow(df), visits = ibis_fake$visits, hcc_conds = ibis_fake$hcc_conds,
                                   ibis = ibis_fake$ibis),
                  chains = 4, iter = 4000)

tidy(p20stnsim, exponentiate = T)
ext <- extract.samples(p20stnsim)
precis(ext, depth = 2)
## posterior_interval(as.matrix(p20stnsim))

########

diff_alpha <- exp(ext2$alpha[, 2]) - exp(ext2$alpha[, 1])
diff_beta <- exp(ext2$beta[, 2]) - exp(ext2$beta[, 1])
precis(list(diff_alpha = diff_alpha, diff_beta = diff_beta))

### check against rstanarm
### Finally, the simplest model does agree in both cases
#
p.2.rstn <- stan_glm(visits ~  ibis + hcc_conds, ibis_fake,
                               family = poisson,
                               prior_intercept = normal(0, 1, autoscale = T),
                                prior = normal(1, 1, autoscale = T),
                               prior_aux = exponential(1, autoscale = TRUE),
                               chains = 4, iter = 2000
)

p.2.rstn_df <- as.data.frame(p.2.rstn)
mean(p.2.rstn_df$`(Intercept)`)
p.2.rstn$coefficients
p.2.rstn$stan_summary

tidy(p.2.rstn)

posterior_interval(p.2.rstn)
posterior_predict(p.2.rstn, newdata = data.frame(ibis = 1, hcc_conds = 3)) %>%
    posterior_interval(prob = 0.95)

posterior_predict(p.2.rstn, newdata = data.frame(ibis = 0, hcc_conds = 5)) %>%
    mcmc_hist()


precis(p.2.rstn_df, depth = 2)

### overall mea
 ibis_fake %>% group_by(ibis) %>% summarize(ave_visits = mean(visits))

 ### odds ratio, Bayes factor
 pp <- posterior_predict(p.2.rstn, newdata = data.frame(ibis = c(0,1), hcc_conds = c(5,5))) %>%
     as.data.frame()
 colnames(pp) <- c("ibis_no", "ibis_yes")

 pp %>%  summarize(quantile_no = quantile(ibis_no, 0.75),
                      quantile_yes = quantile(ibis_yes, 0.75))

    pp %>% pivot_longer(c(1,2), names_to = "ibis", values_to = "visits") %>%  ggplot() +
        geom_histogram(aes(visits, fill = ibis), position = "dodge")

 ## likelihood ratio
mean(pp$ibis_yes > 3)/mean(pp$ibis_no > 3)

 ## Non-ibis users  with 10 hcc who do not use ibis are 47% more likely to have more than four visits than
 ## users,  while non users
 ## with 2 hcc are twice as likely to have more than four visits than users.

oneway.test(ibis_fake$visits ~ ibis_fake$ibis) %>% summary()
lm(visits ~ ibis, ibis_fake) %>% summary()

###
###
y_rep <- as.matrix(hcc_pois_sim, pars = "y_rep")
 samps <- extract.samples(hcc_pois_sim)
sampsdf <- as.data.frame(hcc_pois_sim)

sampsdf$`y_rep[1]`


