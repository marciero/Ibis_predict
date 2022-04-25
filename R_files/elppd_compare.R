## Example computing elppd  1) in generated quantities, 2) offline/manually by computing the likelihoods,
## 3) inspecting lp__ comparing.  from Xulong Wang example (see pdf).
##

stan_code <- "
data {
  int<lower=0> N;
  real y[N];
}
parameters {
  real mu;
  real<lower=machine_precision()> sigma;
}
model {
  y ~ normal(mu, sigma);
}
generated quantities {
  real lpd;
  lpd = normal_lpdf(y | mu, sigma);
}
"
dat <- list(N = 20, y = rnorm(20, 1, 1)) # pseudo-data
model <- stan(model_code = stan_code, data =  dat) # compilation myfit <- sampling(model, data = dat)

model %>% as.data.frame() %>% View()

 y<- extract(model)
## from Xulong Wang example (see pdf).
my_lpd = map(1:4000, function(i)
    sum(dnorm(dat$y, mean = y$mu[i], sd = y$sigma[i], log = TRUE)))
## We reproduce in a data frame. ## Below, lpd and lpd2 are the same. lp3 is the coefficient adjustment but is not exact.
## See also Sec 30.1 stan users manual. Constants are not important for model comparisons with
## same number of  observations
modeldf <- model %>% as.data.frame()
modeldf <- modeldf %>% group_by(1:n()) %>%
    mutate(lp2 =  sum(dnorm(dat$y, mean = mu, sd = sigma, log = TRUE))) %>% ungroup() %>%
    mutate(lp3 = lp2 -20*log(1/sqrt(2*pi)))





