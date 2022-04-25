// Model hcc conditions count nb glm with intercept and age predictor.
// For use with fake data added to augusta data
data {
    int<lower=1> N;
    vector[N] age;
    int hcc_conds[N];
}
transformed data {
    vector[N] agenorm = (age - mean(age))/sd(age);
}
parameters {
    real alpha;
    real beta;
    real<lower = 0> recip_phi;
}
transformed parameters {
    real phi = inv(recip_phi);
}
model {
   alpha ~ normal(log(3.5), 1);
   beta ~ normal(0.5, 1);  // for every sd in age mult hcc by exp(0.5)
   // beta ~ normal(log(1.75)/10, 0.2);  // from intercept only dgp
   recip_phi ~ normal(0,1);
   hcc_conds ~ neg_binomial_2_log(alpha + beta * agenorm, phi);
}
generated quantities {
    int y_rep[N];
    for (n in 1:N) {
    y_rep[n] = neg_binomial_2_log_rng(alpha + beta * agenorm[n], phi);
    }
}
