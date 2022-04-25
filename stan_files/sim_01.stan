// model mod_01 from sim_01.R. For troubleshooting and to try  cmndstan
// Only generated data using mean and sd of actual hcc scores
data {
    int<lower=1> N;
}
model {
}
generated quantities {
    real age[N];
    real age_adj[N];
    int hcc_conds[N];
    real alpha = normal_rng(log(4), 0.5);
    real beta = normal_rng(log(1.75)/10, 0.2);
    real inv_phi = fabs(normal_rng(0, 1));

    for (n in 1:N) {
    age[n] = normal_rng(70, 8);
    age_adj[n] = (age[n] - 70)/8;
    hcc_conds[n] =
              neg_binomial_2_log_rng(alpha + beta * age_adj[n], inv(inv_phi));
    }
}
