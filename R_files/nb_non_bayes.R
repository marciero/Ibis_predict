#########
#########
#########  Non bayes
## Note that tidymodels does not currently support nb regression

pois_mod <- glm(hcc_conds ~ age, data = augusta_unique_hcc,
              family = poisson(link = "log"))

summary(nb_mod)
library(MASS)
pois_mod <- glm(hcc_conds ~ age, data = augusta_unique_hcc,
                family = negative.binomial(theta = 5))

summary(pois_mod)

pois_mod <- glm.nb(hcc_conds ~ age, data = augusta_unique_hcc,
                )

summary(pois_mod)
