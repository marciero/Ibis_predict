library(tidyverse)
library(rstan)
library(rstanarm)
library(bayesplot)
library(bayesrules)
library(broom)
library(broom.mixed)
#library(rethinking)
library(janitor)

augusta_raw <- read_csv("data/augusta_2020.csv")

## get rid of spaces and also remove decimals in ICD codes, as the lookup table has none.
augusta <- augusta_raw %>% rename_with(~ str_replace_all(., " ", "_"), "PATIENT NBR": "CPT 15") %>%
    mutate(across(DIAG_1:DIAG_15, ~str_remove(., "\\.")))

nrow(augusta)
n_distinct(augusta$PATIENT_NBR)


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


augusta_long_hcc <- augusta_long %>% left_join(hcc, by = c("code" = "Code"))

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

