
library(brms)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

global_allometry = data.table::fread("./indir/GlobalAllometricDatabase.csv")
global_allometry = global_allometry %>% filter(Biogeographic_zone == "Nearctic")
global_allometry$Biome %>% table

fit1 <- brm(D ~ t2(sqrt(H)) + t2(CD) + (1|g|Biome) +
              (1|r|Functional_type)
              #(1|s|Biogeographic_zone)
            , data = global_allometry, family = lognormal())
fit1
bayes_R2(fit1)
predictions=predict(fit1)

plot(fit1)
crown_marginal_effects = conditional_effects(fit1, "CD")
plot(crown_marginal_effects)
