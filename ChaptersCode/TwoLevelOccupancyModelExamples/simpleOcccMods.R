library(rstan)
library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## @knitr simSimpleData
## Simulate true occupancy states

set.seed(1234)
## number of sampling units 
n_sampling_units <- 250
## visits per sampling unit 
n_surveys <- 10
p <- 0.6
psi <- 0.4

## simulate site-level occupancy
z <- rbinom( n_sampling_units, 1, psi)

## simulate  sample-level detections 
y <- matrix( NA, n_sampling_units, n_surveys);
for (site in 1:n_sampling_units){
    y[site,] <- rbinom( n_surveys, 1, z[site] * p)
}


simple_model_input <- list(n_sampling_units = n_sampling_units,
                           n_surveys = n_surveys,
                           y = y)

## @knitr fitOcc
## Fit the model in stan

model_wide <- stan_model('./occupancy.stan')

fit_wide <- sampling(model_wide, simple_model_input)
fit_wide

## @knitr fitOccMu
## Fit the model using parameters on the mu scale
model_wide_mu <- stan_model('./occupancy_mu.stan')

fit_wide_mu <- sampling(model_wide_mu, data = simple_model_input)
fit_wide_mu

## @knitr convertToBin

## Now use summizations of the data
y_sum <- apply(y, 1, sum)
k <- rep(dim(y)[2], dim(y)[1])

stan_sumation_data <-
    list(
        n_sampling_units = n_sampling_units,
        y = y_sum,
        k = k)

## @knitr fitOccMuBi
stan_model_mu_binomial  <- stan_model('./occupancy_mu_binomial.stan')

fit_mu_binomial <- sampling(stan_model_mu_binomial,
                            data = stan_sumation_data)

fit_mu_binomial

## @knitr convertToLong
## Convert to long format
y_df <- 
    y %>%
    as.data.frame() %>%
    mutate(site = factor(1:n_sampling_units),
           z_obs = ifelse(rowSums(y) >0, 1, 0))

y_long <-
    y_df %>%
    gather(survey, y,  -site, -z_obs) %>%
    arrange(site) %>%
    mutate(index = 1:nrow(y_long))

head(y_long, 12)
tail(y_long, 12)

y_long_summary <-
    y_long %>%
    group_by(site) %>%
    summarize(any_seen = ifelse(sum(y) >0, 1, 0),
              N = n())

X_psi  <-  model.matrix( ~ 1, data = y_long_summary)
n_site <- y_long %>% pull(site) %>% unique() %>% length()
m_psi  <- ncol(X_psi)

total_surveys <- nrow(y_long)

X_p <- model.matrix( ~ 1, data = y_long)
m_p <- ncol(X_p)
    
start_end <-
    y_long %>%
    group_by(site) %>%
    summarize(start_idx = min(index),
              end_idx = max(index))

site = y_long %>% pull(site) %>% as.numeric()

stan_d <- list(
    n_site = n_site,
    m_psi = m_psi,
    X_psi = X_psi,
    total_surveys = total_surveys, 
    m_p = m_p, 
    X_p = X_p, 
    site = site,
    y = y_long %>% pull(y),
    start_idx = start_end %>% pull( start_idx),
    end_idx = start_end %>% pull( end_idx),  
    any_seen = y_long_summary %>% pull(any_seen), 
    n_survey = y_long_summary %>% pull(N)
    )

## @knitr fitLongForm
## Fit long form model 
stan_model_bern <- stan_model("./occupancy_bernoulli.stan")

fit_bern <- sampling(stan_model_bern,
                     data= stan_d,
                     chains = 4, iter = 2000)

print(fit_bern, pars = c("beta_psi", "beta_p"))

##              mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
## beta_psi[1] -0.52       0 0.13 -0.78 -0.60 -0.51 -0.43 -0.27  3909    1
## beta_p[1]    0.39       0 0.07  0.25  0.34  0.39  0.43  0.51  3393    1
