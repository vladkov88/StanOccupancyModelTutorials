### Code for chapter 1

x <- 1:10
y <- rbinom(length(x), size = 1, prob = 0.5)
            
## @knitr glmExample
glm(y ~ x, family = 'binomial')

## @knitr glmProbit
glm( y ~ x, family = binomial(link = 'probit'))

## @knitr logisticRegression
library(rstan)
options(mc.cores = parallel::detectCores())

### simulate and examine raw data
n_obs = 20
x_label <- factor(rep(c("a", "b"), times = n_obs)) # used to Id variables
x_prob  <- rep(c(0.25, 0.75), times = n_obs) # used to simulate data
y <- rbinom(n = n_obs * 2, size = 1, p = x_prob)
aggregate(y, by = list(x_label), FUN = mean)

### First, we fit a GLM using base R
glm_out <- glm(y ~ x_label, family = 'binomial')
summary(glm_out)

### Second, we fit the model using Stan:
stan_data <- list(N = length(y),
                  x = as.numeric(x_label) - 1,
                  y = y)

stan_model <- stan_model(file = "logistic_stan.stan")

stan_out <- sampling(stan_model, data = stan_data, chains = 4, iter = 500)

### Useful functions for looking at Stan outputs include:
print(stan_out)
traceplot(stan_out)
traceplot(stan_out, inc_warmup = TRUE)
plot(stan_out)
pairs(stan_out)


## @knitr logit_target
## Complie the model
stan_model_target <- stan_model("logistic_target_stan.stan")

## Sample from the model
stan_out_target <- sampling(stan_model_target,
                            data = stan_data, chains = 4,
                            iter = 1000)

## Useful functions for looking at Stan outputs include:
print(stan_out_target, pars = c("alpha", "beta", "lp__"))
traceplot(stan_out_target, pars = c("alpha", "beta", "lp__"))
traceplot(stan_out_target, inc_warmup = TRUE, pars = c("alpha", "beta", "lp__"))
plot(stan_out_target, pars = c("alpha", "beta", "lp__"))


## @knitr demoGLMinput
## Simulate data in long format
set.seed(1223)
xSim <- rep(c(0.25, 0.75), each = 14)
x <- rep(c("a", "b"), each = 14)
y <- rbinom(n = length(xSim), size = 1, prob = xSim)
data_long <-
    data.frame(x = x,
               y = factor(y, labels = c("fail", "success")))

### cast the data to wide format using tidyverse
library(tidyverse)
data_wide <-
    data_long %>%
    group_by(x,y) %>%
    summarize(N = n()) %>% 
    spread( y, N) %>%
    mutate(Total = success + fail,
           success_proportion = success / (success + fail))
data_wide

### Compare the three methods for fitting a logistic regression in R
glm(y ~ x, family = 'binomial', data = data_long)

glm(cbind(fail, success) ~ x, family = 'binomial', data = data_wide)

glm(success_proportion ~ x, family = 'binomial', data = data_wide,
    weights = Total)

## @knitr demoMM
df <- data.frame(city = c("La Crosse", "St. Louis", "Cairo"))
model.matrix( ~ city, data = df)

## @knitr demoMM1
df <- data.frame(city = c("La Crosse", "St. Louis", "Cairo"))
model.matrix( ~ city -1, data = df)

## @knitr demoMMm1b
df <- expand.grid(city = c("La Crosse", "St. Louis", "Cairo"),
                  month = c("May", "June"))
model.matrix( ~ city + month -1, data = df)
model.matrix( ~ month + city -1, data = df)

## @knitr numericMatrix
df1 <- expand.grid(month = c( 5, 6, 7))
model.matrix( ~ month, data = df1)

## @knitr factorMatrix
df2 <- expand.grid(month = factor(c( 5, 6, 7)))
model.matrix( ~ month, data = df2)

## @knitr m1
set.seed(12351513)
df_occ <-
    data.frame(occ = factor(sample(rep(c("yes", "no"), each = 10))),
               site = factor(rep(c("lake", "river"), each = 10)))

## @knitr m1glm
summary(glm(occ ~ site, data = df_occ, family = 'binomial'))

## @knitr m1mm
### base line 
as.numeric(factor(df_occ$occ))
### need -1 to create a vecor of zeros and ones
as.numeric(factor(df_occ$occ)) - 1
