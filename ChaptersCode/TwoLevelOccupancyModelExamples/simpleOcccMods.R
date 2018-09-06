library(rstan)
library(data.table)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## @knitr simSimpleData
## Simulate true occupancy states

set.seed(1234)
nSites <- 250
nSurveys <- 10
p <- 0.6
psi <- 0.4

## simulate site-level occupancy
z <- rbinom( nSites, 1, psi)

## simulate  sample-level detections 
y <- matrix( NA, nSites, nSurveys);
for (site in 1:nSites)
  y[site,] <- rbinom( nSurveys, 1, z[site] * p);


simpleModelInput <- list(nSites = nSites,
                         nSurveys = nSurveys,
                         y = y)

## @knitr fitOcc
## Fit the model in stan

fitWide <- stan('../simpleModelExamples/occupancy.stan',
                data = simpleModelInput)

fitWide

## @knitr fitOccMu

## Fit the model using parameters on the mu scale
fitWideMu <- stan('../simpleModelExamples/occupancyMu.stan',
                  data = simpleModelInput)

fitWideMu

## @knitr convertToBin

## Now use summizations of the data
ySum <- apply(y, 1, sum)
k <- rep(dim(y)[2], dim(y)[1])

stanSumationData <- list(
    nSites = nSites,
    y = ySum,
    k = k)

## @knitr fitOccMuBi

fitWideMuBinomial <- stan('../simpleModelExamples/occupancyMuBinomial.stan',
                          data = stanSumationData)

fitWideMuBinomial

## @knitr convertToLong
## Convert to long format
library(data.table)
yDT <- data.table(y)
yDT[ , site := factor(1:nSites)]
yDT[ , z := z]

yDTlong <- melt(yDT, id.vars = c("site", "z"),
                variable.name = "visit",
                value.name = "y")
yDTlong <- yDTlong[ order(site)]
yDTlong[ , index := 1:nrow(yDTlong)]
yDTlong[ , zObs := ifelse(sum(y) >= 1, 1, 0), by = .(site)]

yDTlongSummary <-
    yDTlong [ , .(any_seen = ifelse(sum(y) >0, 1, 0), .N), by = site]

X_psi <-  model.matrix( ~ 1, data = yDTlongSummary)
n_site <- yDTlong[ , length(unique(site))]
m_psi <- dim(X_psi)[2]

total_surveys <- dim(yDTlong)[1]

X_p <- model.matrix( ~ 1, data = yDTlong)
m_p <- dim(X_p)[2]
startStop <- yDTlong[ , .(start_idx = min(index),
                          end_idx = max(index)), by = site]
site = yDTlong[ , as.numeric(site)]

stan_d <- list(
    n_site = n_site,
    m_psi = m_psi,
    X_psi = X_psi,
    total_surveys = total_surveys, 
    m_p = m_p, 
    X_p = X_p, 
    site = site,
    y = yDTlong[ , y],
    start_idx = startStop[ , start_idx], 
    end_idx = startStop[ , end_idx], 
    any_seen = yDTlongSummary[ , any_seen], 
    n_survey = yDTlongSummary[ , N])

## @knitr fitLongForm
## Fit long form model 
fitWideMuLong <- stan('../simpleModelExamples/bernoulli-occupancy.stan',
                      data= stan_d,
                      chains = 4, iter = 2000)

print(fitWideMuLong, pars = c("beta_psi", "beta_p"))

