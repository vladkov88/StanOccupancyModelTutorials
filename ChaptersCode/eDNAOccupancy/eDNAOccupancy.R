## Code for Chapter on eDNA Occupancy modeling

## Load required packages
library(rstan)
library(data.table)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(12345)
## Simulate data 
nSites       = 100
nSamples     = 10
k            = 8
totalSamples = nSites * nSamples
psi   = 0.8
theta = 0.6
p     = 0.3

siteLevelData <- data.table(
    site  = 1:nSites,
    zSim  = rbinom(nSites, 1, psi))
siteLevelData

sampleLevelData <-  data.table(
    site   = rep(1:nSites, each = nSamples),
    sample = rep(1:nSamples, times = nSites),
    theta = theta, 
    p = p,
    k = k
)

setkey(siteLevelData, "site")
setkey(sampleLevelData, "site")

dataUse <- siteLevelData[ sampleLevelData]
dataUse
## was eDNA captured in sample?
dataUse[ , aSim := zSim * rbinom(n = nrow(dataUse), size = 1, p = theta)]
## was eDNA detected in sample?
dataUse[ , y := rbinom(n = nrow(dataUse), size = k, p = aSim * p)]
## Calculate observed a and z
dataUse[ ,  aObs := ifelse(y > 0, 1, 0)]
dataUse[ ,  zObs := ifelse(sum(aSim) > 0, 1, 0), by = .(site)]
dataUse[ aObs != aSim,]



psi
dataUse[ , mean(zSim)]
dataUse[ , mean(zObs)]

theta
dataUse[ zObs > 0, mean(aSim)]
dataUse[ zObs > 0, mean(aObs)]

p
dataUse[ aObs > 0, mean(y/k)]


siteSummary <- dataUse[ , .(a = ifelse(sum(y) >0, 1, 0), .N, sumY = sum(y),
                            z = ifelse(sum(zObs) >0, 1, 0)), by = site]
siteSummary
## Create predictor matricies 
Xpsi <- model.matrix( ~ 1, data = siteSummary)
nPsiCoef = dim(Xpsi)[2]

Vp <- model.matrix( ~ 1, data = dataUse)
nPCoef <- dim(Vp)[2]
## nThetaCoef
## Wtheta
dataUse[ , index:=1:nrow(dataUse)]
startStop <- dataUse[ , .(start_idx = min(index),
                          end_idx = max(index)), by = site]
y <- dataUse[ , y]
aObs <- dataUse[ , aObs]
k <- dataUse[ , k]
zObs = siteSummary[ , z]

startIndex <- startStop[ , start_idx]
endIndex   <- startStop[ , end_idx]

stanData <- list(nSites = nSites,
                 nPsiCoef = nPsiCoef,
                 Xpsi = Xpsi,
                 totalSamples = nrow(dataUse),
                 nPCoef = nPCoef,
                 nThetaCoef = nPCoef,
                 Vp = Vp,
                 Wtheta = Vp,
                 y = y,
                 k = k,
                 zObs = siteSummary[ , z],
                 startIndex = startIndex,
                 endIndex = endIndex,
                 nSamples = siteSummary[ , N],
                 aObs = dataUse[, aObs]
                 )


fitWideMuLong <- stan('./eDNAoccupancy.stan',
                      data= stanData,
                      chains = 4, iter = 2000)

traceplot(fitWideMuLong, pars = c("beta_psi", "beta_theta", "beta_p",  "lp__"), inc_warmup = TRUE)
traceplot(fitWideMuLong, pars = c("beta_psi", "beta_theta", "beta_p",  "lp__"), inc_warmup = FALSE)

## fitWideMuLong
plot(fitWideMuLong, pars = c("beta_psi", "beta_theta", "beta_p")) 

gc()
print(fitWideMuLong, pars = c("beta_psi", "beta_theta", "beta_p",  "lp__"))

psi
dataUse[ , mean(zSim)]
dataUse[ , mean(zObs)]
plogis(1.32)


theta
dataUse[ zObs > 0, mean(aSim)]
dataUse[ zObs > 0, mean(aObs)]
plogis(0.43)

p
dataUse[ aObs > 0, mean(y/k)]
plogis(-0.77)

    



