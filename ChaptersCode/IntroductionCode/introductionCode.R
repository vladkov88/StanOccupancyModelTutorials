### Code for chapter 1

x <- 1:10
y <- rbinom(length(x))
            
## @knitr glmExample
glm(y ~ x, family = 'binomial')

## @knitr glmProbit
glm( y ~ x, family = binomial(link = 'probit'))

## @knitr logisticRegression
library(rstan)

### simulate and examine raw data
nObs = 20
xLabel <- factor(rep(c("a", "b"), times = nObs)) # used to Id variables
xProb  <- rep(c(0.25, 0.75), times = nObs) # used to simulate data
y <- rbinom(n = nObs * 2, size = 1, p = xProb)
aggregate(y, by = list(xLabel), FUN = mean)

### First, we fit a GLM using base R
glmOut <- glm(y ~ xLabel, family = 'binomial')
summary(glmOut)

### Second, we fit the model using Stan:
stanData <- list(N = length(y),
                 x = as.numeric(xLabel) - 1,
                 y = y)

stanOut <- stanOutO <- stan(file = "logisticStan.stan",
                 data = stanData, chains = 4, iter = 500)

### Useful functions for looking at Stan outputs include:
print(stanOut)
traceplot(stanOut)
traceplot(stanOut, inc_warmup = TRUE)
plot(stanOut)
pairs(stanOut)

## @knitr demoGLMinput
## Simulate data in long format
set.seed(1223)
xSim <- rep(c(0.25, 0.75), each = 14)
x <- rep(c("a", "b"), each = 14)
y <- rbinom(n = length(xSim), size = 1, prob = xSim)
dataLong <- data.frame(x = x, y = factor(y, labels = c("fail", "success")))

### cast the data to wide format using reshape2
library(reshape2)
dataWide <- dcast(dataLong, x ~ y)
dataWide$Total <- with(dataWide, success + fail)
dataWide$successProportion <- with(dataWide, success /Total)
dataWide

### Compare the three methods for fitting a logistic regression in R
glm(y ~ x, family = 'binomial', data = dataLong)

glm(cbind(fail, success) ~ x, family = 'binomial', data = dataWide)

glm(successProportion ~ x, family = 'binomial', data = dataWide,
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
dfocc <- data.frame(occ = factor(sample(rep(c("yes", "no"), each = 10))),
                    site = factor(rep(c("lake", "river"), each = 10)))

## @knitr m1glm
summary(glm(occ ~ site, data = dfocc, family = 'binomial'))

## @knitr m1mm
### base line 
as.numeric(factor(dfocc$occ))
### need -1 to create a vecor of zeros and ones
as.numeric(factor(dfocc$occ)) - 1
