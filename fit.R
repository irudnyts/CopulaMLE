library("copula")
library("actuar")
library("stats4")
library("fitdistrplus")

data(loss) # load data
summary(loss) # summary of loss data
boxplot(loss[,c(1, 2)]) # no use :)

sum(loss$censored) # number of censored data
loss[(loss$loss >= loss$limit), c(1,3)]
# there are 186 cases, when loss$loss >= loss$limit:
#   34 right-censored values
#   148 without limit (limit = -99, that is why loss > limit)

#_______________________________Fitting marginals_______________________________
#_________________________LOSS:

# if data is not censored use code below:
# LogNormCoef2 <- fitdist(loss$loss, distr = plnorm)
# ks.test(x = loss$loss, y = plnorm, LogNormCoef2$estimate[[1]],
#        LogNormCoef2$estimate[[2]])
# ParetoCoef <- fitdist(data = loss$loss, distr = ppareto, 
#                       start = list(shape = 1, scale = 1))
# ks.test(x = loss$loss, y = ppareto, ParetoCoef$estimate[[1]], 
#         ParetoCoef$estimate[[2]])


# Densities:

lognormDensity <- function(meanlog, sdlog) {
    sum <- 0
    for(i in 1:length(loss$loss)) {
        if(loss$censored[i] == 0) {
            # case if observation is not censored f(x_i)
            sum <- sum - log(dlnorm(x = loss$loss[i],
                                    meanlog = meanlog, sdlog = sdlog))
        } else {
            # case if obervations censored P(X>x_i)
            sum <- sum - log(1 - plnorm(q = loss$limit[i],
                                        meanlog = meanlog, sdlog = sdlog))
        }
    }
    sum
}

weibullDensity <- function(shape, scale) {
    sum <- 0
    for(i in 1:length(loss$loss)) {
        if(loss$censored[i] == 0) {
            # case if observation is not censored f(x_i)
            sum <- sum - log(dweibull(x = loss$loss[i],
                                      shape = shape, scale = scale))
        } else {
            # case if obervations censored P(X>x_i)
            sum <- sum - log(1 - pweibull(q = loss$limit[i],
                                          shape = shape, scale = scale))
        }
    }
    sum
}

paretoDensity <- function(shape, scale) {
    sum <- 0
    for(i in 1:length(loss$loss)) {
        if(loss$censored[i] == 0) {
            # case if observation is not censored f(x_i)
            sum <- sum - log(dpareto(x = loss$loss[i],
                                      shape = shape, scale = scale))
        } else {
            # case if obervations censored P(X>x_i)
            sum <- sum - log(1 - ppareto(q = loss$limit[i],
                                          shape = shape, scale = scale))
        }
    }
    sum
}


# Coeficients:
# initial values one can figure out using fucntion
# fitdist(data = loss$loss, distr = pweibull)

lognorm.coef <- mle(minuslogl = lognormDensity,
                   start = list(meanlog = 5, sdlog = 5))@coef

weibull.coef <- mle(minuslogl = weibullDensity,
                    start = list(shape = 0.6348689, scale = 28271.23))@coef

pareto.coef <- mle(minuslogl = paretoDensity,
                    start = list(shape = 10, scale = 15000))@coef

# Kolmogorov-Smirnov test
ks.test(x = loss$loss, y = plnorm, lognorm.coef[[1]], lognorm.coef[[2]])$p.value
ks.test(x = loss$loss, y = pweibull, weibull.coef[[1]], weibull.coef[[2]])$p.value
ks.test(x = loss$loss, y = ppareto, pareto.coef[[1]], pareto.coef[[2]])$p.value

# lognorm is the best

#_________________________ALAE:
# the data is not censored, so we can use fitdistr
coef1 <- fitdist(data = loss$alae, distr = plnorm)
ks.test(loss$alae, plnorm, coef1$estimate[[1]], coef1$estimate[[2]])$p.value

coef2 <- fitdist(data = loss$alae, distr = pweibull)
ks.test(loss$alae, pweibull, coef2$estimate[[1]], coef2$estimate[[2]])$p.value

coef3 <- fitdist(data = loss$alae, distr = ppareto, 
                 start = list(scale = 10, shape = 10))
ks.test(loss$alae, ppareto, coef3$estimate[[1]], coef3$estimate[[2]])$p.value

# lognorm as well is the best

# TODO:
# 1. hist + qqplot
# 2. ecdf vs lognorm cdf 
# 3. copula mle