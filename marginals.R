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
                    start = list(shape = 10, scale = 150000))@coef

# Kolmogorov-Smirnov test
ks.test(x = loss$loss, y = plnorm, lognorm.coef[[1]], lognorm.coef[[2]])$p.value
ks.test(x = loss$loss, y = pweibull, weibull.coef[[1]], weibull.coef[[2]])$p.value
ks.test(x = loss$loss, y = ppareto, pareto.coef[[1]], pareto.coef[[2]])$p.value

# lognorm is the best

# one can use more fancy and sophisticated fuction:

performKsTest <- function(pdf, cdf, param1, param2) {
    # Perform ks.test for given loss$loss taking into account data trancation
    # for given pdf, cdf, and initial parameters.
    #
    # Args:
    #   pdf: the probability density fuctions (dlnorm, ...)
    #   cdf: the comulative distribution fuction (plnorm, ...)
    #   param1: initial value for the first parameter 
    #   param2: initial value for the second parameter
    #
    # Returns:
    #   output of ks.test
    #
    # Further development: (for general case)
    #   1. case for any number of parameters using ...
    #   2. pass as parameters sample, limits and if data censored or not
    
    density <- function(arg1, arg2) {
        sum <- 0
        for(i in 1:length(loss$loss)) {
            if(loss$censored[i] == 0) {
                # case if observation is not censored f(x_i)
                sum <- sum - log(pdf(x = loss$loss[i], arg1, arg2))
            } else {
                # case if obervations censored P(X>x_i)
                sum <- sum - log(1 - cdf(q = loss$limit[i], arg1, arg2))
            }
        }
        sum
    }
    
    parameters <- mle(minuslogl = density, 
                      start = list(arg1 = param1, arg2 = param2))@coef
    
    ks.test(x = loss$loss, y = cdf, parameters[[1]], parameters[[2]])
    
}
performKsTest(pdf = dlnorm, cdf = plnorm, param1 = 5, param2 = 5)
performKsTest(pdf = dweibull, cdf = pweibull, param1 = 0.6348689,
              param2 = 28271.23)
performKsTest(pdf = dpareto, cdf = ppareto, param1 = 10, param2 = 150000)


#_________________________ALAE:
# the data is not censored, so we can use fitdistr
coef1 <- fitdist(data = loss$alae, distr = plnorm)
ks.test(loss$alae, plnorm, coef1$estimate[[1]], coef1$estimate[[2]])$p.value

coef2 <- fitdist(data = loss$alae, distr = pweibull)
ks.test(loss$alae, pweibull, coef2$estimate[[1]], coef2$estimate[[2]])$p.value

coef3 <- fitdist(data = loss$alae, distr = ppareto, 
                 start = list(scale = 10, shape = 10))
ks.test(loss$alae, ppareto, coef3$estimate[[1]], coef3$estimate[[2]])$p.value


# alae.lnorm.param <- fitdist(data = loss$alae, distr = plnorm)
# ks.test(loss$alae, plnorm, alae.lnorm.param$estimate[[1]],
#         alae.lnorm.param$estimate[[2]])$p.value
# 
# alae.weibull.param <- fitdist(data = loss$alae, distr = pweibull)
# ks.test(loss$alae, pweibull, alae.weibull.param$estimate[[1]], 
#         alae.weibull.param$estimate[[2]])$p.value
# 
# alae.pareto.param <- fitdist(data = loss$alae, distr = ppareto, 
#                              start = list(scale = 10, shape = 10))
# ks.test(loss$alae, ppareto, alae.pareto.param$estimate[[1]],
#         alae.pareto.param$estimate[[2]])$p.value

# lognorm as well is the best

# TODO:
# 1. hist + qqplot
# 2. ecdf vs lognorm cdf 
# 3. copula mle