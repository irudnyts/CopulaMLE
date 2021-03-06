library("copula")
library("actuar")
library("stats4")
library("fitdistrplus")

options(scipen=999) # disable scientific representation

data(loss) # load data
# summary(loss) # summary of loss data
# boxplot(loss[,c(1, 2)]) # no use :)

sum(loss$censored) # number of censored data
loss[(loss$loss >= loss$limit), c(1,3)]
# there are 186 cases, when loss$loss >= loss$limit:
#   34 right-censored values
#   148 without limit (limit = -99, that is why loss > limit)

#_______________________________Fitting marginals_______________________________
#_________________________LOSS:

fitAndTest <- function(sample, limit, censored, pdf, cdf, param1, param2) {
    # Estimates parameters and perform ks.test for given loss$loss taking into
    # account data censoring for given pdf, cdf, and initial parameters.
    #
    # Args:
    #   sample: a numeric vector of random sample
    #   limit: a numeric vector of limit (-99 means no limit)
    #   censored: a numeric vector: 1 means censored (limit reached)
    #       and 0 otherwise.
    #   pdf: the probability density fuctions (dlnorm, ...)
    #   cdf: the comulative distribution fuction (plnorm, ...)
    #   param1: initial value for the first parameter 
    #   param2: initial value for the second parameter
    #
    # Returns:
    #   list of parameters and output of ks.test
    #
    # Further development: (for general case)
    #   1. case for any number of parameters using ...
    #   2. pass as parameters sample, limits and if data censored or not
    
    logL <- function(arg1, arg2) {
        sum <- 0
        for(i in 1:length(sample)) {
            if(censored[i] == 0) {
                # case if observation is not censored f(x_i)
                sum <- sum - log(pdf(x = sample[i], arg1, arg2))
            } else {
                # case if obervations censored P(X>x_i)
                sum <- sum - log(1 - cdf(q = limit[i], arg1, arg2))
            }
        }
        sum
    }
    
    parameters <- mle(minuslogl = logL, 
                      start = list(arg1 = param1, arg2 = param2))@coef
    
    list(ks.test(x = sample, y = cdf, parameters[[1]], parameters[[2]]),
         parameters)
}

loss.lnorm <- fitAndTest(sample = loss$loss, limit = loss$limit,
                         censored = loss$censored,
                         pdf = dlnorm, cdf = plnorm,
                         param1 = 5, param2 = 5)

loss.weibull <- fitAndTest(sample = loss$loss, limit = loss$limit,
                           censored = loss$censored,
                           pdf = dweibull, cdf = pweibull, 
                           param1 = 0.6348689, param2 = 28271.23)

loss.pareto <- fitAndTest(sample = loss$loss, limit = loss$limit, 
                          censored = loss$censored,
                          pdf = dpareto, cdf = ppareto, param1 = 10, 
                          param2 = 150000)
# Histogram
hist(loss$loss, breaks = 400, freq = F, xlim = c(0, 250000),
     ylim = c(0, 0.00009),
     main = "Histogram and fitted densities for loss", xlab = "loss")
lines(seq(from = 10, to = 250000, by = 50),
      dlnorm(x = seq(from = 10, to = 250000, by = 50),
             loss.lnorm[[2]][[1]], loss.lnorm[[2]][[2]]), col = 2)
lines(seq(from = 10, to = 250000, by = 50),
      dweibull(x = seq(from = 10, to = 250000, by = 50),
               loss.weibull[[2]][[1]], loss.weibull[[2]][[2]]), col = 3)
lines(seq(from = 10, to = 250000, by = 50),
      dpareto(x = seq(from = 10, to = 250000, by = 50),
              loss.pareto[[2]][[1]], loss.pareto[[2]][[2]]), col = 4)
legend("topright", legend  = c("Histogram","Lognormal", "Weibull", "Pareto"),
       col = 1:4, lwd = rep(1, 4))

# QQ plot for lnorm, pareto and weibull (for loss data)
qqplot2(loss$loss,
        qF = function(p)
            qlnorm(p, loss.lnorm[[2]][[1]], loss.lnorm[[2]][[2]]),
        main = "Q-Q Plot")

qqplot2(loss$loss,
        qF = function(p)
            qweibull(p, loss.weibull[[2]][[1]], loss.weibull[[2]][[2]]),
        main = "Q-Q Plot")

qqplot2(loss$loss,
        qF = function(p)
            qpareto(p, loss.pareto[[2]][[1]], loss.pareto[[2]][[2]]),
        main = "Q-Q Plot")

#_________________________ALAE:

alae.lnorm <- fitAndTest(sample = loss$alae, limit = rep(-99, nrow(loss)),
                         censored = rep(0, nrow(loss)),
                         pdf = dlnorm, cdf = plnorm,
                         param1 = 5, param2 = 5)


alae.weibull <- fitAndTest(sample = loss$alae, limit = rep(-99, nrow(loss)),
                         censored = rep(0, nrow(loss)),
                         pdf = dweibull, cdf = pweibull,
                         param1 = 0.741633, param2 = 9994.507167)

alae.pareto <- fitAndTest(sample = loss$alae, limit = rep(-99, nrow(loss)),
                           censored = rep(0, nrow(loss)),
                           pdf = dpareto, cdf = ppareto,
                           param1 = 10, param2 = 10)

# Histogram
hist(loss$alae, breaks = 200, freq = F, ylim = c(0, 0.0002), xlim = c(0, 40000),
     main = "Histogram and fitted densities for alae", xlab = "alae")
lines(seq(from = 10, to = 501863, by = 50),
      dlnorm(x = seq(from = 10, to = 501863, by = 50),
             alae.lnorm[[2]][[1]], alae.lnorm[[2]][[2]]), col = 2)
lines(seq(from = 10, to = 501863, by = 50),
      dweibull(x = seq(from = 10, to = 501863, by = 50),
               alae.weibull[[2]][[1]], alae.weibull[[2]][[2]]), col = 3)
lines(seq(from = 10, to = 501863, by = 50),
      dpareto(x = seq(from = 10, to = 501863, by = 50),
              alae.pareto[[2]][[1]], alae.pareto[[2]][[2]]), col = 4)
legend("topright", legend  = c("Histogram","Lognormal", "Weibull", "Pareto"),
       col = 1:4, lwd = rep(1, 4))

# QQ plot for lnorm, pareto and weibull (for alae data)

qqplot2(loss$alae,
        qF = function(p)
            qlnorm(p, alae.lnorm[[2]][[1]], alae.lnorm[[2]][[2]]),
        main = "Q-Q Plot")

qqplot2(loss$alae,
        qF = function(p)
            qweibull(p, alae.weibull[[2]][[1]], alae.weibull[[2]][[2]]),
        main = "Q-Q Plot")

qqplot2(loss$alae,
        qF = function(p)
            qpareto(p, alae.pareto[[2]][[1]], alae.pareto[[2]][[2]]),
        main = "Q-Q Plot")

#_________________________Empirical cdf:
loss.ecdf <- ecdf(loss$loss)
alae.ecdf <- ecdf(loss$alae)

#________________________________Copula fitting_________________________________

pseudo.sample.parametric <- data.frame(
    loss = plnorm(loss$loss, loss.lnorm[[2]][[1]], loss.lnorm[[2]][[2]]),
    alae = plnorm(loss$alae, alae.lnorm[[2]][[1]], alae.lnorm[[2]][[2]]))

pseudo.sample.non.parametric <- data.frame(loss = loss.ecdf(loss$loss),
                                           alae = alae.ecdf(loss$alae))

par(mfrow = c(1,2))
plot(pseudo.sample.non.parametric$loss, pseudo.sample.non.parametric$alae,
     xlab = "loss", ylab = "alae", main = "Non-parametric pseudo data")
plot(pseudo.sample.parametric$loss, pseudo.sample.parametric$alae,
     xlab = "loss", ylab = "alae", main = "Parametric pseudo data")
par(mfrow = c(1,1))

# estimatiom tau, rho_s, rho for parametric pseudo data
cor(pseudo.sample.parametric$loss, pseudo.sample.parametric$alae,
    method = c("kendall"))
cor(pseudo.sample.parametric$loss, pseudo.sample.parametric$alae,
    method = c("spearman"))
cor(pseudo.sample.parametric$loss, pseudo.sample.parametric$alae,
    method = c("pearson"))

# estimatiom tau, rho_s, rho for non-parametric pseudo data
cor(pseudo.sample.non.parametric$loss, pseudo.sample.non.parametric$alae,
    method = c("kendall"))
cor(pseudo.sample.non.parametric$loss, pseudo.sample.non.parametric$alae,
    method = c("spearman"))
cor(pseudo.sample.non.parametric$loss, pseudo.sample.non.parametric$alae,
    method = c("pearson"))

# all these values are the same. so, for simplicity we will be dealing with
# only parametric data

# the problem, which may occure, F(max(x)) (of empirical cdf) will give value 1,
# wich will give the density equals to Inf 

# estimate parameters for particular copulas
copula.mle <- function(sample, copula1, copula2, initial, lower, upper) {
    # Maximum likelihood estimator for mixed copula
    #
    # Args:
    #   sample: a matrix of 2-dimenssional random sample from copula (i.e.
    #       values between 0 and 1, columns are univariate random sample from
    #       uniform distribution)
    #   copula1: the class of the first copula (gumbelCopula etc.)
    #   copula2: the class of the second copula
    #   intital: the vector of initial values for p, copula1 and copula2
    #       parameters
    #   lower: the vector of lower boundaries of parameters for p, copula1 and
    #       copula2 respectivly
    #   upper: the vector of upper boundaries of parameters for p, copula1 and
    #       copula2 respectivly
    #
    # Returns:
    #   the vector of parameters p, alpha and beta respectivly
    
    # definition of mixed copula logL fucntion
    logL <- function(arg) {
        sum(log(arg[1] * dCopula(u = sample, copula = 
                                     copula1(param = arg[2], dim = 2))
                + (1-arg[1]) * dCopula(u = sample, copula = 
                                           copula2(param = arg[3], dim = 2))))
    }
    # perform the optimization and subset estimated parameters
    optim(par = initial, fn = logL, control = list(fnscale = -1),
          lower = lower, upper = upper)$par
}

#___________________________Estimation and simulating___________________________

gumbel.gumbel.coef <- copula.mle(sample = as.matrix(pseudo.sample.parametric), 
                                 copula1 = gumbelCopula, copula2 = gumbelCopula,
                                 initial = c(0.2, 1.2, 1.2),
                                 lower = c(0, 1, 1), upper = c(1, Inf, Inf))

n <- 10000 # number of observation
p <- gumbel.gumbel.coef[1] # mixing parameter
a <- gumbel.gumbel.coef[2] # parameter of 1st copula
b <- gumbel.gumbel.coef[3] # parameter of second copula
I <- rbinom(n = n, size = 1, prob = p) # Bernoili r.v.
random <- I * rCopula(n = n, copula = gumbelCopula(param = a, dim = 2)) +
    (1 - I) * rCopula(n = n, copula = gumbelCopula(param = b, dim = 2))
cor(random[ ,1], random[ ,2], method = "kendall")
cor(random[ ,1], random[ ,2], method = "spearman")
cor(random[ ,1], random[ ,2], method = "pearson")

frank.frank.coef <- copula.mle(sample = as.matrix(pseudo.sample.parametric), 
                                 copula1 = frankCopula, copula2 = frankCopula,
                                 initial = c(0.2, 1.2, 1.2),
                                 lower = c(0, 0, 0), upper = c(1, Inf, Inf))

n <- 10000 # number of observation
p <- frank.frank.coef[1] # mixing parameter
a <- frank.frank.coef[2] # parameter of 1st copula
b <- frank.frank.coef[3] # parameter of second copula
I <- rbinom(n = n, size = 1, prob = p) # Bernoili r.v.
random <- I * rCopula(n = n, copula = frankCopula(param = a, dim = 2)) +
    (1 - I) * rCopula(n = n, copula = frankCopula(param = b, dim = 2))
cor(random[ ,1], random[ ,2], method = "kendall")
cor(random[ ,1], random[ ,2], method = "spearman")
cor(random[ ,1], random[ ,2], method = "pearson")

clayton.clayton.coef <- copula.mle(sample = as.matrix(pseudo.sample.parametric), 
                               copula1 = claytonCopula, copula2 = claytonCopula,
                               initial = c(0.2, 1.2, 1.2),
                               lower = c(0, 0, 0), upper = c(1, Inf, Inf))

n <- 10000 # number of observation
p <- clayton.clayton.coef[1] # mixing parameter
a <- clayton.clayton.coef[2] # parameter of 1st copula
b <- clayton.clayton.coef[3] # parameter of second copula
I <- rbinom(n = n, size = 1, prob = p) # Bernoili r.v.
random <- I * rCopula(n = n, copula = claytonCopula(param = a, dim = 2)) +
    (1 - I) * rCopula(n = n, copula = claytonCopula(param = b, dim = 2))
cor(random[ ,1], random[ ,2], method = "kendall")
cor(random[ ,1], random[ ,2], method = "spearman")
cor(random[ ,1], random[ ,2], method = "pearson")

fgm.fgm.coef <- copula.mle(sample = as.matrix(pseudo.sample.parametric), 
                                   copula1 = fgmCopula, copula2 = fgmCopula,
                                   initial = c(0.5, 0, 0),
                                   lower = c(0, -1, -1), upper = c(1, 1, 1))

n <- 10000 # number of observation
p <- fgm.fgm.coef[1] # mixing parameter
a <- fgm.fgm.coef[2] # parameter of 1st copula
b <- fgm.fgm.coef[3] # parameter of second copula
I <- rbinom(n = n, size = 1, prob = p) # Bernoili r.v.
random <- I * rCopula(n = n, copula = fgmCopula(param = a, dim = 2)) +
    (1 - I) * rCopula(n = n, copula = fgmCopula(param = b, dim = 2))
cor(random[ ,1], random[ ,2], method = "kendall")
cor(random[ ,1], random[ ,2], method = "spearman")
cor(random[ ,1], random[ ,2], method = "pearson")

fgm.clayton.coef <- copula.mle(sample = as.matrix(pseudo.sample.parametric), 
                           copula1 = fgmCopula, copula2 = claytonCopula,
                           initial = c(0.2, 0.2, 0.2),
                           lower = c(0, -1, 0), upper = c(1, 1, Inf))

n <- 10000 # number of observation
p <- fgm.clayton.coef[1] # mixing parameter
a <- fgm.clayton.coef[2] # parameter of 1st copula
b <- fgm.clayton.coef[3] # parameter of second copula
I <- rbinom(n = n, size = 1, prob = p) # Bernoili r.v.
random <- I * rCopula(n = n, copula = fgmCopula(param = a, dim = 2)) +
    (1 - I) * rCopula(n = n, copula = claytonCopula(param = b, dim = 2))
cor(random[ ,1], random[ ,2], method = "kendall")
cor(random[ ,1], random[ ,2], method = "spearman")
cor(random[ ,1], random[ ,2], method = "pearson")

gumbel.frank.coef <- copula.mle(sample = as.matrix(pseudo.sample.parametric), 
                               copula1 = gumbelCopula, copula2 = frankCopula,
                               initial = c(1, 1.2, 4),
                               lower = c(0, 1, 0), upper = c(1, Inf, Inf))

n <- 10000 # number of observation
p <- gumbel.frank.coef[1] # mixing parameter
a <- gumbel.frank.coef[2] # parameter of 1st copula
b <- gumbel.frank.coef[3] # parameter of second copula
I <- rbinom(n = n, size = 1, prob = p) # Bernoili r.v.
random <- I * rCopula(n = n, copula = gumbelCopula(param = a, dim = 2)) +
    (1 - I) * rCopula(n = n, copula = frankCopula(param = b, dim = 2))
cor(random[ ,1], random[ ,2], method = "kendall")
cor(random[ ,1], random[ ,2], method = "spearman")
cor(random[ ,1], random[ ,2], method = "pearson")

gumbel.coef <- fitCopula(gumbelCopula(dim = 2), pseudo.sample.parametric,
                         method = "mpl")
random <- rCopula(n = n, copula = gumbelCopula(param = gumbel.coef@estimate,
                                               dim = 2))
cor(random[ ,1], random[ ,2], method = "kendall")
cor(random[ ,1], random[ ,2], method = "spearman")
cor(random[ ,1], random[ ,2], method = "pearson")




### !!! ecdf gives back value 1, which gives Inf of log likelihood
### pseudo.sample.non.parametric[1268,]
