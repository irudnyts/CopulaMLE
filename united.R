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
    
    density <- function(arg1, arg2) {
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
    
    parameters <- mle(minuslogl = density, 
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
     main = "Histofram and fitted densities for loss", xlab = "loss")
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
     main = "Histofram and fitted densities for alae", xlab = "alae")
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

#_________________________Empirical cdf:
loss.ecdf <- ecdf(loss$loss)
alae.ecdf <- ecdf(loss$alae)

#________________________________Copula fitting_________________________________
#_________________________Pseudo sample (parametrically):

pseudo.sample.parametric <- data.frame(
    loss = plnorm(loss$loss, loss.lnorm[[2]][[1]], loss.lnorm[[2]][[2]]),
    alae = plnorm(loss$alae, alae.lnorm[[2]][[1]], alae.lnorm[[2]][[2]]))

pseudo.sample.non.parametric <- data.frame(loss = loss.ecdf(loss$loss),
                                           alae = alae.ecdf(loss$alae))

par(mfrow = c(1,1))
plot(pseudo.sample.non.parametric$loss, pseudo.sample.non.parametric$alae)
plot(pseudo.sample.parametric$loss, pseudo.sample.parametric$alae)


#_________________________Pseudo sample (non-parametrically):


# TODO:
# 1. hist + qqplot
# 2. ecdf vs lognorm cdf 
# 3. copula mle
# 4. calculate tau, Spearman's and Pearson's rho