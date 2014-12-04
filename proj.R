#install.packages("copula")
library("copula")
set.seed(1)

n <- 1000 # number of observation
p <- 0.3 # mixing parameter
a <- 5 # parameter of 1st copula
b <- 10 # parameter of second copula
I <- rbinom(n = n, size = 1, prob = p) # Bernoili r.v.

#_________________________________Gumbel copula_________________________________
data <-  I * rCopula(n = n, copula = gumbelCopula(param = a, dim = 2)) +
    (1 - I) * rCopula(n = n, copula = gumbelCopula(param = b, dim = 2))

density <- function(arg) {
    sum(log(arg[1] * dCopula(u = data, copula = 
                                 gumbelCopula(param = arg[2], dim = 2))
         + (1 - arg[1]) * dCopula(u = data, copula = 
                                 gumbelCopula(param = arg[3], dim = 2))))
}

optim(par = c(0, 1, 1), fn = density, control = list(fnscale = -1),
      lower = c(0,1,1), upper = c(1, Inf, Inf))

#_________________________________Frank copula__________________________________
data <-  I * rCopula(n = n, copula = frankCopula(param = a, dim = 2)) +
    (1 - I) * rCopula(n = n, copula = frankCopula(param = b, dim = 2))

density <- function(arg) {
    sum(log(arg[1] * dCopula(u = data, copula = 
                                 frankCopula(param = arg[2], dim = 2))
            + (1 - arg[1]) * dCopula(u = data, copula = 
                                     frankCopula(param = arg[3], dim = 2))))
}

optim(par = c(0.1, 1, 1), fn = density, control = list(fnscale = -1),
      lower = c(0,1,1), upper = c(1, Inf, Inf))

#________________________________Clayton copula_________________________________
data <-  I * rCopula(n = n, copula = claytonCopula(param = a, dim = 2)) +
    (1 - I) * rCopula(n = n, copula = claytonCopula(param = b, dim = 2))

density <- function(arg) {
    sum(log(arg[1] * dCopula(u = data, copula = 
                                 claytonCopula(param = arg[2], dim = 2))
            + (1 - arg[1]) * dCopula(u = data, copula = 
                                     claytonCopula(param = arg[3], dim = 2))))
}

optim(par = c(0.1, 1, 1), fn = density, control = list(fnscale = -1),
      lower = c(0, 1, 1), upper = c(1, Inf, Inf))
#________________________________Mixed 2 copulas________________________________
data <-  I * rCopula(n = n, copula = claytonCopula(param = a, dim = 2)) +
    (1 - I) * rCopula(n = n, copula = gumbelCopula(param = b, dim = 2))

density <- function(arg) {
    sum(log(arg[1] * dCopula(u = data, copula = 
                                 claytonCopula(param = arg[2], dim = 2))
            + (1 - arg[1]) * dCopula(u = data, copula = 
                                     gumbelCopula(param = arg[3], dim = 2))))
}
optim(par = c(0.1, 1, 1), fn = density, control = list(fnscale = -1),
      lower = c(0,1,1), upper = c(1, Inf, Inf))
#_______________________________General function________________________________

copula.mle <- function(sample, copula1, copula2, lower, upper) {
    # Maximum likelihood estimator for mixed copula
    #
    # Args:
    #   sample: a matrix of 2-dimenssional random sample from copula (i.e.
    #       values between 0 and 1, columns are univariate random sample from
    #       uniform distribution)
    #   copula1: the class of the first copula (gumbelCopula etc.)
    #   copula2: the class of the second copula
    #   lower: the vector of lower boundaries of parameters for copula1 and
    #       copula2 respectivly
    #   upper: the vector of upper boundaries of parameters for copula1 and
    #       copula2 respectivly
    #
    # Returns:
    #   the vector of parameters p, alpha and beta respectivly
    
    # definition of mixed copula density fucntion
    density <- function(arg) {
        sum(log(arg[1] * dCopula(u = data, copula = 
                                     copula1(param = arg[2], dim = 2))
                + (1-arg[1]) * dCopula(u = data, copula = 
                                     copula2(param = arg[3], dim = 2))))
    }
    # perform the optimization and subset estimated parameters
    optim(par = c(0, 1, 1), fn = density, control = list(fnscale = -1),
          lower = c(0, lower), upper = c(1, upper))$par
}

# example Gumbel-Gumbel
data <-  I * rCopula(n = n, copula = gumbelCopula(param = a, dim = 2)) +
    (1 - I) * rCopula(n = n, copula = gumbelCopula(param = b, dim = 2))
copula.mle(sample = data, copula1 = gumbelCopula, copula2 = gumbelCopula, 
     lower = c(1, 1), upper = c(Inf, Inf))
# example Frank-Gumbel
data <-  I * rCopula(n = n, copula = frankCopula(param = a, dim = 2)) +
    (1 - I) * rCopula(n = n, copula = gumbelCopula(param = b, dim = 2))
copula.mle(sample = data, copula1 = frankCopula, copula2 = gumbelCopula, 
    lower = c(1, 1), upper = c(Inf, Inf))