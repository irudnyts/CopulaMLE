#install.packages("copula")
library("copula")
set.seed(1)

n <- 1000
p <- 0.3
a <- 5
b <- 10
I <- rbinom(n = n, size = 1, prob = p)
I.2d <- matrix(rep(I, 2), ncol = 2, byrow = F)

#_________________________________Gumbel copula_________________________________
data <-  I.2d * rCopula(n = n, copula = gumbelCopula(param = a, dim = 2)) +
    (1 - I.2d) * rCopula(n = n, copula = gumbelCopula(param = b, dim = 2))

density <- function(arg) {
    sum(log(arg[1] * dCopula(u = data, copula = 
                                 gumbelCopula(param = arg[2], dim = 2))
         + (1-arg[1])*dCopula(u = data, copula = 
                                 gumbelCopula(param = arg[3], dim = 2))))
}

optim(par = c(0.1, 1, 1), fn = density, control=list(fnscale=-1),
      lower = c(0,1,1), upper = c(1, Inf, Inf))
#_________________________________Frank copula__________________________________
data <-  I.2d * rCopula(n = n, copula = frankCopula(param = a, dim = 2)) +
    (1 - I.2d) * rCopula(n = n, copula = frankCopula(param = b, dim = 2))

density <- function(arg) {
    sum(log(arg[1] * dCopula(u = data, copula = 
                                 frankCopula(param = arg[2], dim = 2))
            + (1-arg[1])*dCopula(u = data, copula = 
                                     frankCopula(param = arg[3], dim = 2))))
}

optim(par = c(0.1, 1, 1), fn = density, control=list(fnscale=-1),
      lower = c(0,1,1), upper = c(1,10,20))
#________________________________Clayton copula_________________________________
data <-  I.2d * rCopula(n = n, copula = claytonCopula(param = a, dim = 2)) +
    (1 - I.2d) * rCopula(n = n, copula = claytonCopula(param = b, dim = 2))

density <- function(arg) {
    sum(log(arg[1] * dCopula(u = data, copula = 
                                 claytonCopula(param = arg[2], dim = 2))
            + (1-arg[1])*dCopula(u = data, copula = 
                                     claytonCopula(param = arg[3], dim = 2))))
}

optim(par = c(0.1, 1, 1), fn = density, control=list(fnscale=-1),
      lower = c(0,1,1), upper = c(1,10,20))
#________________________________Mixed 2 copulas________________________________
data <-  I.2d * rCopula(n = n, copula = claytonCopula(param = a, dim = 2)) +
    (1 - I.2d) * rCopula(n = n, copula = gumbelCopula(param = b, dim = 2))

density <- function(arg) {
    sum(log(arg[1] * dCopula(u = data, copula = 
                                 claytonCopula(param = arg[2], dim = 2))
            + (1-arg[1])*dCopula(u = data, copula = 
                                     gumbelCopula(param = arg[3], dim = 2))))
}
optim(par = c(0.1, 1, 1), fn = density, control=list(fnscale=-1),
      lower = c(0,1,1), upper = c(1,10,20))
#_______________________________General function________________________________

pmle <- function(sample, copula1, copula2, lower, upper) {
    
    density <- function(arg) {
        sum(log(arg[1] * dCopula(u = data, copula = 
                                     copula1(param = arg[2], dim = 2))
                + (1-arg[1])*dCopula(u = data, copula = 
                                     copula2(param = arg[3], dim = 2))))
    }
    optim(par = c(0.1, 1, 1), fn = density, control=list(fnscale=-1),
          lower = c(0, lower), upper = c(1, upper))
}

pmle(sample = data, copula1 = gumbelCopula, copula2 = gumbelCopula, 
     lower = c(1, 1), upper = c(Inf, Inf))



# testing fucntions in R
# # max likelyhood estimator
# x <- rpois(n = 1000, lambda = 2)
# optimize(f = LL, interval = c(0,200), maximum = TRUE)
# LL <- function(lambda) {
#         sum(dnorm(x, lambda, log = TRUE))
# }
# 
# # test optimizer for 2 dim fnct
# BB <- function(x) {
#     x[1] * x[2]
# }
# optim(par = c(1,1),fn = BB, control=list(fnscale=-1))
# 
# # likelyhood for simple copula
# x <- rCopula(n = 10000, copula = gumbelCopula(param = 2, dim = 2))
# 
# LL <- function(lam) {
#     sum(dCopula(u = x, copula = gumbelCopula(param = lam, dim = 2), log = TRUE))
# }
# 
# optimize(f = LL, maximum = TRUE, interval = c(0, 10))

