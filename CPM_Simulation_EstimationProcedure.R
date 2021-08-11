# SIMULATION STUDY #
# CAROLINE BIRDROW #

# Set working directory
path <- "/home/birdroci/rsch/project/CPM"  
setwd(path)

# Set seed for reproducibility
set.seed(ITER)

##G-computation with a CPM vs. a Parametric Model 
library("rms")

## Define useful functions
expit <- function(x) {exp(x)/(1 + exp(x))}

pull.one <- function(x, y, fhat)
{
  return(y[max(which(fhat < x))])
}

## Simulation parameters
alpha0 <- -0.3
alpha1 <- 0.7
alpha2 <- 0.5
alpha3 <- 0.7

sigma0 <- -0.6
sigma1 <- 0.25
sigma2 <- 0.2

beta0 <- -2.5
beta1 <- 0.2
beta2 <- 0.1
beta3 <- 0.4
beta4 <- 0.3
beta5 <- 0.6
beta6 <- 0.5
sigmaY <- 0.3

s <- 1
s.d <- 0.3

c <- 4
d <- 1

if(c == 1){
  A1 = 0
  A2 = 0
  A3 = 0
  sa1 = 0
  sa2 = 0
  sa3 = 0
}

if(c == 2){
  A1 = 0
  A2 = 1
  A3 = 1
  sa1 = 0
  sa2 = 1
  sa3 = 1
}

if(c == 3){
  A1 = 0
  A2 = 0
  A3 = 1
  sa1 = 0
  sa2 = 0
  sa3 = 1
}

if(c == 4){
  A1 = 1
  A2 = 1
  A3 = 1
  sa1 = 1
  sa2 = 1
  sa3 = 1
}

## Simulation parameters
N <- 1000
M <- 5000 #Monte-Carlo iterations within standardization
nsim <- 5

## Question 1: Does mean(res) reflect true outcome, my0?

res <- matrix(0, nrow = nsim, ncol = 1)

for (j in 1:nsim){
## Generate data
## First Time Point
l1 <- rnorm(N, mean = 0, sd = s)
a1 <- rbinom(N, 1, expit(alpha0 + alpha1 * l1))

## Second Time Point
eps.l2 <- rnorm(N, mean = 0, sd = s.d)
l2 <- sigma0 + sigma1 * l1 + sigma2 * a1 + eps.l2
l2[l2 < -1] <- 0

if(d == 1){
l2[l2 > -1] <- exp(l2[l2 > -1])
}

if(d == 2){
l2[l2 > -1] <- (l2[l2 > -1])^2
}

if(d == 3){
l2[l2 > -1] <- abs(l2[l2 > -1])
}

a2 <- rbinom(N, 1, expit(alpha0 + alpha1 * l1 + alpha2 * a1 + alpha3 * l2))

## Third Time Point
eps.l3 <- rnorm(N, mean = 0, sd = s.d)
l3 <- sigma0 + sigma1 * l2 + sigma2 * a2 + eps.l3
l3[l3 < -1] <- 0

if(d == 1){
  l3[l3 > -1] <- exp(l3[l3 > -1])
}

if(d == 2){
  l3[l3 > -1] <- (l3[l3 > -1])^ 2
}

if(d == 3){
  l3[l3 > -1] <- abs(l3[l3 > -1])
}

a3 <- rbinom(N, 1, expit(alpha0 + alpha1 * l2 + alpha2 * a2 + alpha3 * l3))

## Final Outcome
eps.Y <- rnorm(N, mean = 0, sd = s)
Ystar <- beta0 + beta1*l1 + beta2*a1 + beta3*l2 + beta4*a2 + beta5*l3 + beta6*a3 + eps.Y
Y <- Ystar
Y[Ystar < -3] <- 0

if(d == 1){
Y[Ystar > -3] <- exp(Ystar[Ystar > -3])
}

if(d == 2){
Y[Ystar > -3] <- (Ystar[Ystar > -3])^2
}

if(d == 3){
Y[Ystar > -3] <- abs(Ystar[Ystar > -3])
}

## Simulated data set:
dat <- data.frame(cbind(l1, a1, l2, a2, l3, a3, Y))
names(dat) <- c("l1", "a1", "l2", "a2", "l3", "a3", "Y")
    
## Estimate sigmas using a CPM
zz <- orm(l2 ~ l1 + a1, data = dat, family = "probit")
alpha.hat1 <- as.numeric(zz$coef[-c((length(zz$coef) - 1):length(zz$coef))])
sigma1.hat1 <- as.numeric(zz$coef[length(zz$coef) - 1]) 
sigma2.hat1 <- as.numeric(zz$coef[length(zz$coef)]) 

## Estimate sigmas using a CPM
xx <- orm(l3 ~ l2 + a2, data = dat, family = "probit")
alpha.hat2 <- as.numeric(xx$coef[-c((length(xx$coef) - 1):length(xx$coef))])
sigma1.hat2 <- as.numeric(xx$coef[length(xx$coef) - 1]) 
sigma2.hat2 <- as.numeric(xx$coef[length(xx$coef)]) 

## Estimate betas using a CPM
reg.y <- orm(Y ~ l1 + a1 + l2 + a2 + l3 + a3, data = dat, family = 'probit')
alpha.hat3 <- as.numeric(reg.y$coef[-c((length(reg.y$coef) - 5):length(reg.y$coef))])
beta1.hat <- as.numeric(reg.y$coef[length(reg.y$coef) - 5]) 
beta2.hat <- as.numeric(reg.y$coef[length(reg.y$coef) - 4]) 
beta3.hat <- as.numeric(reg.y$coef[length(reg.y$coef) - 3]) 
beta4.hat <- as.numeric(reg.y$coef[length(reg.y$coef) - 2]) 
beta5.hat <- as.numeric(reg.y$coef[length(reg.y$coef) - 1]) 
beta6.hat <- as.numeric(reg.y$coef[length(reg.y$coef)]) 

## Goal: Approximate integral of E[Y|A = a, L = l] with respect to f_L(l)

## First for a = 0, a1 = 0, and a2 = 0
tmp1 <- matrix(0, nrow = length(alpha.hat1), ncol = M)
tmp2 <- matrix(0, nrow = length(alpha.hat2), ncol = M)
tmp3 <- matrix(0, nrow = length(alpha.hat3), ncol = M)

sl1 <- sample(dat$l1, size = M, replace = TRUE)

for (i in 1:length(alpha.hat1))
{
  tmp1[i,] <- as.numeric(1 - pnorm(alpha.hat1[i] + sigma1.hat1*sl1 + sigma2.hat1*sa1, mean = 0, sd = 1))
}
tmp1 <- rbind(0, tmp1)

sl2 <- rep(0, M)
for (i in 1:M)
{
  F.hat1 <- tmp1[,i]
  sl2[i] <- pull.one(x = runif(1, 0, 1), y = sort(unique(dat$l2)), fhat = F.hat1)
}
sl2

for (n in 1:length(alpha.hat2))
{
  tmp2[n,] <- as.numeric(1 - pnorm(alpha.hat2[n] + sigma1.hat2*sl2 + sigma2.hat2*sa2, mean = 0, sd = 1))
}

tmp2 <- rbind(0, tmp2)

sl3 <- rep(0, M)
for (i in 1:M)
{
  F.hat2 <- tmp2[,i]
  sl3[i] <- pull.one(x = runif(1, 0, 1), y = sort(unique(dat$l3)), fhat = F.hat2)
}
sl3

for (k in 1:length(alpha.hat3))
{
  tmp3[k,] <- as.numeric(1 - pnorm(alpha.hat3[k] + beta1.hat*sl1  + beta2.hat*sa1 + beta3.hat*sl2 + beta4.hat*sa2 + beta5.hat*sl3 + beta6.hat*sa3, mean = 0, sd = 1))
}

tmp3 <- rbind(0, tmp3)

sy0 <- rep(0, M)
for (k in 1:M)
{
  F.hat3 <- tmp3[,k]
  sy0[k] <- pull.one(x = runif(1, 0, 1), y = sort(unique(dat$Y)), fhat = F.hat3)
}
m0 <- mean(sy0)
res[j, 1] <- m0 

}

#Save data 
outcome.estimates <- as.data.frame(res)

filename1 <- paste("outcome.estimates-", ITER, ".csv")

write.csv(outcome.estimates, filename1)