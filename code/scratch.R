#univariate example:
library(depmixS4)

lambdas = c(3, 15)
deltas = c(0.3, 0.7)
cdf = cumsum(deltas)
unifvec = runif(1000)
d1 = rpois(sum(unifvec < cdf[1]),lambdas[1])
d2 = rpois(sum(unifvec >= cdf[1]),lambdas[2])
x = c(d1,d2) 
x.df = as.data.frame(x)

mx <- depmix(x~1,nstates =2,family=poisson(),data=x.df)
set.seed(1)
fmx <- fit(mx)
exp(summary(fmx))


# this creates data with a single change point with Poisson data
set.seed(3)
y1 <- rpois(50,10)
y2 <- rpois(50,30)
ydf <- data.frame(y=c(y1,y2))

# fit models with 1 to 3 states
m1 <- depmix(y~1,ns=1,family=poisson(),data=ydf)
set.seed(1)
fm1 <- fit(m1)
m2 <- depmix(y~1,ns=2,family=poisson(),data=ydf)
set.seed(1)
fm2 <- fit(m2)
m3 <- depmix(y~1,ns=3,family=poisson(),data=ydf)
set.seed(1)
fm3 <- fit(m3,em=em.control(maxit=500))

# plot the BICs to select the proper model
plot(1:3,c(BIC(fm1),BIC(fm2),BIC(fm3)),ty="b")

exp(summary(fm2))

y <- rpois(1000, lambda = 15)
respst <- c(0,1,2,1)
trst <- c(0.9,0.1,0.1,0.9)
df <- data.frame(y=y)
mod <- depmix(y~1,data=df,respst=respst,trst=trst,inst=c(0.5,0.5),nti=1000,nst=2)
mod <- depmix(y~1,data=df,nti=1000,nst=1, family = poisson())
mod <- simulate(mod)


#bivariate example

lambda_buy = c(3, 15)
lambda_sell = c(10, 40)
delta_buy = c(0.3, 0.6, 0.1)
delta_sell = c(0.1, 0.9)

simx <- function(delta, lambdas)
{
  cdf = cumsum(delta)
  unifvec = runif(1000)
  d1 = rpois(sum(unifvec < cdf[1]),lambdas[1])
  d2 = rpois(sum(unifvec >= cdf[1]),lambdas[2])
  x = c(d1,d2) 
  return (x)
}

pois.HMM.generate_sample <- function(n,m,lambda_buy,lambda_sell,gamma,delta=NULL)
{
  x = matrix(nrow = n, ncol = 2)
  if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))
  lambda_combo = expand.grid(lambda_buy, lambda_sell)#generates all the pairwise combinations of lambdas
  mvect <- 1:m
  state <- numeric(n)
  state[1] <- sample(mvect,1,prob=delta)
  for (i in 2:n)
    state[i]<-sample(mvect,1,prob=gamma[state[i-1],])
  for (i in 1:n)
  {
    x[i,] = rpois(2, t(as.vector(lambda_combo)[state[i],]))
  }
  x
}
m = length(lambda_buy)
n = length(lambda_sell)

mn <- m*n
gamma = matrix(data = rep(1/mn, mn*mn), nrow = mn, ncol = mn)

buy = pois.HMM.generate_sample(10, mn, lambda_buy, lambda_sell, gamma)
buy = simx(lambda_buy, delta_buy)
sell = simx(lambda_sell, delta_sell)
x = cbind(buy, sell)
bi.x.df = data.frame(x)

plot.ts(x, xlab = NULL, ylab = NULL)

bm1 <- depmix(x~1,ns=2,family=poisson(),data=bi.x.df)
set.seed(1)
bfm1 <- fit(bm1)
exp(summary(bfm1))


PMat2 <- function(t, lambda_buy, lambda_sell) {
  buy_t <- t[1]
  sell_t <- t[2]
  m = length(lambda_buy)
  n = length(lambda_sell)
  ret = matrix(data=0, nrow = m*n, ncol = m*n)
  counter = 1
  for(i in 1:m)
  {
    for(j in 1:n)
    {
      ret[counter,counter]  = dpois(buy_t, lambda_buy[i]) * dpois(sell_t, lambda_sell[j])
      counter = counter + 1
    }
  }
  return(ret)
}

PMat <- function(x, lambda) diag(dpois(x, lambda))

bi.pois.hmm.alpha <- function(x, gamma, delta = NULL)
{
  n = length(x[,1])
  alpha = vector(mode="list", length=n)
  alpha0= statdist(gamma)
  
  for(i in 1:n)
  {
    if(i==1)
    {
      prevAlpha = alpha0
    }
    else
    {
      prevAlpha = alpha[[i-1]]
    }
    browser()
    alpha[[i]] = prevAlpha %*% gamma %*% PMat2(x[i,],lambda_buy, lambda_sell)
  }
  return(alpha[[n]])
}


alpha = bi.pois.hmm.alpha(x, gamma)
PMat2(x[1,], lambda_buy, lambda_sell)
