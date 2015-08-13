norm.HMM.generate_sample <- function(n,m,mu,sigma,gamma,delta=NULL)
{
  browser()
  if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))
  mvect <- 1:m
  state <- numeric(n)
  state[1] <- sample(mvect,1,prob=delta)
  for (i in 2:n)
    state[i]<-sample(mvect,1,prob=gamma[state[i-1],])
  x <- rnorm(n,mu[state],sigma[state])
  x
}

statdist <- function(gamma){
  m = dim(gamma)[1]
  matrix(1,1,m) %*% solve(diag(1,m) - gamma + matrix(1,m,m))
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

PMat <- function(t, lambda_buy, lambda_sell) {
  buy_t <- t[1]
  sell_t <- t[2]
  m = length(lambda_buy)
  n = length(lambda_sell)
  pmat = matrix(data=0, nrow = m*n, ncol = m*n)
  counter = 1
  for(i in 1:m)
  {
    for(j in 1:n)
    {
      pmat[counter,counter]  = dpois(buy_t, lambda_buy[i]) * dpois(sell_t, lambda_sell[j])
      counter = counter + 1
    }
  }
  pmat
}

##
# generate LOGof forward(alpha) and backward(beta) probabilities  p60
#         alpha(mxn) is JOINT prob of being in state(1:m) at time(1:n), AND known observations x up to that time
#         beta (mxn) is CONDL prob of being in state(1:m) at time(1:n), conditional on observations x AFTER that time
#
#  ... generic function: the only line changed from the pois.HMM.lalphabeta routine is:  allprobs   <- outer(x,lambda,dpois)
##
pois.HMM.lalphabeta<-function(x,m,lambda_buy,lambda_sell,gamma,delta=NULL)  
{                                                           
  #browser()
  #if(is.null(delta)){
    delta = solve(t(diag(m)-gamma+1),rep(1,m))   
  #}
  n          <- dim(x)[1]
  lalpha     <- lbeta <- matrix(NA,m,n)                       
  # allprobs <- outer(x,lambda,dpois)                     # <- OLD
  #allprobs   <- drop(outer(x,mu,dnorm,sd=sigma %x% rep(1,n)))  # nxm array of probabilities for the n observations under the m distributions
  # alpha
  foo        <- delta %*% PMat(x[1,], lambda_buy, lambda_sell)
  # check that a delta argument does not encounter zero-prob states, which would cause all foo to be zero
  if(sum(foo) == 0) foo <- delta
  sumfoo     <- sum(foo)                                    
  lscale     <- log(sumfoo)                                 
  foo        <- foo/sumfoo                                   
  lalpha[,1] <- log(foo)+lscale

  for (i in 2:n)                                             
  {                                                        
    foo        <- foo %*% gamma %*% PMat(x[i,], lambda_buy, lambda_sell)    
    sumfoo     <- sum(foo)                                   
    lscale     <- lscale+log(sumfoo)                         
    foo        <- foo/sumfoo                                 
    lalpha[,i] <- log(foo)+lscale                            # since alpha = foo * scale, in logs this is log(alpha) = log(foo) + log(scale)
  }
  
  # beta 													   
  lbeta[,n]  <- rep(0,m)                                     # beta(T) is the m-vector 1 so log(1)=0 so initialise at the zero-vector
  foo        <- rep(1/m,m)                                   # rescaling beta by its sum=m
  lscale     <- log(m)                                       # p46 
  for (i in (n-1):1)                                         
  {                                                        
    foo        <- gamma %*% (PMat(x[i,], lambda_buy, lambda_sell) %*% foo)           # p61    
    lbeta[,i]  <- log(foo)+lscale                            
    sumfoo     <- sum(foo)                                   
    foo        <- foo/sumfoo                                 
    lscale     <- lscale+log(sumfoo)                         
  }                                                        
  list(la=lalpha,lb=lbeta)                                   
}

pois.HMM.lalphabeta.uni<-function(x,m,lambda,gamma,delta=NULL)  
{                                                           
  if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))   
  n          <- length(x)                                    
  lalpha     <- lbeta<-matrix(NA,m,n)                       
  allprobs   <- outer(x,lambda,dpois)                        
  foo        <- delta*allprobs[1,]                           
  sumfoo     <- sum(foo)                                    
  lscale     <- log(sumfoo)                                 
  foo        <- foo/sumfoo                                   
  lalpha[,1] <- log(foo)+lscale                              
  for (i in 2:n)                                             
  {                                                        
    foo        <- foo%*%gamma*allprobs[i,]                   
    sumfoo     <- sum(foo)                                   
    lscale     <- lscale+log(sumfoo)                         
    foo        <- foo/sumfoo                                 
    lalpha[,i] <- log(foo)+lscale                            
  }                                                        
  lbeta[,n]  <- rep(0,m)                                     
  foo        <- rep(1/m,m)                                   
  lscale     <- log(m)                                       
  for (i in (n-1):1)                                         
  {                                                        
    foo        <- gamma%*%(allprobs[i+1,]*foo)               
    lbeta[,i]  <- log(foo)+lscale                            
    sumfoo     <- sum(foo)                                   
    foo        <- foo/sumfoo                                 
    lscale     <- lscale+log(sumfoo)                         
  }                                                        
  list(la=lalpha,lb=lbeta)                                   
}  

PMat.uni <- function(x, lambda) diag(dpois(x, lambda))

pois.HMM.EM <- function(x,m,lambda_buy,lambda_sell,gamma,delta,            
                        maxiter=1000,tol=1e-6,...)         
{
  n                <- length(x)#dim(x)[1]
  lambda_buy.next  <- lambda_buy
  lambda_sell.next <- lambda_sell
  gamma.next       <- gamma
  delta.next       <- delta
  for (iter in 1:maxiter)                                    
  {                                                        
    lallprobs    <- outer(x,lambda_buy,dpois,log=TRUE)           
    fb  <-  pois.HMM.lalphabeta.uni(x,m,lambda_buy,gamma)#,lambda_sell, delta=delta)   
    la  <-  fb$la                                            
    lb  <-  fb$lb                                            
    c   <-  max(la[,n])                                      
    llk <- c+log(sum(exp(la[,n]-c)))
    browser()
    for (j in 1:m)                                           
    {                                                       
      for (k in 1:m)                                         
      {                                                      
        orig <- gamma[j,k] * sum( exp(la[j,1:(n-1)] + lallprobs[2:n,k] + lb[k,2:n] - llk) )
        new <- gamma[j,k] * sum( exp(la[j,1:(n-1)] + PMat.uni(x, lambda_buy) + lb[k,2:n] - llk) )
        gamma.next[j,k] <- orig
      }                                                      
      lambda.next[j] <- sum(exp(la[j,]+lb[j,]-llk)*x)/         
        sum(exp(la[j,]+lb[j,]-llk))            
    }                                                       
    gamma.next <- gamma.next/apply(gamma.next,1,sum)         
    delta.next <- exp(la[,1]+lb[,1]-llk)                     
    delta.next <- delta.next/sum(delta.next)                 
    crit       <- sum(abs(lambda-lambda.next)) +             
      sum(abs(gamma-gamma.next)) +               
      sum(abs(delta-delta.next))                 
    if(crit<tol)                                             
    {                                                      
      np     <- m*m+m-1                                      
      AIC    <- -2*(llk-np)                                  
      BIC    <- -2*llk+np*log(n)                             
      return(list(lambda=lambda,gamma=gamma,delta=delta,    
                  mllk=-llk,AIC=AIC,BIC=BIC))                     
    }                                                      
    lambda     <- lambda.next                               
    gamma      <- gamma.next                                 
    delta      <- delta.next                                
  }                                                        
  print(paste("No convergence after",maxiter,"iterations"))  
  NA                                                         
}                                                           


lambda_buy = c(3, 15)
lambda_sell = c(10, 40)
delta_buy = c(0.3, 0.6, 0.1)
delta_sell = c(0.1, 0.9)

m = length(lambda_buy)
n = length(lambda_sell)

mn <- m*n
gamma = matrix(data = rep(1/mn, mn*mn), nrow = mn, ncol = mn)

# Generate synthetic data
x = pois.HMM.generate_sample(10, mn, lambda_buy, lambda_sell, gamma)
# Fit normal-HMM with the EM algorith,
#res <- norm.HMM.EM(x,m,mu,sigma,gamma,delta)
res <- pois.HMM.lalphabeta(x,mn,lambda_buy,lambda_sell,gamma,delta)

#test EM with unvariate sample first
gamma.uni = matrix(data = rep(1/m, m*m), nrow = m, ncol = m)
pois.HMM.EM(x[,1],m,lambda_buy,lambda_sell,gamma.uni,delta_buy)


