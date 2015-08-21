statdist <- function(gamma){
  m = dim(gamma)[1]
  matrix(1,1,m) %*% solve(diag(1,m) - gamma + matrix(1,m,m))
}

bi.pois.HMM.generate_sample <- function(n,m,lambda_buy,lambda_sell,gamma,delta=NULL)
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

PMatAll <- function(x, lambda_buy, lambda_sell) {
  buy <- x[,1]
  sell <- x[,2]
  m_buy = length(lambda_buy)
  m_sell = length(lambda_sell)
  n = dim(x)[1]
  pmat = matrix(data=0, nrow = n, ncol = m_buy*m_sell)
  counter = 1
  for(i in 1:m_buy)
  {
    for(j in 1:m_sell)
    {
      pmat[,counter]  = dpois(buy, lambda_buy[i]) * dpois(sell, lambda_sell[j])
      counter = counter + 1
    }
  }
  pmat
}

##
# generate LOGof forward(alpha) and backward(beta) probabilities  p60
#         alpha(mxn) is JOINT prob of being in state(1:m) at time(1:n), AND known observations x up to that time
#         beta (mxn) is CONDL prob of being in state(1:m) at time(1:n), conditional on observations x AFTER that time
##

# ---- bi.pois.HMM.lalphabeta ----
bi.pois.HMM.lalphabeta<-function(x,m,lambda_buy,lambda_sell,gamma,delta=NULL)  
{                                                           
  #browser()
  if(is.null(delta)){
    delta = solve(t(diag(m)-gamma+1),rep(1,m))   
  }
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

# ---- bi.pois.HMM.EM ----
bi.pois.HMM.EM <- function(x,m_buy,m_sell,lambda_buy,lambda_sell,gamma,delta,            
                        maxiter=1000,tol=1e-6,...)         
{
  #CDDL helper functions
  vhat <- function (i,j) {
    gamma[i,j] * sum(exp(la[i,1:(t-1)] + lallprobs[2:t,j] + lb[j,2:t] - llk))
  }
  
  uhat <- function (j,t) {
    exp(la[j,t]+lb[j,t]-llk)
  }
  
  t                <- dim(x)[1]   # num of observations
  m                <- m_buy * m_sell
  buy              <- x[,1]
  sell             <- x[,2]
  lambda_buy.next  <- lambda_buy
  lambda_sell.next <- lambda_sell
  gamma.next       <- gamma
  delta.next       <- delta
  for (iter in 1:maxiter)                                    
  {                                                        
    lallprobs    <- log(PMatAll(x, lambda_buy, lambda_sell))
    fb  <-  bi.pois.HMM.lalphabeta(x,m,lambda_buy,lambda_sell,gamma,delta)   
    la  <-  fb$la                                            
    lb  <-  fb$lb                                            
    c   <-  max(la[,t])                                      
    llk <- c+log(sum(exp(la[,t]-c)))
       
    #initiliase a state index lookup table for resolving (i,j) -> stateIndex
    stateEnv<-new.env() 
    counter = 1
    for (i in 1:m_buy)               
    {
      for (j in 1:m_sell)
      {
        stateEnv[[paste(i,j)]] = counter
        counter = counter + 1
      }
    }
    
    #calculate gamma
    for (ij in 1:m)     
    {
      #because gamma is mn * mn, ij represents (i,j) the combination of buy and sell states  
      for (kl in 1:m)  
      {
        #similarly kl here represents (k,l) 
        numerator <- vhat(ij,kl)

        denominator = 0
        for (k in 1:m_buy)     
        {
          for (l in 1:m_sell)  
          {
            klprime = stateEnv[[paste(k,l)]]
            denominator = denominator + vhat(ij, klprime)
          }
        }
        gamma.next[ij,kl] <- numerator / denominator
      }
    }
     gamma.next <- gamma.next/apply(gamma.next,1,sum)#"stochastisize"
    print(gamma.next)
    if(sum(gamma.next == 0) > 0)
    {
      #browser()
    }
    
    #calculate lambda_buy_i
    for (i in 1:m_buy)               
    {
      numerator <- 0
      denominator <- 0
      for (j in 1:m_sell)
      {
        ij = stateEnv[[paste(i,j)]]
        numerator <- numerator + sum(uhat(ij)*buy)
        denominator <- denominator + sum(uhat(ij))
      }
      
      lambda_buy.next[i] <- numerator / denominator            
    }
    print(lambda_buy)
    
    #calculate lambda_sell_j
    for (j in 1:m_sell)                                         
    {
      numerator <- 0
      denominator <- 0
      for (i in 1:m_buy)               
      {
        ij = stateEnv[[paste(i,j)]]
        numerator <- numerator + sum(uhat(ij)*sell)
        denominator <- denominator + sum(uhat(ij))
      }
      
      lambda_sell.next[j] <- numerator / denominator
    }
    print(lambda_sell)
    
    #calculate delta
    delta.next <- uhat(t = 1)
    delta.next <- delta.next/sum(delta.next)                 
    
    crit       <- sum(abs(lambda_buy-lambda_buy.next)) +             
                  sum(abs(lambda_sell-lambda_sell.next)) +
                  sum(abs(gamma-gamma.next)) +               
                  sum(abs(delta-delta.next))                 

    if(is.na(crit))
    {
      browser()
    }
    
    print(paste("Iteration:",iter," Crit:", crit, "LLK:", llk))  
    if(crit<tol)                                             
    {                                                      
      np     <- m*m+m-1                                      
      AIC    <- -2*(llk-np)                                  
      BIC    <- -2*llk+np*log(t)                             
      return(list(lambda_buy=lambda_buy,lambda_sell=lambda_sell,gamma=gamma,delta=delta,    
                  mllk=-llk,AIC=AIC,BIC=BIC))                     
    }                                                      
    lambda_buy  <- lambda_buy.next                               
    lambda_sell <- lambda_sell.next                               
    gamma       <- gamma.next                                 
    delta       <- delta.next                                
  }                                                        
  print(paste("No convergence after",maxiter,"iterations"))  
  NA                                                         
}          

# ---- test ----

lambda_buy = c(1,20,30)
lambda_sell = c(10,20)

m_buy = length(lambda_buy)
m_sell = length(lambda_sell)

mn <- m_buy * m_sell

# gamma <- matrix(runif(mn), nrow = mn, ncol = mn, byrow = TRUE)
gamma <- matrix(rep(1,mn), nrow = mn, ncol = mn, byrow = TRUE)
gamma = gamma/apply(gamma,1,sum)#create stochastic transition matrix

# Generate synthetic data
# set.seed(1)
# n <- 10
# x = bi.pois.HMM.generate_sample(n, mn, lambda_buy,lambda_sell, gamma)
# delta = c(0.3, 0.3, 0.1, 0.1, 0.1, 0.1)
# print(bi.pois.HMM.EM(x,m_buy,m_sell,c(10,11,12), c(13,14),gamma,delta))


