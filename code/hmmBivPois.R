statdist <- function(gamma){
  m = dim(gamma)[1]
  matrix(1,1,m) %*% solve(diag(1,m) - gamma + matrix(1,m,m))
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
  # browser()
  if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))   
  n          <- dim(x)[1]
  lalpha     <- lbeta<-matrix(NA,m,n)                       
  allprobs   <- PMatAll(x,lambda_buy, lambda_sell)
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

# ---- bi.pois.HMM.EM ----

bi.pois.HMM.EM <- function(x,m_buy,m_sell,lambda_buy,lambda_sell,gamma,delta,            
                           maxiter=1000,tol=1e-6,...)         
{                                                          
  # browser()
  n                <- dim(x)[1]   # num of observations
  m                <- m_buy * m_sell
  lambda_buy.next  <- lambda_buy
  lambda_sell.next <- lambda_sell
  gamma.next       <- gamma                                    
  delta.next       <- delta                                   
  
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
  
  for (iter in 1:maxiter)                                    
  {                                                        
    # if(iter==80) browser()
    lallprobs    <- log(PMatAll(x, lambda_buy, lambda_sell))
    fb  <-  bi.pois.HMM.lalphabeta(x,m,lambda_buy,lambda_sell,gamma,delta=delta)   
    la  <-  fb$la                                            
    lb  <-  fb$lb                                            
    c   <-  max(la[,n])                                      
    llk <- c+log(sum(exp(la[,n]-c)))                         
    # browser()
    for (j in 1:m)                                           
    {                                                       
      for (k in 1:m)                                         
      {                                                      
        gamma.next[j,k] <- gamma[j,k]*sum(exp(la[j,1:(n-1)]+   
                                                lallprobs[2:n,k]+lb[k,2:n]-llk))  
      }                                                      
    }                                                       
    gamma.next <- gamma.next/apply(gamma.next,1,sum)         
    # print(gamma.next)
    
    uhat <- function (j,t) {
      exp(la[j,t]+lb[j,t]-llk)
    }
    
    buy              <- x[,1]
    sell             <- x[,2]
    
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
    
    delta.next <- exp(la[,1]+lb[,1]-llk)                     
    delta.next <- delta.next/sum(delta.next)                 
    
    crit       <- sum(abs(lambda_buy-lambda_buy.next)) +             
      sum(abs(lambda_sell-lambda_sell.next)) +
      sum(abs(gamma-gamma.next)) +               
      sum(abs(delta-delta.next))                 
    
    print(paste("Iteration:",iter," Crit:", crit, "LLK:", llk))  
    
    if(is.na(crit))                                             
    {                                                      
      AIC    <- NA
      BIC    <- NA
      return(list(lambda_buy=lambda_buy,lambda_sell=lambda_sell,gamma=gamma,delta=delta,    
                  mllk=-llk,AIC=AIC,BIC=BIC))                     
    }                                                      
    else if(crit<tol)                                             
    {                                                      
      np     <- m*m+m-1                                      
      AIC    <- -2*(llk-np)                                  
      BIC    <- -2*llk+np*log(n)                             
      return(list(lambda_buy=lambda_buy,lambda_sell=lambda_sell,gamma=gamma,delta=delta,    
                  mllk=-llk,AIC=AIC,BIC=BIC))                     
    }                                                      
    lambda_buy  <- lambda_buy.next                               
    lambda_sell <- lambda_sell.next                               
    gamma      <- gamma.next                                 
    delta      <- delta.next                                
  }                                                        
  print(paste("No convergence after",maxiter,"iterations"))  
  NA                                                         
}                                                           

# ---- state decoding ----

bi.pois.HMM.state_probs <-                                    
  function(x,m,lambda_buy,lambda_sell,gamma,delta=NULL,...)                 
  {                                                           
    if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m)) 
    n          <- dim(x)[1]                                    
    fb         <- bi.pois.HMM.lalphabeta(x,m,lambda_buy,lambda_sell,gamma,        
                                      delta=delta)                               
    la         <- fb$la                                        
    lb         <- fb$lb                                        
    c          <- max(la[,n])                                  
    llk        <- c+log(sum(exp(la[,n]-c)))                   
    stateprobs <- matrix(NA,ncol=n,nrow=m)                     
    for (i in 1:n) stateprobs[,i]<-exp(la[,i]+lb[,i]-llk)     
    stateprobs                                                
  }                                                           

bi.pois.HMM.local_decoding <- function(x,m,lambda_buy,lambda_sell,gamma,delta=NULL,...)                 
  {  
    n   <- dim(x)[1]
    stateprobs <- bi.pois.HMM.state_probs(x,m,lambda_buy,lambda_sell,gamma,delta=delta)      
    ild <- rep(NA,n)                                           
    for (i in 1:n) ild[i]<-which.max(stateprobs[,i])           
    ild                                                        
  }                                                           

# ---- prediction ----

bi.pois.HMM.state_prediction <- function(x,m,lambda_buy,lambda_sell,gamma,delta=NULL,H=1,...)
  {                                                           
    if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))  
    n          <- dim(x)[1]
    fb         <- bi.pois.HMM.lalphabeta(x,m,lambda_buy,lambda_sell,gamma,delta=delta)                  
    la         <- fb$la                                       
    c          <- max(la[,n])                                 
    llk        <- c+log(sum(exp(la[,n]-c)))                    
    statepreds <- matrix(NA,ncol=H,nrow=m)                     
    foo1       <- exp(la[,n]-llk)                              
    foo2       <- diag(m)                                      
    for (i in 1:H)                                             
    {                                                        
      foo2           <- foo2%*%gamma                           
      statepreds[,i] <- foo1%*%foo2                            
    }                                                        
    statepreds                                                 
  }                                                           


# ---- test ----

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

lambda_buy = c(20,30)
lambda_sell = c(5,20)

m_buy = length(lambda_buy)
m_sell = length(lambda_sell)

mn <- m_buy * m_sell

# gamma <- matrix(runif(mn), nrow = mn, ncol = mn, byrow = TRUE)
gamma <- matrix(rep(1,mn), nrow = mn, ncol = mn, byrow = TRUE)
gamma = gamma/apply(gamma,1,sum)#create stochastic transition matrix

# Generate synthetic data
# set.seed(1)
# n <- 50
# x = bi.pois.HMM.generate_sample(n, mn, lambda_buy,lambda_sell, gamma)
delta = matrix(rep(1, mn), ncol = mn)
delta = delta/apply(delta, 1, sum)
# set.seed(1)
n <- 10
# x = bi.pois.HMM.generate_sample(n, mn, lambda_buy,lambda_sell, gamma)
# delta = c(0.3, 0.3, 0.1, 0.1, 0.1, 0.1)
# print(bi.pois.HMM.EM(x,m_buy,m_sell,c(10,11,12), c(13,14),gamma,delta))
# model <- bi.pois.HMM.EM(x,m_buy,m_sell,lambda_buy, lambda_sell,gamma,delta)
# print(model)
# states = bi.pois.HMM.local_decoding(x, mn, model$lambda_buy, model$lambda_sell,model$gamma)

