norm.HMM.generate_sample <- function(n,m,mu,sigma,gamma,delta=NULL)
{
  if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))
  mvect <- 1:m
  state <- numeric(n)
  state[1] <- sample(mvect,1,prob=delta)
  for (i in 2:n)
    state[i]<-sample(mvect,1,prob=gamma[state[i-1],])
  x <- rnorm(n,mu[state],sigma[state])
  x
}

##
# generate LOGof forward(alpha) and backward(beta) probabilities  p60
#         alpha(mxn) is JOINT prob of being in state(1:m) at time(1:n), AND known observations x up to that time
#         beta (mxn) is CONDL prob of being in state(1:m) at time(1:n), conditional on observations x AFTER that time
#
#  ... generic function: the only line changed from the pois.HMM.lalphabeta routine is:  allprobs   <- outer(x,lambda,dpois)
##
norm.HMM.lalphabeta<-function(x,m,mu,sigma,gamma,delta=NULL)  
{                                                           
  browser()
  if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))   
  n          <- length(x)                                    
  lalpha     <- lbeta <- matrix(NA,m,n)                       
  # allprobs <- outer(x,lambda,dpois)                     # <- OLD
  allprobs   <- drop(outer(x,mu,dnorm,sd=sigma %x% rep(1,n)))  # nxm array of probabilities for the n observations under the m distributions
  
  # alpha
  foo        <- delta*allprobs[1,]                           
  # check that a delta argument does not encounter zero-prob states, which would cause all foo to be zero
  if(sum(foo) == 0) foo <- delta
  sumfoo     <- sum(foo)                                    
  lscale     <- log(sumfoo)                                 
  foo        <- foo/sumfoo                                   
  lalpha[,1] <- log(foo)+lscale
  #browser()
  for (i in 2:n)                                             
  {                                                        
    foo        <- foo %*% gamma * allprobs[i,]                   
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
    foo        <- gamma %*% (allprobs[i+1,] * foo)           # p61    
    lbeta[,i]  <- log(foo)+lscale                            
    sumfoo     <- sum(foo)                                   
    foo        <- foo/sumfoo                                 
    lscale     <- lscale+log(sumfoo)                         
  }                                                        
  list(la=lalpha,lb=lbeta)                                   
}

# case m=2
m= 2
mu = c(-3,3)
sigma = c(1,2)
gamma= rbind(c(0.9,0.1),c(0.1,0.9))
delta= statdist(gamma)
n=200
# Generate synthetic data
x <- norm.HMM.generate_sample(n,m,mu,sigma,gamma,delta)
# Fit normal-HMM with the EM algorith,
#res <- norm.HMM.EM(x,m,mu,sigma,gamma,delta)
res <- norm.HMM.lalphabeta(x,m,mu,sigma,gamma,delta)
