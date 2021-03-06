# Copy of R-code given in the Appendix 
# A.2 More on Poisson?HMMs, including estimation by EM

pois.HMM.generate_sample <-                               
 function(n,m,lambda,gamma,delta=NULL)                     
{                                                           
 if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))    
 mvect <- 1:m                                               
 state <- numeric(n)                                        
 state[1] <- sample(mvect,1,prob=delta)                     
 for (i in 2:n)                                             
   state[i]<-sample(mvect,1,prob=gamma[state[i-1],])        
 x <- rpois(n,lambda=lambda[state])                      
 x                                                          
}                                                          


pois.HMM.lalphabeta<-function(x,m,lambda,gamma,delta=NULL)  
{                                                           
 # browser()
 if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))   
 n          <- length(x)                                    
 lalpha     <- lbeta<-matrix(NA,m,n)                       
 allprobs   <- outer(x,lambda,dpois)/100
 browser()
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


pois.HMM.EM <- function(x,m,lambda,gamma,delta,            
                        maxiter=1000,tol=1e-6,...)         
{                                                          
 # browser()
 n              <- length(x)                                        
 lambda.next    <- lambda                                   
 gamma.next     <- gamma                                    
 delta.next     <- delta                                   
 for (iter in 1:maxiter)                                    
   {                                                        
   allprobs    <- outer(x,lambda,dpois)/100
   lallprobs    <- log(allprobs)
   fb  <-  pois.HMM.lalphabeta(x,m,lambda,gamma,delta=delta)   
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
   lambda.next[j] <- sum(exp(la[j,]+lb[j,]-llk)*x)/         
                     sum(exp(la[j,]+lb[j,]-llk))            
   }                                                       
   gamma.next <- gamma.next/apply(gamma.next,1,sum)         
   print(gamma.next)
   delta.next <- exp(la[,1]+lb[,1]-llk)                     
   delta.next <- delta.next/sum(delta.next)                 
   crit       <- sum(abs(lambda-lambda.next)) +             
                 sum(abs(gamma-gamma.next)) +               
                 sum(abs(delta-delta.next))                 
   
   print(paste("Iteration:",iter," Crit:", crit, "LLK:", llk))  
   
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
   print(lambda)
   }                                                        
 print(paste("No convergence after",maxiter,"iterations"))  
 NA                                                         
}                                                           

statdist <- function(gamma){
  m = dim(gamma)[1]
  matrix(1,1,m) %*% solve(diag(1,m) - gamma + matrix(1,m,m))
}



pois.HMM.viterbi<-function(x,m,lambda,gamma,delta=NULL,...) 
{                                                           
 if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))   
 n         <- length(x)                                     
 poisprobs <- outer(x,lambda,dpois)                         
 xi        <- matrix(0,n,m)                                
 foo       <- delta*poisprobs[1,]                          
 xi[1,]    <- foo/sum(foo)                                 
 for (i in 2:n)                                             
   {                                                        
   foo    <- apply(xi[i-1,]*gamma,2,max)*poisprobs[i,]      
   xi[i,] <- foo/sum(foo)                                   
   }                                                        
 iv<-numeric(n)                                             
 iv[n]     <-which.max(xi[n,])                             
 for (i in (n-1):1)                                         
   iv[i] <- which.max(gamma[,iv[i+1]]*xi[i,])               
 iv                                                         
}                                                          


pois.HMM.state_probs <-                                    
 function(x,m,lambda,gamma,delta=NULL,...)                 
{                                                           
 if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m)) 
 n          <- length(x)                                    
 fb         <- pois.HMM.lalphabeta(x,m,lambda,gamma,        
                 delta=delta)                               
 la         <- fb$la                                        
 lb         <- fb$lb                                        
 c          <- max(la[,n])                                  
 llk        <- c+log(sum(exp(la[,n]-c)))                   
 stateprobs <- matrix(NA,ncol=n,nrow=m)                     
 for (i in 1:n) stateprobs[,i]<-exp(la[,i]+lb[,i]-llk)     
 stateprobs                                                
}                                                           


pois.HMM.local_decoding <-                                  
 function(x,m,lambda,gamma,delta=NULL,...)                 
{  
 n   <- length(x)    #inserted by IM                                                      
 stateprobs <-                                              
    pois.HMM.state_probs(x,m,lambda,gamma,delta=delta)      
 ild <- rep(NA,n)                                           
 for (i in 1:n) ild[i]<-which.max(stateprobs[,i])           
 ild                                                        
}                                                           


pois.HMM.state_prediction <-                               
 function(x,m,lambda,gamma,delta=NULL,H=1,...)              
{                                                           
 if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))  
 n          <- length(x)                                    
 fb         <- pois.HMM.lalphabeta(x,m,                     
                 lambda,gamma,delta=delta)                  
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


pois.HMM.forecast <- function(x,m,lambda,gamma, 
       delta=NULL,xrange=NULL,H=1,...)          
{                                               
 if(is.null(delta))                             
       delta<-solve(t(diag(m)-gamma+1),rep(1,m))   
 if(is.null(xrange))                            
       xrange<-qpois(0.001,min(lambda)):        
               qpois(0.999,max(lambda))         
 n        <- length(x)                          
 allprobs <- outer(x,lambda,dpois)              
 allprobs <- ifelse(!is.na(allprobs),allprobs,1)   
 foo      <- delta*allprobs[1,]                 
 sumfoo   <- sum(foo)                           
 lscale   <- log(sumfoo)                        
 foo      <- foo/sumfoo                         
 for (i in 2:n)                                 
   {                                            
   foo    <- foo%*%gamma*allprobs[i,]           
   sumfoo <- sum(foo)                          
   lscale <- lscale+log(sumfoo)                 
   foo    <- foo/sumfoo                         
   }                                            
 xi       <- matrix(NA,nrow=m,ncol=H)           
 for (i in 1:H)                                 
   {                                            
   foo    <- foo%*%gamma                        
   xi[,i] <- foo                                
   }                                            
 allprobs <- outer(xrange,lambda,dpois)         
 fdists   <- allprobs%*%xi[,1:H]                
 list(xrange=xrange,fdists=fdists)              
}                                               


pois.HMM.conditionals <-                                    
 function(x,m,lambda,gamma,delta=NULL,xrange=NULL,...)      
{                                                           
 if(is.null(delta))                                        
   delta  <- solve(t(diag(m)-gamma+1),rep(1,m))             
 if(is.null(xrange))                                        
   xrange <-qpois(0.001,min(lambda)):                       
            qpois(0.999,max(lambda))                        
 n      <- length(x)                                        
 fb     <- pois.HMM.lalphabeta(x,m,lambda,gamma,delta=delta)  
 la     <- fb$la                                            
 lb     <- fb$lb                                           
 la     <- cbind(log(delta),la)                             
 lafact <- apply(la,2,max)                                  
 lbfact <- apply(lb,2,max)                                  
 w      <- matrix(NA,ncol=n,nrow=m)                         
 for (i in 1:n)                                             
   {                                                        
   foo   <- (exp(la[,i]-lafact[i])%*%gamma)*               
             exp(lb[,i]-lbfact[i])                          
   w[,i] <- foo/sum(foo)                                    
   }                                                        
 allprobs <- outer(xrange,lambda,dpois)                     
 cdists   <- allprobs%*%w                                   
 list(xrange=xrange,cdists=cdists)                         
}                                                           


pois.HMM.pseudo_residuals <-                                
 function(x,m,lambda,gamma, delta=NULL,...)                 
{                                                           
 if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))   
 n        <- length(x)                                      
 cdists   <- pois.HMM.conditionals(x,m,lambda, gamma,       
                  delta=delta,xrange=0:max(x))$cdists       
 cumdists <- rbind(rep(0,n),apply(cdists,2,cumsum))         
 ul <- uh <- rep(NA,n)                                      
 for (i in 1:n)                                             
   {                                                        
   ul[i]  <- cumdists[x[i]+1,i]                             
   uh[i]  <- cumdists[x[i]+2,i]                             
   }                                                        
 um       <- 0.5*(ul+uh)                                   
 npsr     <- qnorm(rbind(ul,um,uh))                         
 npsr                                                       
 }                                                           

lambda <- c(10, 20,30, 50, 80, 70, 90, 100, 10)
m <- length(lambda)
#gamma <- matrix(c(0.9,0.1,0.3,0.7), nrow = m, ncol = m, byrow = TRUE)
gamma <- matrix(rep(1,m*m), nrow = m, ncol = m, byrow = TRUE)
gamma = gamma/apply(gamma,1,sum)#create stochastic transition matrix
n <- 100
set.seed(1)
x = pois.HMM.generate_sample(n,m,lambda, gamma)
delta = statdist(gamma)
pois.HMM.EM(x,m,lambda,gamma,delta)
