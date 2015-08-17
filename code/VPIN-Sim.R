# ---- vpin.sim ----

finv <- function(u,p)
{
#  return (qnorm(u, mean = p))
#  return (qpois(u, p))
  return (rpois(1, p))
}

numSim = 1000
Vpin = numeric(numSim)

V = 50 #size of each volume bucket
n = 100 #number of volume buckets

alpha = 0.5 #probability of information event
delta = 0.4 #probability of bad news event
mu = 40 #arrival rate of informed trader
epsilon = (V - alpha * mu)/2 #arrival rate of uninformed trader

s = 1
while(s < numSim){

  Vbuy = numeric(n) #buy volume buckets
  Vsell = numeric(n) #sell volume buckets
  
  j = 1
 
 while(j <= n){
    
    u1 = runif(1)
    u2 = runif(1)
    u3 = runif(1)
    
    if(u1 < alpha){ #we have an information event 
      if(u2 < delta){ #its a bad news event
        
        #only uninformed traders buy when there's bad news
        Vbuy[j] = finv(u3, epsilon)  
        
        #both informed and uninformed traders sell when there's bad news
        Vsell[j] = finv(u3, mu + epsilon)
      }
      else{ #its a good news event
        
        #both informed and uninformed traders buy when there's good news
        Vbuy[j] = finv(u3, mu + epsilon)
        
        #only uninformed traders sell when there's good news
        Vsell[j] = finv(u3, epsilon)
      }
    }
    else{
      #no information event
      
      #uninformed traders buy and sell in equal quantities 
      Vbuy[j] = finv(u3, epsilon)
      Vsell[j] = Vbuy[j]
    }
    
    j = j + 1
  }
  
  Vpin[s] = sum(abs(Vbuy - Vsell)) / sum(Vbuy + Vsell)
  
  s = s + 1
}
meanVpin = sum(Vpin)/numSim
varVpin =  sum(Vpin^2)/(numSim - 1) - (sum(Vpin)/(numSim * (numSim - 1)))^2

source('jump.R')

testdata <- matrix(c(Vbuy, Vsell),byrow=T,ncol=2)
temp <- jump(testdata,y=c(1.5,2,2.5),rand=10,trace=F,plotjumps =FALSE)

source('hmmBivPois.R')
lambda_buy = c(10,55)
lambda_sell = c(5,42)

m_buy = length(lambda_buy)
m_sell = length(lambda_sell)
mn <- m_buy * m_sell

delta_buy = matrix(rep(1, mn), ncol = mn)
delta_buy = delta_buy/apply(delta_buy, 1, sum)

gamma <- matrix(rep(1,mn), nrow = mn, ncol = mn, byrow = TRUE)
gamma = gamma/apply(gamma,1,sum)#create stochastic transition matrix

result <- bi.pois.HMM.EM(testdata,m_buy,m_sell,lambda_buy,lambda_sell,gamma,delta_buy)
print(result)

