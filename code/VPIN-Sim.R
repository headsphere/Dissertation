numSim = 100000
Vpin = numeric(numSim)

V = 100 #size of volume buckets
Vb = numeric(V) #buy volume buckets
Vs = numeric(V) #sell volume buckets
n = 50 #number of buckets

alpha = 0.3 #probability of information event
delta = 0.5 #probability of bad news event
mu = 10 #arrival rate of informed trader
epsilon = (V - alpha * mu)/2 #arrival rate of uninformed trader

s = 1

while(s < numSim){
  
  j = 1
  
  while(j < n){
    
    u1 = runif(1)
    u2 = runif(1)
    u3 = runif(1)
    
    if(u1 < alpha){
      #we have an information event 
      
      if(u2 < delta)
      {
        #its a bad news event
        Vb[j] = finv(u3, epsilon)
        Vs[j] = finv(u3, mu + epsilon)
      }
      else{
        #its a good news event
        Vb[j] = finv(u3, mu + epsilon)
        Vs[j] = finv(u3, epsilon)
      }
    }
    else{
      #no information event
      Vb[j] = finv(u3, epsilon)
      Vs[j] = Vb[j]
    }
    
    j = j + i
  }
  
  num = 0
  for(j in 1:n)
  {
    num = num + abs(Vb[j] - Vs[j])
  }
  
  denom = 0
  for(j in 1:n)
  {
    denom = denom + (Vb[j] + Vs[j])
  }
  Vpin[s] = num/denom
}
meanVpin = sum(Vpin)/numSim
varVpin =  sum(Vpin^2)/(numSim - 1) - (sum(Vpin)/(numSim * (numSim - 1)))^2

finv <- function(u,p)
{
  #TODO - WHAT IS THIS?!?!  
}