# ---- vpin.sim ----

generate.trades.sim <- function (n, alpha, delta, epsilon, mu) {
  # numSim = 1000
  # Vpin = numeric(numSim)
  # s = 1
  # while(s < numSim){
  # 
    Vbuy = numeric(n) #buy volume buckets
    Vsell = numeric(n) #sell volume buckets
    j = 1
   
   while(j <= n){
      
      u1 = runif(1)
      u2 = runif(1)
      
      if(u1 < alpha){ #we have an information event 
        if(u2 < delta){ #its a bad news event
          
          #only uninformed traders buy when there's bad news
          Vbuy[j] = rpois(1, epsilon)  
          
          #both informed and uninformed traders sell when there's bad news
          Vsell[j] = rpois(1, mu + epsilon)
        }
        else{ #its a good news event
          
          #both informed and uninformed traders buy when there's good news
          Vbuy[j] = rpois(1, mu + epsilon)
          
          #only uninformed traders sell when there's good news
          Vsell[j] = rpois(1, epsilon)
        }
      }
      else{
        #no information event
        
        #uninformed traders buy and sell in equal quantities 
        Vbuy[j] = rpois(1, epsilon)
        Vsell[j] = Vbuy[j]
      }
      
      j = j + 1
    }
    
  #   Vpin[s] = sum(abs(Vbuy - Vsell)) / sum(Vbuy + Vsell)
  #   
  #   s = s + 1
  # }
  # meanVpin = sum(Vpin)/numSim
  # varVpin =  sum(Vpin^2)/(numSim - 1) - (sum(Vpin)/(numSim * (numSim - 1)))^2
    
    return(data.frame(Buckets = 1:n, Buy=Vbuy, Sell=Vsell))
}

# ---- vpin.sim.test ----

V = 55 #size of each volume bucket
n = 10000 #number of volume buckets

alpha = 0.28 #probability of information event
delta = 0.33 #probability of bad news event
mu = 80 #arrival rate of informed trader
epsilon = (V - alpha * mu)/2 #arrival rate of uninformed trader

trades = generate.trades.sim(n, alpha, delta, epsilon, mu)
Vbuy = trades$Vbuy
Vsell = trades$Vsell

source('jump.R')

testdata <- matrix(c(Vbuy, Vsell),byrow=T,ncol=2)
temp <- jump(testdata,y=c(1.5,2,2.5),rand=10,trace=F,plotjumps =FALSE)

colnames(testdata) <- c("buy", "sell")
(cl <- kmeans(testdata, 9))
plot(testdata, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex = 2)

source('hmmBivPois.R')
lambda_sell = c(mu + epsilon, epsilon,epsilon)
lambda_buy = c(epsilon, mu + epsilon, epsilon)
# lambda_buy = c(10,55)
# lambda_sell = c(5,42)

m_buy = length(lambda_buy)
m_sell = length(lambda_sell)
mn <- m_buy * m_sell

delta_vec = c(alpha*delta, alpha*(1-delta), 1-alpha)
delta_vec = as.vector(delta_vec %o% delta_vec)
delta_buy = matrix(rep(1, mn), ncol = mn)
delta_buy = delta_buy/apply(delta_buy, 1, sum)

m_buy = length(lambda_buy)
m_sell = length(lambda_sell)
mn <- m_buy * m_sell

delta_buy = matrix(rep(1, mn), ncol = mn)
delta_buy = delta_buy/apply(delta_buy, 1, sum)

gamma <- matrix(rep(1,mn), nrow = mn, ncol = mn, byrow = TRUE)
gamma = gamma/apply(gamma,1,sum)

model <- bi.pois.HMM.EM(testdata,m_buy,m_sell,lambda_buy,lambda_sell,gamma,delta_buy)
print(model)

states = bi.pois.HMM.local_decoding(testdata, mn, model$lambda_buy, model$lambda_sell,gamma)
stateProbs = bi.pois.HMM.state_probs(testdata, mn, model$lambda_buy, model$lambda_sell,gamma)


