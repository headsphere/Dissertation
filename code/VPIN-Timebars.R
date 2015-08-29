source('VPIN-Sim.R')

V = 55 #size of each volume bucket
n = 100 #number of time bars
alpha = 0.28 #probability of information event
delta = 0.33 #probability of bad news event
mu = 80 #arrival rate of informed trader
epsilon = (V - alpha * mu)/2 #arrival rate of uninformed trader

trades = generate.trades.sim(n, alpha, delta, mu, epsilon) 

# trades = data.frame(Buckets=seq(1:10), 
#                     Buy=c(1,2,2,3,11,3,5,2,2,1), 
#                     Sell=c(2,4,5,7,3,6,2,4,4,1), 
#                     States=c(1,2,3,2,3,3,2,1,2,1))

trades=data.frame(trades, Total=trades$Buy+trades$Sell)
expanded = data.frame()
for(i in 1:length(trades[,1]))
{
  repeats = trades[i,]$Total
  state = trades[i,]$States
  expanded = rbind(expanded, data.frame(Timebar=seq(1:repeats)*0+i, TimeState=state, Volume=repeats))
}

n <- length(expanded[,1])
noOfVolBars <- floor(n / V)
VolBars = data.frame(VolBar = rep(NA,noOfVolBars), 
                     State = rep(NA,noOfVolBars), 
                     Buy = rep(NA,noOfVolBars), 
                     Sell = rep(NA,noOfVolBars), 
                     Imbalance = rep(NA,noOfVolBars), 
                     VPIN = rep(NA,noOfVolBars), 
                     stringsAsFactors = FALSE)
for(i in 1:noOfVolBars)
{
  end = i*V
  start = end-V+1
  volData = expanded[start:end,]
  
  mode = names(sort(-table(volData$TimeState)))[1]

  timeBars = unique(volData$Timebar)
  timeBarWeights = numeric()
  j = 1
  for(t in min(timeBars):max(timeBars))
  {
    timeBarWeights[j] = sum(volData$Timebar == t)/V
    j = j + 1
  }
  buy = trades[timeBars,]$Buy %*% timeBarWeights
  sell = trades[timeBars,]$Sell %*% timeBarWeights
  
  imbalance <- abs(buy-sell)
  VolBars[i,] = c(i, mode, buy, sell,imbalance, imbalance/V)
}
#mode = expanded$State[which.max(expanded$State)]
print(VolBars)

data.df = data.frame(vpin = as.numeric(VolBars$VPIN), state=factor(VolBars$State), stringsAsFactors = FALSE)
boxplot(vpin ~ state, data = data.df, ylab = "VPIN Value")

model.lm = lm(vpin ~ state, data = data.df)
summary(model.lm)

model.aov = aov(vpin ~ state, data = data.df)
summary(model.aov)