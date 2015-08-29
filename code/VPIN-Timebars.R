### Simulated input
source('VPIN-Sim.R')
V = 55 #size of each volume bucket
n = 100 #number of time bars
alpha = 0.28 #probability of information event
delta = 0.33 #probability of bad news event
mu = 80 #arrival rate of informed trader
epsilon = (V - alpha * mu)/2 #arrival rate of uninformed trader
trades = generate.trades.sim(n, alpha, delta, mu, epsilon) 

### Empirical input
states = readRDS("../data/decoded-states.RDS")
deser = readRDS("../data/SPY_BuysSells.RDS")
V = 55
n = length(states)
trades = data.frame(Buckets=seq(1:n), 
                    Buy=coredata(deser)[,2], 
                    Sell=coredata(deser)[,3],
                    States=states,
                    TimeBar=index(deser))

### Test Input from spreadsheet example
# V=10
# n=10
# trades = data.frame(Buckets=seq(1:10), 
#                     Buy=c(1,2,2,3,11,3,5,2,2,1), 
#                     Sell=c(2,4,5,7,3,6,2,4,4,1), 
#                     States=c(1,2,3,2,3,3,2,1,2,1))

library(data.table)
trades=data.frame(trades, Total=trades$Buy+trades$Sell)
totalVol = sum(trades$Total)
expanded = data.table(Timebar= integer(totalVol),
                      TimeState=integer(totalVol),
                      Volume=integer(totalVol),
                      TimebarStart= integer(totalVol))
startBlock = 1
for(i in 1:n)
{
  print(paste(i,"of",n,"trades processing"))
  repeats = trades[i,]$Total
  state = trades[i,]$States
  timeBucket = trades[i,]$TimeBar
  endBlock = startBlock + repeats -1
  expanded[startBlock:endBlock,1 := i]
  expanded[startBlock:endBlock,2 := state] 
  expanded[startBlock:endBlock,3 := repeats] 
  expanded[startBlock:endBlock,4 := as.ITime(timeBucket)] #just stores the integer repesentation of the seconds component of the day
  startBlock = endBlock + 1
}

#as.POSIXct("2010-05-06",tz = "GMT") + expanded$Timebar[1]

noOfVolBars <- floor(totalVol / V)
VolBars = data.frame(VolBar = rep(NA,noOfVolBars), 
                     State = rep(NA,noOfVolBars), 
                     Buy = rep(NA,noOfVolBars), 
                     Sell = rep(NA,noOfVolBars), 
                     Imbalance = rep(NA,noOfVolBars), 
                     VPIN = rep(NA,noOfVolBars),
                     TimeStart=rep(as.POSIXct(NA,"",tz = "GMT"),noOfVolBars),
                     stringsAsFactors = FALSE)
for(i in 1:noOfVolBars)
{
  print(paste(i,"of",noOfVolBars,"volbars processing"))
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

  VolBars[i,] = c(i, mode, buy, sell,imbalance, imbalance/V,as.POSIXct(NA,"",tz = "GMT"))
  VolBars$TimeStart[i] = as.POSIXct(volData$TimebarStart[1], origin = "2010-05-06 00:00:00", tz = "GMT")
}
#mode = expanded$State[which.max(expanded$State)]
print(VolBars)

plot.xts(zoo(VolBars$VPIN, VolBars$TimeStart))

data.df = data.frame(vpin = as.numeric(VolBars$VPIN), 
                     state=factor(VolBars$State), 
                     stringsAsFactors = FALSE)

boxplot(vpin ~ state, data = data.df, ylab = "VPIN Value")

model.lm = lm(vpin ~ state, data = data.df)
summary(model.lm)

model.aov = aov(vpin ~ state, data = data.df)
summary(model.aov)
