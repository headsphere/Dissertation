# ---- taq.load ----
library(zoo)
library(xts)
library(highfrequency)
library(Defaults)
library(TTR)
library(quantmod)
library(timeDate)

conversionRequired = FALSE;
setwd("/Users/nick/Documents/RStudioProjects/Dissertation/code")
datasource = "/Users/nick/Documents/RStudioProjects/Dissertation/data/raw";
datadestination = "/Users/nick/Documents/RStudioProjects/Dissertation/data/xts";
datadestination_cleaned = "/Users/nick/Documents/RStudioProjects/Dissertation/data/xts_cleaned";

from = "2010-05-06";
to = "2010-05-06";

tickers = c("SPY")

tradeColNames <- c("SYMBOL","DATE","TIME","PRICE","SIZE","G127","CORR","COND","EX")

quoteColNames <- c("SYMBOL","DATE","TIME","BID","OFR","BIDSIZ","OFRSIZ","MODE","EX","MMID")

if(conversionRequired){
  #TODO: check the source of the G127 NA errors
  print("Beginning conversion")
  suppressWarnings(
    convert(from, to, datasource, datadestination, trades=TRUE,
            quotes=TRUE,ticker=tickers, dir=FALSE, extension="csv",
            header=TRUE,tradecolnames=tradeColNames,quotecolnames=quoteColNames,
            format="%Y%m%d %H:%M:%S")
  )
  print("Completed conversion")
}

# options("digits.secs"=3); #Shows milliseconds
print("Loading Trades")
tdata = TAQLoad(tickers=tickers, from=from, to=to, trades=T, quotes=F, 
                datasource=datadestination)

print("Loading Quotes")
qdata = TAQLoad(tickers=tickers, from=from, to=to, trades=F, quotes=T, 
                datasource=datadestination)

print("Cleaning up quotes")
qdata = quotesCleanup(qdataraw = qdata,exchanges="T", report = FALSE, maxi=25) 

print("Cleaning up trades")
tdataAfterFinalCleanup = tradesCleanupFinal(qdata = qdata, tdata = tdata)

print("Matching trades and quotes")
tqdata = matchTradesQuotes(tdataAfterFinalCleanup,qdata)

print("Getting trade directions using Lee-Ready")
tradeDirection = getTradeDirection(tqdata)

core = coredata(tradeDirection)
n = length(tradeDirection)
buy = rep(0, n)
sell = rep(0, n)
buy[which(core[1:n] == 1)] = 1
sell[which(core[1:n] == -1)] = 1

print("Calculating time bars")
buys = zoo(buy, index(tqdata))
sells = zoo(sell, index(tqdata))
buys.ts = period.sum(buys, endpoints(buys, "seconds", 1))
sells.ts = period.sum(sells, endpoints(sells, "seconds", 1))
buysSells.ts = merge.xts(buys.ts,sells.ts) 
buysSells.ts = merge.xts(tqdata$PRICE, buysSells.ts, join = 'inner') 

plot.zoo(buysSells.ts)
saveRDS(buysSells.ts, "../data/SPY_BuysSells.RDS")
