TAQLoadTrades <- function () {
  library(highfrequency)
  
  setwd("/Users/nick/Documents/RStudioProjects/Dissertation/code")
  
  from = "2014-03-04";
  to = "2014-03-04";
  datasource = "/Users/nick/Documents/RStudioProjects/Dissertation/data/raw";
  datadestination = "/Users/nick/Documents/RStudioProjects/Dissertation/data/xts";
  datadestination_cleaned = "/Users/nick/Documents/RStudioProjects/Dissertation/data/xts_cleaned";
  #tickers = c("ACW", "ALG","HZO","ITG")
  tickers = c("HZO")
  tradeColNames <- c("DATE","TIME_M","EX","SYMBOL","SYM_SUFFIX","COND","SIZE","PRICE","TR_STOPIND","CORR","TR_SEQNUM","TR_SOURCE","TR_RF")
  quoteColNames <- c("DATE","TIME_M","EX","SYMBOL","SYM_SUFFIX","BID","BIDSIZ","OFR","OFRSIZ","QU_COND","BIDEX","ASKEX","QU_SEQNUM","NATBBO_IND","NASDBBO_IND","QU_CANCEL","QU_SOURCE")
  
  convert(from, to, datasource, datadestination, trades=TRUE,
          quotes=TRUE,ticker=tickers, dir=FALSE, extension="csv",
          header=TRUE,tradecolnames=tradeColNames,quotecolnames=quoteColNames,
          format="%Y%m%d %H:%M:%OS");
  #             '%d/%m/%Y %H:%M:%OS
  
  #TODO: check the source of the G127 NA errors
  
  options("digits.secs"=3); #Show milliseconds
  tdata = TAQLoad(tickers="HZO", from=from, to=to, trades=T, quotes=F, 
                         datasource=datadestination)
  
  qdata = TAQLoad(tickers="HZO", from=from, to=to, trades=F, quotes=T, 
                         datasource=datadestination)
  
  #TODO: investigate R timezone handling (low priority)
  
  #tdata = autoSelectExchangeTrades(tdata = tdata)
  #qdata = autoSelectExchangeQuotes(qdata = qdata)
  
  tdata = tradesCleanup(tdataraw = tdata,exchanges="D", report = FALSE)
  qdata = quotesCleanup(qdataraw = qdata,exchanges="T", report = FALSE)
  
  # x[1:(length(x)-window)]
  
  tqdata = matchTradesQuotes(tdata,qdata);
  
  return(tdata)
}

#plot(x = tqdata[, "PRICE", "BID", "OFR"])

plot.xts(tqdata[, "PRICE"], type = "bars")
lines(tqdata[, "OFR"], col="red")
lines(tqdata[, "BID"], col="green")
