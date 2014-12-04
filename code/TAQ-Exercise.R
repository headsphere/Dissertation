library(highfrequency)

setwd("/Users/nick/Documents/RStudioProjects/Dissertation/code")

from = "2014-03-04";
to = "2014-03-04";
datasource = "/Users/nick/Documents/RStudioProjects/Dissertation/data/raw";
datadestination = "/Users/nick/Documents/RStudioProjects/Dissertation/data/xts";
tickers = c("ACW", "ALG","HZO","ITG")
tradeColNames <- c("SYMBOL","DATE","TIME_M","EX","TR_SCOND","SIZE","PRICE")
quoteColNames <- c("DATE","TIME_M","EX","SYM_ROOT","SYM_SUFFIX","BID","BIDSIZ","ASK","ASKSIZ","QU_COND","BIDEX","ASKEX","QU_SEQNUM","NATBBO_IND","NASDBBO_IND","QU_CANCEL","QU_SOURCE")

convert(from, to, datasource, datadestination, trades=FALSE,
        quotes=TRUE,ticker=tickers, dir=FALSE, extension="csv",
        header=TRUE,tradecolnames=tradeColNames,quotecolnames=quoteColNames,
        format="%Y%m%d %H:%M:%S.%OS");