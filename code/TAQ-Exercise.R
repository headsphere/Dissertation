library(highfrequency)

setwd("/Users/nick/Documents/RStudioProjects/Dissertation/code")

from = "2014-03-04";
to = "2014-03-04";
datasource = "/Users/nick/Documents/RStudioProjects/Dissertation/data/raw";
datadestination = "/Users/nick/Documents/RStudioProjects/Dissertation/data/xts";
tickers = c("ACW", "ALG","HZO","ITG")
tradeColNames <- c("DATE","TIME_M","EX","SYM_ROOT","SYM_SUFFIX","TR_SCOND","SIZE","PRICE","TR_STOPIND","TR_CORR","TR_SEQNUM","TR_SOURCE","TR_RF")
quoteColNames <- c("DATE","TIME_M","EX","SYM_ROOT","SYM_SUFFIX","BID","BIDSIZ","ASK","ASKSIZ","QU_COND","BIDEX","ASKEX","QU_SEQNUM","NATBBO_IND","NASDBBO_IND","QU_CANCEL","QU_SOURCE")

convert(from, to, datasource, datadestination, trades=TRUE,
        quotes=TRUE,ticker=tickers, dir=FALSE, extension="csv",
        header=TRUE,tradecolnames=tradeColNames,quotecolnames=quoteColNames,
        format="%Y%m%d %H:%M:%OS");
#             '%d/%m/%Y %H:%M:%OS

#TODO: check the source of the G127 NA errors

options("digits.secs"=3); #Show milliseconds
sample_tdata = TAQLoad(tickers="HZO", from=from, to=to, trades=T, quotes=F, 
                       datasource=datadestination)

sample_qdata = TAQLoad(tickers="HZO", from=from, to=to, trades=F, quotes=T, 
                       datasource=datadestination)

#head(sample_tdata)

sample_qdata = sample_qdata[!is.na(index(sample_qdata))]
sample_tdata = sample_tdata[!is.na(index(sample_tdata))]

#TODO: investigate R timezone handling (low priority)

tqdata = matchTradesQuotes(sample_tdata,sample_qdata);

plot.xts(x = tqdata[, "PRICE", "BID", "OFR"])