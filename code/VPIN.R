#### VPIN calculation #########################################################
#### install.packages('fasttime',repos='http://www.rforge.net/')
require(data.table)
require(fasttime)
require(plyr)
# Assuming TAQ data is arranged in 1 year stock csv files

stock = TAQLoadTrades()

# Now we have an xts data frame called 'stock' with a DateTime index and...  two columns: Price and Volume
# Vbucket=Number of volume buckets in an average volume day (Vbucket=50)
VPIN = function(stock, Vbucket) {
  stock$dP1 = diff(stock[, "Price"], lag = 1, diff = 1, na.pad = TRUE)
  ends = endpoints(stock, "minutes")
  timeDF = period.apply(stock[, "dP1"], INDEX = ends, FUN = sum)
  timeDF$Volume = period.apply(stock[, "Size"], INDEX = ends, FUN = sum)
  Vbar = mean(period.apply(timeDF[, "Size"], INDEX = endpoints(timeDF, "days"), FUN = sum))/Vbucket
  timeDF$Vfrac = timeDF[, "Size"]/Vbar
  timeDF$CumVfrac = cumsum(timeDF[, "Vfrac"])
  timeDF$Next = (timeDF[, "CumVfrac"] - floor(timeDF[, "CumVfrac"]))/timeDF[, "Vfrac"]
  timeDF[timeDF[, "Next"] < 1, "Next"] = 0
  timeDF$Previous = lag(timeDF[, "dP1"]) * lag(timeDF[, "Next"])
  timeDF$dP2 = (1 - timeDF[, "Next"]) * timeDF[, "dP1"] + timeDF[, "Previous"]
  timeDF$Vtick = floor(timeDF[, "CumVfrac"])
  timeDF[, "Vtick"] = timeDF[, "Vtick"] - diff(timeDF[, "Vtick"])
  timeDF[1, "Vtick"] = 0
  timeDF = as.data.frame(timeDF)
  timeDF[, "DateTime"] = row.names(timeDF)
  timeDF = ddply(as.data.frame(timeDF), .(Vtick), last)
  timeDF = as.xts(timeDF[, c("Size", "dP2", "Vtick")], order.by = fastPOSIXct(timeDF$DateTime, tz = "GMT"))
  timeDF[1, "dP2"] = 0
  timeDF$sigma = rollapply(timeDF[, "dP2"], Vbucket, sd, fill = NA)
  timeDF$sigma = na.fill(timeDF$sigma, "extend")
  timeDF$Vbuy = Vbar * pnorm(timeDF[, "dP2"]/timeDF[, "sigma"])
  timeDF$Vsell = Vbar - timeDF[, "Vbuy"]
  timeDF$OI = abs(timeDF[, "Vsell"] - timeDF[, "Vbuy"])
  timeDF$VPIN = rollapply(timeDF[, "OI"], Vbucket, sum)/(Vbar * Vbucket)
  timeDF = timeDF[, c("VPIN")]
  return(timeDF)
}
out = VPIN(stock, 50) 