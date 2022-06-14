# loopbox.R
loopesc <- function(nn){
  ss <- 0
  for (i in 1:nn) {
    xx <- exp(sin(cos(1.0*i)))
    ss <- ss + xx
  }
  xx
}

require("microbenchmark")
nn <- 10000
library(compiler)
le <- cmpfun(loopesc)
tcompiled <- microbenchmark(le(nn), unit='us')$time
tcoldstart <- microbenchmark(loopesc(nn), unit = 'us')$time
twarmstart <- microbenchmark(loopesc(nn), unit = 'us', control=list(warmup=10))$time
tdf <- data.frame(tcoldstart, twarmstart, tcompiled)
boxplot(x=as.list(log(tdf)), main="log times for 10000 loops in microsecs")

