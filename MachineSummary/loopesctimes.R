loopesc <- function(nn){
  ss <- 0
  for (i in 1:nn) {
    xx <- exp(sin(cos(1.0*i)))
    ss <- ss + xx
  }
  xx
}

require("microbenchmark")
# nvals<-c(1000, 10000, 100000, 1000000)
# tvals<-rep(NA, length(nvals))
# i<-0
# for (nn in nvals) {
#   i <- i + 1
#   cat(nn, "\n")
nn <- 10000
tt <- microbenchmark(loopesc(nn), unit = 'us')
print(tt)
ht<-hist(tt$time)
plot(ht, main="Histogram of loopesc times")
tmp<-readlines("next?")
graphics.off()
tt2 <- microbenchmark(loopesc(nn), unit = 'us', control=list(warmup=2))
print(tt2)
ht2<-hist(tt2$time)
plot(ht2, main="Histogram of loopesc times after 2 warmup cycles")
# }
# plot(log10(nvals), log(tvals))
