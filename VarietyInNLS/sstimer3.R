# crossprod timer
library(microbenchmark)
suml<-function(vv) {
    ss<-0.0
    for (i in 1:length(vv)) {ss<-ss+vv[i]^2}
    ss  
}
sums<-function(vv) {
  ss<-sum(vv^2)
  ss  
}
sumc<-function(vv) {
  ss<-as.numeric(crossprod(vv))
  ss  
}


## ll <- c(100, 1000, 10000, 100000, 1000000, 10000000)
ll <- c(100, 1000, 10000, 100000, 1000000)

cat(" n  \t  t(forloop) : ratio \t  t(sum) : ratio \t t(crossprod) \t all.equal \n")
for (nn in ll ){
   set.seed(1234)
   vv <- runif(nn)
   tsuml<-microbenchmark(sl<-suml(vv), unit="us", setup=suml(vv))
   tsums<-microbenchmark(ss<-sums(vv), unit="us", setup=sums(vv))
   tsumc<-microbenchmark(sc<-sumc(vv), unit="us", setup=sumc(vv))
   ml<-mean(tsuml$time)
   ms<-mean(tsums$time)
   mc<-mean(tsumc$time)
   cat(nn,"\t",ml," : ",ml/mc,"\t",ms," : ",ms/mc,"\t",mc,"\t",all.equal(sl, ss, sc),"\n")   
}
