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

csuml<-compiler::cmpfun(suml)
csumc<-compiler::cmpfun(sumc)
csums<-compiler::cmpfun(sums)

# ll <- c(100, 1000, 10000, 100000, 1000000, 10000000)
ll <- c(100, 1000, 10000, 100000, 1000000)

cat(" n  \t  t(forloop) : ratio \t  t(sum) : ratio \t t(crossprod) \t all.equal \n")
for (nn in ll ){
   set.seed(1234)
   vv <- runif(nn)
   tsuml<-microbenchmark(sl<-csuml(vv), unit="us")
   tsums<-microbenchmark(ss<-csums(vv), unit="us")
   tsumc<-microbenchmark(sc<-csumc(vv), unit="us")
   ml<-mean(tsuml$time)
   ms<-mean(tsums$time)
   mc<-mean(tsumc$time)
   cat(nn,"\t",ml," : ",ml/mc,"\t",ms," : ",ms/mc,"\t",mc,"\t",all.equal(sl, ss, sc),"\n")   
}
