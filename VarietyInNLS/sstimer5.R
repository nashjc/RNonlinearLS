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
   tsuml<-microbenchmark(sl<-suml(vv), unit="us")
   tsums<-microbenchmark(ss<-sums(vv), unit="us")
   tsumc<-microbenchmark(sc<-sumc(vv), unit="us")
   tcsuml<-microbenchmark(sl<-csuml(vv), unit="us")
   tcsums<-microbenchmark(ss<-csums(vv), unit="us")
   tcsumc<-microbenchmark(sc<-csumc(vv), unit="us")
   tsumld<-microbenchmark({sl<-0;for(ii in 1:nn){sl<-sl+vv[ii]^2}}, unit="us")
   tsumsd<-microbenchmark(ss<-sum(vv^2), unit="us")
   tsumcd<-microbenchmark(sc<-as.numeric(crossprod(vv)), unit="us")
   ml<-mean(tsuml$time)
   ms<-mean(tsums$time)
   mc<-mean(tsumc$time)
   mcl<-mean(tcsuml$time)
   mcs<-mean(tcsums$time)
   mcc<-mean(tcsumc$time)
   mld<-mean(tsumld$time)
   msd<-mean(tsumsd$time)
   mcd<-mean(tsumcd$time)
   cat(nn,"\t",ml," : ",ml/mc,"\t",ms," : ",ms/mc,"\t",mc,"\t",all.equal(sl, ss, sc),"\n")   
   cat("\t",mcl," : ",mcl/mcc,"\t",mcs," : ",mcs/mcc,"\t",mcc,"\n")   
   cat("\t",mld," : ",mld/mcd,"\t",mcs," : ",mcs/mcd,"\t",mcd,"\n")   
}
