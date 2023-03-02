?microbenchmark
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopesctimes.R")
library(benchmarkme)
get_sys_details
sessionInfo
sessionInfo()
?COMPILE
library(compiler)
?compile
cmpfun(loopesc)
le<-cmpfun(loopesc)
microbenchmark(le)
microbenchmark(le, unit=us)
microbenchmark(le, unit="us")
microbenchmark(loopesc, unit="us")
nn
microbenchmark(le(nn), unit="us")
microbenchmark(loopesc(nn), unit="us")
microbenchmark()
?microbenchmark
boxplot.microbenchmark(t2)
boxplot.microbenchmark(tt2)
library(microbenchmark)
?boxplot.microbenchmark
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
#  print(te)
#  boxplot(te)
tcoldstart <- microbenchmark(loopesc(nn), unit = 'us')$time
#  print(tt)
#  boxplot(tt)
#  ht<-hist(tt$time)
#  plot(ht, main="Histogram of loopesc times")
twarmstart <- microbenchmark(loopesc(nn), unit = 'us', control=list(warmup=2))$time
tdf <- data.frame(tcompiled, tcoldstart, twarmstart)
str(tdf)
names(tdf)
?data.frame
boxplot(x=as.list(tdf))
boxplot(x=as.list(log(tdf))
)
?boxplot
tdf <- data.frame(tcompiled, tcoldstart, twarmstart)
boxplot(x=as.list(log(tdf)), main="log times for 10000 loops in microsecs")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
source("~/current/GSoC2021/RNonlinearLS/MachineSummary/loopbox.R")
?nlminb
?qt
curve(qt(p), from=-4, to=4)
curve(qt(x), from=-4, to=4)
curve(qt(x, df=5), from=-4, to=4)
curve(qt(x, df=5), from=0, to=1)
install.packages("pacman")
source("~/current/GSoC2021/RNonlinearLS/DNase1Models.R")
DNase1 <- subset(DNase, Run == 1) # select the data
## using a selfStart model - do not specify the starting parameters
fm1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
summary(fm1)
## using conditional linearity - leave out the Asym parameter
fm2 <- nls(density ~ 1/(1 + exp((xmid - log(conc))/scal)),
data = DNase1, start = list(xmid = 0, scal = 1),
algorithm = "plinear")
summary(fm2)
## using conditional linearity AND Asym does NOT work
fm2a <- try(nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
data = DNase1, start = list(Asym=3, xmid = 0, scal = 1),
algorithm = "plinear", trace = TRUE))
source("~/current/GSoC2021/RNonlinearLS/DNase1Models.R")
DNase1 <- subset(DNase, Run == 1) # select the data
## using a selfStart model - do not specify the starting parameters
fm1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
summary(fm1)
## using conditional linearity - leave out the Asym parameter
fm2 <- nls(density ~ 1/(1 + exp((xmid - log(conc))/scal)),
data = DNase1, start = list(xmid = 0, scal = 1),
algorithm = "plinear")
summary(fm2)
## without conditional linearity
fm3 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
data = DNase1,
start = list(Asym = 3, xmid = 0, scal = 1))
summary(fm3)
## using Port's nl2sol algorithm
fm4 <- try(nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1),
algorithm = "port"))
summary(fm4)
## using conditional linearity AND Asym does NOT work
fm2a <- try(nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
data = DNase1, start = list(Asym=3, xmid = 0, scal = 1),
algorithm = "plinear", trace = TRUE))
summary(fm2a)
# bhobbsX.R
#### bounds with formula specification of problem using Hobbs Weed problem
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weedframe <- data.frame(y=weed, tt=tt)
st <- c(b1=1, b2=1, b3=1) # a default starting vector (named!)
wmods <- y ~ 100*b1/(1+10*b2*exp(-0.1*b3*tt)) ## Scaled model
require(minpack.lm)
require(nlsr)
traceval<-TRUE
# Hobbs scaled problem with bounds, formula specification
anlxb1 <- nlxb(wmods, start=start1, data=weedframe, lower=c(0,0,0),
upper=c(2,6,3), trace=traceval)
print(anlxb1)
# Hobbs scaled problem with bounds, formula specification
anlxb1 <- nlxb(wmods, start=st, data=weedframe, lower=c(0,0,0),
upper=c(2,6,3), trace=traceval)
print(anlxb1)
## nlsLM seems NOT to work properly with bounds
anlsLM1b <- nlsLM(wmods, start=st, data=weedframe, lower=c(0,0,0),
upper=c(2,6,3), trace=traceval)
print(anlsLM1b)
source("~/current/GSoC2021/RNonlinearLS/RefactoringNLS/bhobbsXsmall.R")
source("~/current/GSoC2021/RNonlinearLS/RefactoringNLS/bhobbsXsmall.R")
summary(anlxb1)
print(anlxb1)
str(anlxb1)
source("~/current/GSoC2021/RNonlinearLS/RefactoringNLS/bhobbsXsmall.R")
source("~/current/GSoC2021/RNonlinearLS/RefactoringNLS/bhobbsXsmall.R")
source("~/current/GSoC2021/RNonlinearLS/RefactoringNLS/bhobbsX.R")
source("~/current/GSoC2021/RNonlinearLS/RefactoringNLS/bhobbsXsmall.R")
source("~/current/GSoC2021/RNonlinearLS/RefactoringNLS/bhobbsXsmall.R")
# Meyer_10 function as a least squares problem
y <- c(34780, 28610, 23650, 19630, 16370, 13720, 11540, 9744, 8261, 7030,
6005, 5147, 4427, 3820, 3307, 2872)
m <- 16
t <- 45 + 5 * (1:m)
df <- data.frame(t, y)
modl <- y ~ x1 * exp(x2/(t + x3))
library(minpack.lm)
library(nlsr)
cat("nls:\n")
anls<-nls(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anls)
anlxb<-nlxb(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anlxb)
anlxb
anlsLM<-nlsLM(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anlsLM)
library(nlsj)
anlsj<-nlsj(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anlsj)
library(nlsralt)
anlxbx<-nlxbx(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anlxbx) # gets to min but slowly
anlxbx
anlxb
anlsp<-nls(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, algorithm="plinear", trace=TRUE)
summary(anlsp) # fails singular gradient
cat("NOTE: modlp NOT modl\n")
anlsp<-nls(formula=modlp, start=c(x1=1, x2=1, x3=1), data=df, algorithm="plinear", trace=TRUE)
modlp <- y ~ exp(x2/(t + x3))
anls<-nls(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anls) # fails singular gradient
cat("NOTE: modlp NOT modl\n")
anlsp<-nls(formula=modlp, start=c(x1=1, x2=1, x3=1), data=df, algorithm="plinear", trace=TRUE)
summary(anlsp) # fails singular gradient
anlsp<-nls(formula=modlp, start=c(x2=1, x3=1), data=df, algorithm="plinear", trace=TRUE)
summary(anlsp) # fails singular gradient
anlsp
library(nlsj)
anlsj<-nlsj(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anlsj) # NOT good
anlsj<-nlsjm(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE, algorithm="marquardt")
anlsjm<-nlsj(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE, algorithm="marquardt")
summary(anlsjm) # NOT good
anlsjm
anlsLM
anlxb
anlsjm2<-nlsj(formula=modl, start=c(x1=1, x2=100, x3=1), data=df, trace=TRUE, algorithm="marquardt")
summary(anlsjm2) # NOT good
anlxb
anlsjm3<-nlsj(formula=modl, start=c(x1=.01, x2=100, x3=1), data=df, trace=TRUE, algorithm="marquardt")
summary(anlsjm3) # NOT good
anlsjm3
anlxb
anlsjm3<-nlsj(formula=modl, start=c(x1=.01, x2=1000, x3=100), data=df, trace=TRUE, algorithm="marquardt")
summary(anlsjm3) # NOT good
anlsjm3
?sys.calls
library(minpack.lm)
?wfct
?on.exit
?missing
?match.call
?n
?sapply
?nls.lm
?cmpfun
source("~/current/GSoC2021/CmpDerivs.R")
knitr::opts_chunk$set(echo = TRUE)
## require(bookdown) # language engine to display text - does not seem necessary
reala <- -3; realb <- 5; realc <- 0.5; realr <- 0.7; realm <- 1 # underlying parameters
x <- 1:11 # x values; 11 data points
y <- reala + realb * realr^(x - realm) + realc * x # linear + exponential function
testdat <- data.frame(x, y) # save in a data frame
strt <- list(a = -3, b = 5, c = 0.5, r = 0.7, m = 1) # give programs a good start
jform <- y ~ a + b * r^(x - m) + c * x # Formula
library(nlsr)
linexp2 <- try(nlxb(jform, data = testdat, start = strt, trace = F))
linexp2 # Note singular values of Jacobian in rightmost column
# We get an error with nls()
linexp <- try(nls(jform, data = testdat, start = strt, trace = F))
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
wts <- 0.5^tt # simple weights
nowt<-nls(weed ~ SSlogis, data=weeddf)
nowt<-nls(weed ~ SSlogis(Asym, xmid, scal), data=weeddf)
nowt<-nls(weed ~ SSlogis(), data=weeddf)
?SSlogis
nowt<-nls(weed ~ SSlogis(tt))
Asym<-1; xmid<-1; scal<-1
nowt<-nls(weed ~ SSlogis(tt, Asym, xmid, scal))
nowt
usewt <- nls(weed ~ SSlogis(tt, Asym, xmid, scal), weights=wts)
usewt
nowt$resid
nowt$resid()
str(nowt)
nowt$m$resid()
usewt <- nls(weed ~ SSlogis(tt, Asym, xmid, scal), weights=wts)
usewt
usewt$m$resid()
library(nlsr2)
nowtx<-nlxb(weed ~ SSlogis(tt, Asym, xmid, scal))
nowtx<-nlxb(weed ~ SSlogis(tt, Asym, xmid, scal), start=c(Asym=1, xmid=1, scal=1))
nowtx<-nlxb(weed ~ SSlogis(tt, Asym, xmid, scal), start=c(Asym=1, xmid=1, scal=1), control=list(japprox="jafwd"))
nowtx
nowtx<-nlxb(weed ~ SSlogis(tt, Asym, xmid, scal), start=c(Asym=1, xmid=1, scal=1), control=list(japprox="jafwd"), trace=TRUE)
inpar<-getInitial(weed ~ SSlogis(tt, Asym, xmid, scal))
inpar<-getInitial(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf)
print(inpar)
str(nowtx)
usewtx<-nlxb(weed ~ SSlogis(tt, Asym, xmid, scal), start=c(Asym=1, xmid=1, scal=1), weights=wts,
control=list(japprox="jafwd"), trace=TRUE)
usewtx
usewt
usewtx
wts
usewtx<-nlxb(weed ~ SSlogis(tt, Asym, xmid, scal), start=inpar, weights=wts,
control=list(japprox="jafwd"), trace=TRUE)
usewtx
usewt$m$resid()
usewtx<-nlxb(weed ~ SSlogis(tt, Asym, xmid, scal), start=inpar, weights=wts,
control=list(japprox="jafwd"), trace=TRUE)
ls()
usewtx
usewtx$resid()
usewtx$resid
nowtx$resid
nmod <- nlsModel(weed ~ SSlogis(tt, Asym, xmid, scal), weights=wts)
nmod <- nlsModel(weed ~ SSlogis(tt, Asym, xmid, scal))
nmod <- nlsModel(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf, weights=wts)
nmod <- nlsModel(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf, start=c(Asym=1, xmid=1, scal=1), weights=wts)
?nlsModel
nmod <- nlsModel(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf, wts=wts)
nmod <- nlsModel(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf, start=c(Asym=1, xmid=1, scal=1), wts=wts)
nmod <- nlsModel(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf, start=c(Asym=1, xmid=1, scal=1), wts=wts)
nmod
nmod$resid()
usewt$m$resid()
usewtx<-nlxb(weed ~ SSlogis(tt, Asym, xmid, scal), start=inpar, weights=wts,
control=list(japprox="jafwd"), trace=TRUE)
?SSlogis
SSlogis
frmlogis <- y ~ Asym/(1 + exp((xmid - x)/scal))
knitr::opts_chunk$set(echo = TRUE)
## require(bookdown) # language engine to display text - does not seem necessary
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
wts <- 0.5^tt # simple weights
frmlogis <- y ~ Asym/(1 + exp((xmid - x)/scal))
Asym<-1; xmid<-1; scal<-1
nowt<-nls(weed ~ SSlogis(tt, Asym, xmid, scal))
nowt
nowt$m$resid()
inpar<-getInitial(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf)
print(inpar)
library(nlsr2)
nowtx<-nlxb(frmlogis, start=c(Asym=1, xmid=1, scal=1), trace=TRUE)
nowtx
nowtx<-nlxb(frmlogis, start=c(Asym=1, xmid=1, scal=1), trace=TRUE)
frmlogis <- weed ~ Asym/(1 + exp((xmid - tt)/scal))
nowtx<-nlxb(frmlogis, start=c(Asym=1, xmid=1, scal=1), trace=TRUE)
nowtx
usewt <- nls(weed ~ SSlogis(tt, Asym, xmid, scal), weights=wts)
usewt
usewt$m$resid()
nowtxi<-nlxb(frmlogis, start=inpar, trace=TRUE)
nowtxi
?getInitial
usewtx<-nlxb(weed ~ SSlogis(tt, Asym, xmid, scal), start=inpar, weights=wts,
control=list(japprox="jafwd"), trace=TRUE)
usewtx<-nlxb(frmlogis,  start=c(Asym=1, xmid=1, scal=1), weights=wts,
control=list(japprox="jafwd"), trace=TRUE)
usewtx
nowtxi<-nlxb(frmlogis, start=inpar, trace=TRUE)
nowtx<-nlxb(frmlogis, start=c(Asym=1, xmid=1, scal=1), trace=TRUE)
usewtx<-nlxb(frmlogis,  start=c(Asym=1, xmid=1, scal=1), weights=wts,
trace=TRUE)
usewtx
usewtxi<-nlxb(frmlogis, start=inpar, weights=wts,
weights=wts, trace=TRUE)
usewtxi<-nlxb(frmlogis, start=inpar, weights=wts, trace=TRUE)
usewtxi
usewtxn<-nlxb(frmlogis,  start=c(Asym=1, xmid=1, scal=1), weights=wts, trace=TRUE, control=list(japprox="jafwd"))
usewtxn<-nlxb(frmlogis,  start=c(Asym=1, xmid=1, scal=1), weights=wts, trace=TRUE, control=list(japprox="jand"))
usewtxn
usewtxn<-nlxb(frmlogis,  start=c(Asym=1, xmid=1, scal=1), weights=wts, trace=TRUE, control=list(japprox="jacentral"))
usewtxn
usewtxn<-nlxb(frmlogis,  start=c(Asym=1, xmid=1, scal=1), weights=wts, trace=TRUE, control=list(japprox="jaback"))
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
wts <- 0.5^tt # simple weights
frmlogis <- weed ~ Asym/(1 + exp((xmid - tt)/scal))
Asym<-1; xmid<-1; scal<-1
nowt<-nls(weed ~ SSlogis(tt, Asym, xmid, scal))
nowt
nowt$m$resid()
inpar<-getInitial(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf)
print(inpar)
library(nlsr2)
nowtx<-nlxb(frmlogis, start=c(Asym=1, xmid=1, scal=1), trace=TRUE)
nowtx
nowtxi<-nlxb(frmlogis, start=inpar, trace=TRUE)
nowtxi
usewt <- nls(weed ~ SSlogis(tt, Asym, xmid, scal), weights=wts)
usewt
usewt$m$resid()
nmod <- nlsModel(frmlogis, data=weeddf, start=c(Asym=1, xmid=1, scal=1), wts=wts)
usewtx<-nlxb(frmlogis,  start=c(Asym=1, xmid=1, scal=1), weights=wts, trace=TRUE)
usewtx
# NOTE jafwd and jaback fail!
usewtxn<-nlxb(frmlogis,  start=c(Asym=1, xmid=1, scal=1), weights=wts, trace=TRUE, control=list(japprox="jacentral"))
usewtxn
usewtxi<-nlxb(frmlogis, start=inpar, weights=wts, trace=TRUE)
usewtxi
# usewtx$resid
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
wts <- 0.5^tt # simple weights
frmlogis <- weed ~ Asym/(1 + exp((xmid - tt)/scal))
Asym<-1; xmid<-1; scal<-1
nowt<-nls(weed ~ SSlogis(tt, Asym, xmid, scal))
nowt
nowt$m$resid()
inpar<-getInitial(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf)
print(inpar)
library(nlsr2)
nowtx<-nlxb(frmlogis, start=c(Asym=1, xmid=1, scal=1), trace=TRUE)
nowtx
nowtSSx<-nlxb(weed ~ SSlogis(tt, Asym=1, xmid=1, scal=1), start=c(Asym=1, xmid=1, scal=1), trace=TRUE)
nowtSSx<-nlxb(weed ~ SSlogis(tt, Asym=1, xmid=1, scal=1), start=c(Asym=1, xmid=1, scal=1), trace=TRUE,
control=list(japprox="jacentral"))
nowtSsx
nowtSSx
nowtSSx<-nlxb(weed ~ SSlogis(tt, Asym=1, xmid=1, scal=1), start=c(Asym=1, xmid=1, scal=1), trace=TRUE,
control=list(japprox="jacentral"))
nowtSSx
nowtSSx<-nlxb(weed ~ SSlogis(tt, Asym=1, xmid=1, scal=1), start=c(Asym=1, xmid=1, scal=1), trace=TRUE,
control=list(japprox="jand"))
nowtSSx<-nlxb(weed ~ SSlogis(tt, Asym=1, xmid=1, scal=1), start=c(Asym=100, xmid=2, scal=1), trace=TRUE,
control=list(japprox="jand"))
nowtSSx
?droplevels
library(nlsr2)
knitr::opts_chunk$set(echo = TRUE)
## require(bookdown) # language engine to display text - does not seem necessary
Asym<-1; xmid<-1; scal<-1
strt <- c(Asym=Asym, xmid=xmid, scal=scal)
modlogis<-nlsr2::model2rjfun(frmlogis, strt)
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
wts <- 0.5^tt # simple weights
frmlogis <- weed ~ Asym/(1 + exp((xmid - tt)/scal))
Asym<-1; xmid<-1; scal<-1
strt <- c(Asym=Asym, xmid=xmid, scal=scal)
modlogis<-model2rjfun(frmlogis, strt)
Asym<-1; xmid<-1; scal<-1
nowt<-nls(weed ~ SSlogis(tt, Asym, xmid, scal))
nowt
nowt$m$resid() # This has WEIGHTED residual and Jacobian
modlogis(coef(nowt))
inpar<-getInitial(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf)
print(inpar)
modlogis(coef(nowt))
nlsrmod<-modlogis(coef(nowt))
nlsrmod
library(nlsr2)
knitr::opts_chunk$set(echo = TRUE)
## require(bookdown) # language engine to display text - does not seem necessary
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
wts <- 0.5^tt # simple weights
frmlogis <- weed ~ Asym/(1 + exp((xmid - tt)/scal))
Asym<-1; xmid<-1; scal<-1
strt <- c(Asym=Asym, xmid=xmid, scal=scal)
modlogis<-model2rjfun(frmlogis, strt)
Asym<-1; xmid<-1; scal<-1
nowt<-nls(weed ~ SSlogis(tt, Asym, xmid, scal))
nowt
nlsmod0<-nowt$m$resid() # This has WEIGHTED residual and Jacobian
nlsmod0
nlsmod<-nowt$m$resid(coef(nowt)) # This has WEIGHTED residual and Jacobian
nlsrmod<-modlogis(coef(nowt))
nlsrmod
nowt
nowtx<-nlxb(frmlogis, start=c(Asym=1, xmid=1, scal=1), trace=FALSE)
nowtx
# nowtSSx0<-nlxb(weed ~ SSlogis(tt, Asym=1, xmid=1, scal=1), start=c(Asym=1, xmid=1, scal=1), trace=TRUE)
# nowtSsx0
# nowtSSx<-nlxb(weed ~ SSlogis(tt, Asym=1, xmid=1, scal=1), start=c(Asym=100, xmid=2, scal=1), trace=TRUE,
#               control=list(japprox="jand"))
# nowtSSx
nowtxi<-nlxb(frmlogis, start=inpar, trace=TRUE)
nowtxi
inpar<-getInitial(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf)
print(inpar)
nowtx<-nlxb(frmlogis, start=c(Asym=1, xmid=1, scal=1), trace=FALSE)
nowtx
# nowtSSx0<-nlxb(weed ~ SSlogis(tt, Asym=1, xmid=1, scal=1), start=c(Asym=1, xmid=1, scal=1), trace=TRUE)
# nowtSsx0
# nowtSSx<-nlxb(weed ~ SSlogis(tt, Asym=1, xmid=1, scal=1), start=c(Asym=100, xmid=2, scal=1), trace=TRUE,
#               control=list(japprox="jand"))
# nowtSSx
nowtxi<-nlxb(frmlogis, start=inpar, trace=TRUE)
nowtxi
usewt <- nls(weed ~ SSlogis(tt, Asym, xmid, scal), weights=wts)
usewt
usewt$m$resid()
nmod <- nlsModel(frmlogis, data=weeddf, start=c(Asym=1, xmid=1, scal=1), wts=wts)
nmod
nmod$m$resid()
nmod$m$resid(coef(usewt))
nmod$m$resid
str(nmod)
nmod$resid()
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
wts <- 0.5^tt # simple weights
Asym<-1; xmid<-1; scal<-1
nowt<-nls(weed ~ SSlogis(tt, Asym, xmid, scal)) # UNWEIGHTED
nowt
nowt$m$resid() # This has UNWEIGHTED residual and Jacobian. Does NOT take coefficients.
usewt <- nls(weed ~ SSlogis(tt, Asym, xmid, scal), weights=wts)
usewt
usewt$m$resid() # WEIGHTED. Does NOT take coefficients.
nmod <- nlsModel(frmlogis, data=weeddf, start=c(Asym=1, xmid=1, scal=1), wts=wts)
nmod$resid()
nmod$resid()*wts
nmod$resid(coef(usewt))
nmod <- nlsModel(frmlogis, data=weeddf, start=coef(usewt), wts=wts)
nmod$resid()
nmod$resid()*wts
?getInitial
source("~/current/GSoC2021/improvenls/ExamplesVignette/hobbsdata.R", echo=TRUE)
ihobb<-getInitial(weed~100*b1/(1+10*b2*exp(-0.1*b3*tt)), weeddf)
ihobb<-getInitial(weed~SSlogis(TT, Asym, xmid, scal), weeddf)
ihobb<-getInitial(weed~SSlogis(tt, Asym, xmid, scal), weeddf)
ibobb
ihobb
fhobb<-nls(weed~SSlogis(tt, Asym, xmid, scal), weeddf)
fhobb
savehistory("hobbSSlogistry.R")
