# Generated by `rjournal_pdf_article()` using `knitr::purl()`: do not edit by hand
# Please edit ComparingNLStools.Rmd to modify this file

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)


## ----ex01, echo=TRUE----------------------------------------------------------
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
plot(weeddf, main="Hobbs weed infestation data")


## ----ex02set, echo=TRUE-------------------------------------------------------
# model formulas
frmu <- weed ~ b1/(1+b2*exp(-b3*tt))
frms <- weed ~ 100*c1/(1+10*c2*exp(-0.1*c3*tt))
frmt <- weed ~ Asym /(1 + exp((xmid-tt)/scal))
#
# Starting parameter sets
stu1<-c(b1=1, b2=1, b3=1)
sts1<-c(c1=1, c2=1, c3=1)
stt1<-c(Asym=1, xmid=1, scal=1)


## ----ex02fn, echo=TRUE--------------------------------------------------------
# Logistic3U
hobbs.res  <-  function(x){ # scaled Hobbs weeds problem -- residual
  # This variant uses looping
  if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
  y  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
           38.558, 50.156, 62.948, 75.995, 91.972)
  tt  <-  1:12
  res  <-  x[1]/(1+x[2]*exp(-x[3]*tt)) - y
}

hobbs.jac  <-  function(x) { # scaled Hobbs weeds problem -- Jacobian
  jj  <-  matrix(0.0, 12, 3)
  tt  <-  1:12
  yy  <-  exp(-x[3]*tt)
  zz  <-  1.0/(1+x[2]*yy)
  jj[tt,1]   <-   zz
  jj[tt,2]   <-   -x[1]*zz*zz*yy
  jj[tt,3]   <-   x[1]*zz*zz*yy*x[2]*tt
  attr(jj, "gradient") <- jj
  jj
}


## ----shortforms, echo=FALSE---------------------------------------------------
library(nlsr)
library(minpack.lm)


## ----ex02nls------------------------------------------------------------------
unls1<-try(nls(formula=frmu, start=stu1, data=weeddf))
snls1<-try(nls(formula=frms, start=sts1, data=weeddf))
tnls1<-try(nls(formula=frmt, start=stt1, data=weeddf))


## ----ex02nlsr-----------------------------------------------------------------
unlx1<-try(nlxb(formula=frmu, start=stu1, data=weeddf))
print(unlx1) 
snlx1<-try(nlxb(formula=frms, start=sts1, data=weeddf))
pshort(snlx1) # a short-form output
tnlx1<-try(nlxb(formula=frmt, start=stt1, data=weeddf))
pshort(tnlx1) # alternatively print(tnlx1)


## ----ex02minpack--------------------------------------------------------------
unlm1<-try(nlsLM(formula=frmu, start=stu1, data=weeddf))
pnls(unlm1)  # Short form of output
snlm1<-try(nlsLM(formula=frms, start=sts1, data=weeddf))
pnls(snlm1)
tnlm1<-try(nlsLM(formula=frmt, start=stt1, data=weeddf))
pnls(tnlm1) # short form to give sum of squares, else use summary(tnlm1)


## ----ex02gslnls---------------------------------------------------------------
library(gslnls)
ugslnls1<-try(gsl_nls(fn = frmu, data = weeddf,  start = stu1))
pnls(ugslnls1) # to get sum of squares
sgslnls1<-try(gsl_nls(fn = frms, data = weeddf,  start = sts1))
pnls(sgslnls1) # Use summary() to get display
tgslnls1<-try(gsl_nls(fn = frmt, data = weeddf,  start = stt1))
pnls(tgslnls1) 


## ----ex04singval, echo=TRUE, eval=TRUE----------------------------------------
# for nlsLM
if (inherits(tnlm1, "try-error"))  {
   print("Cannot compute solution -- likely singular Jacobian")
 } else {  
   JtnlsLM <- tnlm1$m$gradient() # actually the Jacobian
   svd(JtnlsLM)$d # Singular values
}   
# for gsl_nls
if (inherits(tgslnls1, "try-error")) {
   cat("Cannot compute solution -- likely singular Jacobian")
} else {  
   JtnlsLM <- tgslnls1$m$gradient()
   svd(JtnlsLM)$d # Singular values
}   


## ----ex05, echo=TRUE----------------------------------------------------------
stspecial<- c(Asym = 35.532,  xmid = 43376,  scal = -2935.4)
badstart<-nlxb(formula=frmt, start=stspecial, data=weeddf)
print(badstart)


## ----exhobbsfn, eval=TRUE-----------------------------------------------------
hobnlfb<-nlfb(start=stu1, resfn=hobbs.res, jacfn=hobbs.jac)
pshort(hobnlfb) # use print(hobnlfb) for more detail
hobnlm<-nls.lm(par=stu1, fn=hobbs.res, jac=hobbs.jac)
pnlslm(hobnlm)  
hobgsln<-gsl_nls(start=stu1, fn=hobbs.res, y=rep(0,12))
pnls(hobgsln)
hobgsl<-gsl_nls(start=stu1, fn=hobbs.res, y=rep(0,12), jac=hobbs.jac)
pnls(hobgsl) # using analytic Jacobian


## ----ex10, echo=TRUE----------------------------------------------------------
# Start MUST be feasible i.e. on or within bounds
anlshob1b <- nls(frms, start=sts1, data=weeddf, lower=c(0,0,0),
             upper=c(2,6,3), algorithm='port')
pnls(anlshob1b) #  check the answer (short form)
# nlsLM seems NOT to work with bounds in this example
anlsLM1b <- nlsLM(frms, start=sts1, data=weeddf, lower=c(0,0,0), upper=c(2,6,3))
pnls(anlsLM1b)
# also no warning if starting out of bounds, but gets a good answer!!
st4<-c(c1=4, c2=4, c3=4)
anlsLMob <- nlsLM(frms, start=st4, data=weeddf, lower=c(0,0,0), upper=c(2,6,3))
pnls(anlsLMob)
# Try nlsr::nlxb()
anlx1b <- nlxb(frms, start=sts1, data=weeddf, lower=c(0,0,0), upper=c(2,6,3))
pshort(anlx1b)


## ----ex10m, echo=TRUE---------------------------------------------------------
# Hobbsmaskx.R -- masks with formula specification of the problem
require(nlsr); require(minpack.lm); traceval<-FALSE
stu <- c(b1=200, b2=50, b3=0.3) # a default starting vector (named!)
sts <- c(c1=2, c2=5, c3=3) # a default scaled starting vector (named!)
# fix first parameter
anxbmsk1 <- try(nlxb(frmu, start=stu, data=weeddf, lower=c(200,0,0), 
			upper=c(200, 60, 3), trace=traceval))
print(anxbmsk1)
anlM1 <- try(nlsLM(frmu, start=stu, data=weeddf, lower=c(200,0,0), 
			upper=c(200, 60, 3), trace=traceval))
pnls(anlM1)
anlsmsk1 <- try(nls(frmu, start=stu, data=weeddf, lower=c(200,0,0), 
		upper=c(200, 60, 3),  algorithm="port", trace=traceval))
pnls(anlsmsk1)
# Hobbs scaled problem with bounds, formula specification
anlxmsks1 <- nlxb(frms, start=sts, data=weeddf, lower=c(2,0,0),
                  upper=c(2,6,30))
print(anlxmsks1)
anlshmsk1 <- nls(frms, start=sts, trace=traceval, data=weeddf, lower=c(2,0,0),
             upper=c(2,6,30), algorithm='port')
pnls(anlshmsk1)
anlsLMmsks1 <- nlsLM(frms, start=sts, data=weeddf, lower=c(2,0,0),
                 upper=c(2,6,30))
pnls(anlsLMmsks1)

# Test with all parameters masked
anlxmskall<- try(nlxb(frms, start=sts, data=weeddf, lower=sts, upper=sts))
print(anlxmskall)


## ----nlswtx, echo=TRUE--------------------------------------------------------
wts <- 0.5^tt # simple weights
frmlogis <- weed ~ Asym/(1 + exp((xmid - tt)/scal))
Asym<-1; xmid<-1; scal<-1
nowt<-nls(weed ~ SSlogis(tt, Asym, xmid, scal)) # UNWEIGHTED
rnowt<-nowt$m$resid() # This has UNWEIGHTED residual and Jacobian. Does NOT take coefficients.
attr(rnowt, "gradient")<-NULL; rnowt
usewt <- nls(weed ~ SSlogis(tt, Asym, xmid, scal), weights=wts)
rusewt<-usewt$m$resid() # WEIGHTED. Does NOT take coefficients.
attr(rusewt,"gradient")<-NULL; rusewt
source("nlsModel.R")
nmod0 <- nlsModel(frmlogis, data=weeddf, start=c(Asym=1, xmid=1, scal=1), wts=wts)
rn0<-nmod0$resid() # Parameters are supplied in nlsModel() `start` above.
attr(rn0,"gradient")<-NULL; rn0
nmod <- nlsModel(frmlogis, data=weeddf, start=coef(usewt), wts=wts)
rn<-nmod$resid()
attr(rn,"gradient")<-NULL; rn


## ----tetrarun-----------------------------------------------------------------
time <- c( 1,  2,  3,  4,  6 , 8, 10, 12, 16)
conc <- c( 0.7, 1.2, 1.4, 1.4, 1.1, 0.8, 0.6, 0.5, 0.3)
NLSdata <- data.frame(time,conc)
NLSstart <- c(lrc1 = -2, lrc2 = 0.25, A1 = 150, A2 = 50) # a starting vector (named!)
NLSformula <- conc ~ A1 * exp(-exp(lrc1) * time) + A2 * exp(-exp(lrc2) * time)
tryit <- try(nls(NLSformula, data = NLSdata, start = NLSstart, trace = TRUE))
print(tryit)


## ----log4ways, echo=TRUE, eval=FALSE------------------------------------------
#> DNase1 <- subset(DNase, Run == 1) # select the data
#> ## using a selfStart model - do not specify the starting parameters
#> fm1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
#> summary(fm1)
#> 
#> ## using conditional linearity - leave out the Asym parameter
#> fm2 <- nls(density ~ 1 / (1 + exp((xmid - log(conc)) / scal)),
#>                  data = DNase1, start = list(xmid = 0, scal = 1),
#>                  algorithm = "plinear")
#> summary(fm2)
#> 
#> ## without conditional linearity
#> fm3 <- nls(density ~ Asym / (1 + exp((xmid - log(conc)) / scal)),
#>                  data = DNase1,
#>                  start = list(Asym = 3, xmid = 0, scal = 1))
#> summary(fm3)
#> 
#> ## using Port's nl2sol algorithm
#> fm4 <- try(nls(density ~ Asym / (1 + exp((xmid - log(conc)) / scal)),
#>                  data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1),
#>                  algorithm = "port"))
#> summary(fm4)
#> 
#> ## using conditional linearity AND Asym does NOT work
#> fm2a <- try(nls(density ~ Asym / (1 + exp((xmid - log(conc)) / scal)),
#>                  data = DNase1, start = list(Asym=3, xmid = 0, scal = 1),
#>                  algorithm = "plinear", trace = TRUE))
#> summary(fm2a)

