#parameters used to generate the data
# rm(list=ls())
# https://stats.stackexchange.com/questions/13053/singular-gradient-error-in-nls-with-correct-starting-values
reala=-3
realb=5
realc=0.5
realr=0.7
realm=1
x=1:11 #x values - I have 11 timepoint data
#linear+exponential function
y=reala + realb*realr^(x-realm) + realc*x
#add a bit of noise to avoid zero-residual data
jitter_y = jitter(y,amount=0.2)
testdat=data.frame(x,jitter_y)
strt<-list(a=-3, b=5, c=0.5, r=0.7, m=1)
#try the regression with similar starting values to the the real parameters
jform<-jitter_y~a+b*r^(x-m)+c*x
linexp=try(nls(jform, data=testdat, start=strt, trace=F))
library(nlsr)
jmod<-model2rjfun(jform, strt, jacobian=TRUE, gradient=TRUE)
jmodval<-jmod(strt)
jmodval
jmodJ<-attr(jmodval, "gradient")
cat("Singular values of Jacobian:")
print(svd(jmodJ)$d)
# We can also get the singular values of the Jacobian from nlsr::nlxb()
# BUT note that the values are printed as a column purely to save space
linexp2<-try(nlxb(jform, data=testdat, start=strt, trace=F))
linexp2
