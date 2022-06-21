# bhobbsX.R ## bounded formula specification of problem using Hobbs Weed problem
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
wf <- data.frame(y=weed, tt=tt)
st <- c(b1=1, b2=1, b3=1) # a default starting vector (named!)
wmods <- y ~ 100*b1/(1+10*b2*exp(-0.1*b3*tt)) ## Scaled model
require(minpack.lm)
require(nlsr)
# Hobbs scaled problem with bounds, formula specification
anlxb1 <- nlxb(wmods, start=st, data=wf, lower=c(0,0,0), upper=c(2,6,3))
cat("Using nlsr::nlxb():")
print(anlxb1)
## nlsLM seems NOT to work properly with bounds
anlsLM1b <- nlsLM(wmods, start=st, data=wf, lower=c(0,0,0), upper=c(2,6,3))
cat("\n"); cat("using minpack.lm::nlsLM():")
print(anlsLM1b)
anlsb<-nls(wmods, start=st, data=wf, algorithm="port", lower=c(0,0,0), upper=c(2,6,3))
print(anlsb)
