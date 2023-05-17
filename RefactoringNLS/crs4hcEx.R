# crs4hcEx.R
library(microbenchmark)
library(crsnls)
x <- c(1,2,3,5,7,10)
y <- c(109,149,149,191,213,224)
df <- data.frame(x=x, y=y)
lowerBounds <- c(1, 0.1)
upperBounds <- c(1000, 2)
pnam <- c("b1", "b2")
names(lowerBounds)<-pnam
names(upperBounds)<-pnam
tmod<-microbenchmark(mod <- crs4hc(y ~ b1 * (1-exp(-b2*x)), df, lowerBounds, upperBounds))
mod
cat("evaluations =", mod$numOfEvals,"\n")
# str(mod)

tmde<-microbenchmark(mde <- crs4hce(y ~ b1 * (1-exp(-b2*x)), df, lowerBounds, upperBounds))
mde
cat("evaluations =", mde$numOfEvals,"\n")
# str(mde)
library(nlsr)
tmodn <- microbenchmark(modn <- nlxb(formula = y ~ b1 * (1-exp(-b2*x)), data=df, lower=lowerBounds, 
             upper=upperBounds, start=0.5*(lowerBounds+upperBounds)))
tmodn
modn
summary(tmod)
summary(tmde)
summary(tmodn)
mod$model$estimates
mde$model$estimates
coef(modn)
