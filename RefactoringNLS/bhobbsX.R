# bhobbsX.R
#### bounds with formula specification of problem
## Use the Hobbs Weed problem
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weedframe <- data.frame(y=weed, tt=tt)
st <- c(b1=1, b2=1, b3=1) # a default starting vector (named!)
## Unscaled model
wmodu <- y ~ b1/(1+b2*exp(-b3*tt))
## Scaled model
wmods <- y ~ 100*b1/(1+10*b2*exp(-0.1*b3*tt))
## We can provide the residual and Jacobian as functions
# Scaled Hobbs problem
shobbs.res  <-  function(x){ # scaled Hobbs weeds problem -- residual
  # This variant uses looping
  if(length(x) != 3) stop("shobbs.res -- parameter vector n!=3")
  y  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
           38.558, 50.156, 62.948, 75.995, 91.972)
  tt  <-  1:12
  res  <-  100.0*x[1]/(1+x[2]*10.*exp(-0.1*x[3]*tt)) - y
}

shobbs.jac  <-  function(x) { # scaled Hobbs weeds problem -- Jacobian
  jj  <-  matrix(0.0, 12, 3)
  tt  <-  1:12
  yy  <-  exp(-0.1*x[3]*tt)
  zz  <-  100.0/(1+10.*x[2]*yy)
  jj[tt,1]   <-   zz
  jj[tt,2]   <-   -0.1*x[1]*zz*zz*yy
  jj[tt,3]   <-   0.01*x[1]*zz*zz*yy*x[2]*tt
  attr(jj, "gradient") <- jj
  jj
}

##----- Hobbsbounded unique code starts here ------
require(minpack.lm)
traceval<-FALSE
cat("Infeasible start test\n")
start1inf <- c(b1=4, b2=4, b3=4) # b3 OUT OF BOUNDS for next few tries
## Infeasible start! No warning message!
anlM2i <- try(nlsLM(wmodu, start=start1inf, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                    trace=traceval))

# feasible start i.e. on or within bounds
start1<-st
anlM1 <- try(nlsLM(wmodu, start=start1, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                   trace=traceval))
print(anlM1)

# Hobbs scaled problem with bounds, formula specification

require(nlsr)
anlxb1 <- nlxb(wmods, start=start1, data=weedframe, lower=c(0,0,0),
               upper=c(2,6,3), trace=traceval)
print(anlxb1)

## nlsLM seems NOT to work properly with bounds
anlsLM1b <- nlsLM(wmods, start=start1, data=weedframe, lower=c(0,0,0),
                  upper=c(2,6,3), trace=traceval)
print(anlsLM1b)
newst<-coef(anlsLM1b)
print(newst)
anlsLM1bx <- nlsLM(wmods, start=newst, data=weedframe, lower=c(0,0,0),
                   upper=c(2,6,3))
print(anlsLM1bx)

cat("Single number bounds:\n")
anlsLM1b1 <- try(nlsLM(wmods, start=start1, data=weedframe, lower=0,
                       upper=3))
print(anlsLM1b1)
anlsLM1b1a <- try(nlsLM(wmods, start=start1, data=weedframe, lower=c(0,0,0),
                        upper=c(3,3,3)))
print(anlsLM1b1a)

anlxb1b1<- try(nlxb(wmods, start=start1, data=weedframe, lower=0,
                    upper=3))
print(anlxb1b1)

anlm1b <- nls.lm(par=start1, fn=shobbs.res, jac=shobbs.jac, lower=c(0,0,0),
                 upper=c(2,6,3))
print(anlm1b)


# functional presentation of problem  -- seems to work 
anlm1bx <- nls.lm(par=start1, fn=shobbs.res, jac=shobbs.jac, lower=c(0,0,0),
                  upper=c(2,6,3))
print(anlm1bx)

anlfb1bx <- nlfb(start=start1, resfn=shobbs.res, jacfn=shobbs.jac, lower=c(0,0,0),
                 upper=c(2,6,3))
print(anlfb1bx)
