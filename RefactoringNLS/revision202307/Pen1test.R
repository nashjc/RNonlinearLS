## Penalty function I [More, Garbow, Hillstrom 1981]
alpha <- 1e-5
a <- 1e-5
sqrta <- sqrt(alpha)
### residuals
fres <- function(x) {
  c(sqrt(alpha) * (x - 1), sum(x^2) - 0.25)
}
### jacobian
fjac <- function(x) {
  jj <- rbind(diag(sqrt(alpha), nrow = length(x)), 2 * t(x))
  attr(jj, "gradient") <- jj
  jj
}

ffn = function(par) {
  n <- length(par)
  if (n < 1) {
    stop("Penalty Function I: n must be positive")
  }
  fsum <- 0
  fn1 <- 0
  for (i in 1:n) {
    fi <- sqrta * (par[i] - 1)
    fsum <- fsum + fi * fi
    fn1 <- fn1 + par[i] * par[i]
  }
  
  fn1 <- fn1 - 0.25
  fsum <- fsum + fn1 * fn1
  fsum
}

fgr = function(par) {
  n <- length(par)
  if (n < 1) {
    stop("Penalty Function I: n must be positive")
  }
  grad <- rep(0, n)
  fn1 <- 0
  
  for (i in 1:n) {
    fi <- sqrta * (par[i] - 1)
    grad[i] <- grad[i] + 2 * sqrta * fi
    fn1 <- fn1 + par[i] * par[i]
  }
  fn1 <- fn1 - 0.25
  grad <- grad + 4 * par * fn1
  grad
}

n <- 10
# ftest
x0<-1:n
cat("ffn=",ffn(x0),"\n")
cat("fres:")
fr<-fres(x0)
print(fr)
cat("sum(fr^2)=", sum(fr^2), "\n")
cat("fgr:")
print(fgr(x0))
cat("Jacobian:\n")
JJ<-fjac(x0)
print(JJ)
cat("grad from Jacobian:")
grj <- 2*as.vector(t(JJ)%*%fr)
print(grj)

library(minpack.lm)
library(optimx)
n<-500
x0<-as.numeric(1:n)

tnlsr10<-bench::mark(nlsr10 = nlsr::nlfb(start = x0, resfn = fres, jacfn = fjac),  iterations = 10)
tnlsr10
tnlslm10<-bench::mark(nlslm10 = nls.lm(par = x0, fn = fres, jac = fjac),  iterations = 10)
tnlslm10
tcg10<-bench::mark(cg10 = optimr(par = x0, fn = ffn, gr = fgr, method="Rcgmin"),  iterations = 10)
tcg10
