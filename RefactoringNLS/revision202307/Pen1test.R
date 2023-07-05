# Penalty function I [More, Garbow, Hillstrom 1981]
# sink("Pen1Times.txt", split=TRUE)
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
n<-10
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

library(nlsr)
library(minpack.lm)
library(optimx)

nvals <- c(50, 100, 150, 200, 250, 300)
#  nvals <- c(2, 10, 25, 50)
nn<-length(nvals)
TT <- matrix(NA, nrow=nn, ncol=4)
MM <- TT
for (i in 1:nn){
n <- nvals[i]
TT[i, 1] <- n
MM[i, 1] <- n
x0<-as.numeric(1:n)
tnlsr10<-bench::mark(nlsr10 <- nlsr::nlfb(start = x0, resfn = fres, jacfn = fjac),  iterations = 10, time_unit='ms')
pshort(nlsr10)
TT[i,2] <- tnlsr10$median
MM[i,2] <- tnlsr10$min
tnlslm10<-bench::mark(nlslm10 <- nls.lm(par = x0, fn = fres, jac = fjac),  iterations = 10, time_unit='ms')
pnlslm(nlslm10)
TT[i,3] <- tnlslm10$median
MM[i,3] <- tnlslm10$min
tcg10<-bench::mark(cg10 <- optimr(par = x0, fn = ffn, gr = fgr, method="Rcgmin"),  iterations = 10, time_unit='ms')
proptimr(cg10)
TT[i,4] <- tcg10$median
MM[i,4] <- tcg10$min
}
colnames(TT)<-c("n","nlfb","nls.lm","Rcgmin")
colnames(MM)<-c("n","nlfb","nls.lm","Rcgmin")
print(TT)
print(MM)
sink()
