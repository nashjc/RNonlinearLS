library(nlsr)
## Penalty function I [More, Garbow, Hillstrom 1981]
alpha <- 1e-5
n <- 1000
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
bench::mark(
  "large problem" = nlsr::nlfb(start = 1:n, resfn = fres, jacfn = fjac, trace=TRUE, control=list(femax=1000)),
  iterations = 1
)
