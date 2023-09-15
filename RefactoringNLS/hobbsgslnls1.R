# Example 1: exponential model
# (https://www.gnu.org/software/gsl/doc/html/nls.html#exponential-fitting-example)

weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
# plot(weeddf, main="Hobbs weed infestation data")
# model formulas
frmu <- weed ~ b1/(1+b2*exp(-b3*tt))
frms <- weed ~ 100*c1/(1+10*c2*exp(-0.1*c3*tt))
frmt <- weed ~ Asym /(1 + exp((xmid-tt)/scal))
#
# Starting parameter sets
stu1<-c(b1=1, b2=1, b3=1)
sts1<-c(c1=1, c2=1, c3=1)
stt1<-c(Asym=1, xmid=1, scal=1)

h1fit <- gsl_nls(
  fn = frmu,                        ## model formula
  data = weeddf,                       ## model fit data
  start = stu1                       ## starting values
)
h1fit


## data
set.seed(1)
n <- 50
x <- (seq_len(n) - 1) * 3 / (n - 1)
f <- function(A, lam, b, x) A * exp(-lam * x) + b
y <- f(A = 5, lam = 1.5, b = 1, x) + rnorm(n, sd = 0.25)

## model fit
ex1_fit <- gsl_nls(
  fn = y ~ A * exp(-lam * x) + b,                        ## model formula
  data = data.frame(x = x, y = y),                       ## model fit data
  start = c(A = 0, lam = 0, b = 0)                       ## starting values
)
summary(ex1_fit)                                         ## model summary
predict(ex1_fit, interval = "prediction")                ## prediction intervals

## analytic Jacobian 1
gsl_nls(
  fn = y ~ A * exp(-lam * x) + b,                        ## model formula
  data = data.frame(x = x, y = y),                       ## model fit data
  start = c(A = 0, lam = 0, b = 0),                      ## starting values
  jac = function(par) with(as.list(par),                 ## jacobian
    cbind(A = exp(-lam * x), lam = -A * x * exp(-lam * x), b = 1)
  )
)

## analytic Jacobian 2
gsl_nls(
  fn = y ~ A * exp(-lam * x) + b,                        ## model formula
  data = data.frame(x = x, y = y),                       ## model fit data
  start = c(A = 0, lam = 0, b = 0),                      ## starting values
  jac = TRUE                                             ## automatic derivation
)

## self-starting model
gsl_nls(
  fn =  y ~ SSasymp(x, Asym, R0, lrc),                   ## model formula
  data = data.frame(x = x, y = y)                        ## model fit data
)

# Example 2: Gaussian function
# (https://www.gnu.org/software/gsl/doc/html/nls.html#geodesic-acceleration-example-2)

## data
set.seed(1)
n <- 300
x <- seq_len(n) / n
f <- function(a, b, c, x) a * exp(-(x - b)^2 / (2 * c^2))
y <- f(a = 5, b = 0.4, c = 0.15, x) * rnorm(n, mean = 1, sd = 0.1)

## Levenberg-Marquardt (default)
gsl_nls(
  fn = y ~ a * exp(-(x - b)^2 / (2 * c^2)),             ## model formula
  data = data.frame(x = x, y = y),                      ## model fit data
  start = c(a = 1, b = 0, c = 1),                       ## starting values
  trace = TRUE                                          ## verbose output
)

## Levenberg-Marquardt w/ geodesic acceleration 1
gsl_nls(
  fn = y ~ a * exp(-(x - b)^2 / (2 * c^2)),             ## model formula
  data = data.frame(x = x, y = y),                      ## model fit data
  start = c(a = 1, b = 0, c = 1),                       ## starting values
  algorithm = "lmaccel",                                ## algorithm
  trace = TRUE                                          ## verbose output
)

## Levenberg-Marquardt w/ geodesic acceleration 2
## second directional derivative
fvv <- function(par, v, x) {
  with(as.list(par), {
    zi <- (x - b) / c
    ei <- exp(-zi^2 / 2)
    2 * v[["a"]] * v[["b"]] * zi / c * ei + 2 * v[["a"]] * v[["c"]] * zi^2 / c * ei -
      v[["b"]]^2 * a / c^2 * (1 - zi^2) * ei -
      2 * v[["b"]] * v[["c"]] * a / c^2 * zi * (2 - zi^2) * ei -
      v[["c"]]^2 * a / c^2 * zi^2 * (3 - zi^2) * ei
  })
}

## analytic fvv 1
gsl_nls(
  fn = y ~ a * exp(-(x - b)^2 / (2 * c^2)),             ## model formula
  data = data.frame(x = x, y = y),                      ## model fit data
  start = c(a = 1, b = 0, c = 1),                       ## starting values
  algorithm = "lmaccel",                                ## algorithm
  trace = TRUE,                                         ## verbose output
  fvv = fvv,                                            ## analytic fvv
  x = x                                                 ## argument passed to fvv
)

## analytic fvv 2
gsl_nls(
  fn = y ~ a * exp(-(x - b)^2 / (2 * c^2)),             ## model formula
  data = data.frame(x = x, y = y),                      ## model fit data
  start = c(a = 1, b = 0, c = 1),                       ## starting values
  algorithm = "lmaccel",                                ## algorithm
  trace = TRUE,                                         ## verbose output
  fvv = TRUE                                            ## automatic derivation
)

# Example 3: Branin function
# (https://www.gnu.org/software/gsl/doc/html/nls.html#comparing-trs-methods-example)

## Branin model function
branin <- function(x) {
  a <- c(-5.1 / (4 * pi^2), 5 / pi, -6, 10, 1 / (8 * pi))
  f1 <- x[2] + a[1] * x[1]^2 + a[2] * x[1] + a[3]
  f2 <- sqrt(a[4] * (1 + (1 - a[5]) * cos(x[1])))
  c(f1, f2)
}

## Dogleg minimization w/ model as function
gsl_nls(
  fn = branin,                   ## model function
  y = c(0, 0),                   ## response vector
  start = c(x1 = 6, x2 = 14.5),  ## starting values
  algorithm = "dogleg"           ## algorithm
)
