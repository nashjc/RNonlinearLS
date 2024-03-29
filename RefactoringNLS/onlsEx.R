library(onls)
## 1A. The DNase data from 'nls',
## use all generic functions.
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod1 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
print(mod1)
plot(mod1)
summary(mod1)
predict(mod1, newdata = data.frame(conc = 6))
logLik(mod1)
deviance(mod1)
formula(mod1)
weights(mod1)
df.residual(mod1)
fitted(mod1)
residuals(mod1)
vcov(mod1)
coef(mod1)

DNase2 <- DNase1
DNase2$conc <- DNase2$conc * 2
mod2 <- update(mod1, data = DNase2)
print(mod2)


## 1B. Same model as above, but using the restricted
## predictor range which results in non-orthogonality
## of some points.
onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
     data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1),
     extend = c(0, 0))

## 2. Example from odrpack_guide.pdf, 2.C.i, pages 39ff.
x <- c(0, 0, 5, 7, 7.5, 10, 16, 26, 30, 34, 34.5, 100)
y <- c(1265, 1263.6, 1258, 1254, 1253, 1249.8, 1237, 1218, 1220.6, 
       1213.8, 1215.5, 1212)
DAT <- data.frame(x, y)
mod3 <- onls(y ~ b1 + b2 * (exp(b3 * x) -1)^2, data = DAT, 
             start = list(b1 = 1500, b2 = -50, b3 = -0.1))
deviance_o(mod3) # 21.445 as in page 47
coef(mod3) # 1264.65481/-54.01838/-0.08785 as in page 48

## 3. Example from Algorithm 676: ODRPACK, page 355 + 356.
x <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90, 95, 100, 105)
y <- c(4.14, 8.52, 16.31, 32.18, 64.62, 98.76, 151.13, 224.74, 341.35, 
       423.36, 522.78, 674.32, 782.04, 920.01)
DAT <- data.frame(x, y)
mod4 <- onls(y ~ b1 * 10^(b2 * x/(b3 + x)), data = DAT, 
             start = list(b1 = 1, b2 = 5, b3 = 100))
coef(mod4) # 4.4879/7.1882/221.8383 as in page 363
deviance_o(mod4) # 15.263 as in page 363

## 4. Example with bounds from simple_example.f90
## in https://www.netlib.org/toms/869.zip.
x <- c(0.982, 1.998, 4.978, 6.01)
y <- c(2.7, 7.4, 148.0, 403.0)
DAT <- data.frame(x, y)
mod5 <- onls(y ~ b1 * exp(b2 * x), data = DAT, 
             start = list(b1 = 2, b2 = 0.5), 
             lower = c(0, 0), upper = c(10, 0.9))
coef(mod5) # 1.4397(1.6334)/0.9(0.9) ## Different to reference!
deviance_o(mod5) # 0.1919 (0.2674) => but lower RSS than original ODRPACK!

## 5. Example with a fixed parameter
## => Asym = 3.
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod6 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1), 
             fixed = c(TRUE, FALSE, FALSE))
print(mod6)

## 6. Example to show that one can even conduct
## linear orthogonal regression (Deming regression).
## Comparison to XLstat
## https://help.xlstat.com/6650-run-deming-regression-compare-methods-excel
x <- c(9.8, 9.7, 10.7, 10.9, 12.4, 12.5, 12.8, 12.8, 12.9, 13.3, 
       13.4, 13.5, 13.7, 14.9, 15.2, 15.5)
y <- c(10.1, 11.4, 10.8, 11.3, 11.8, 12.1, 12.3, 13.6, 14.2, 14.4,
       14.6, 15.3, 15.5, 15.8, 16.2, 16.5)
DAT <- data.frame(x, y)
mod7 <- onls(y ~ a + b * x, data = DAT, 
             start = list(a = 2, b = 3))
print(mod7) ## -1.909/1.208 as on webpage
plot(mod7)
