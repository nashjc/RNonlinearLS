DNase1 <- subset(DNase, Run == 1) # select the data
## using a selfStart model - do not specify the starting parameters
fm1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
summary(fm1)

## using conditional linearity - leave out the Asym parameter
fm2 <- nls(density ~ 1/(1 + exp((xmid - log(conc))/scal)),
           data = DNase1, start = list(xmid = 0, scal = 1),
           algorithm = "plinear")
summary(fm2)
## without conditional linearity
fm3 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
           data = DNase1,
           start = list(Asym = 3, xmid = 0, scal = 1))
summary(fm3)

## using Port's nl2sol algorithm
fm4 <- try(nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
               data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1),
               algorithm = "port"))
summary(fm4)

## using conditional linearity AND Asym does NOT work
fm2a <- try(nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
                data = DNase1, start = list(Asym=3, xmid = 0, scal = 1),
                algorithm = "plinear", trace = TRUE))
summary(fm2a)

