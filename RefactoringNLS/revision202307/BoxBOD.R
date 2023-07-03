# BoxBOD.R -- NIST/ITL StRD
y<-c(109, 149, 149, 191, 213, 224)
x<-c(1,2,3,5,7,10)
boddta<-data.frame(x, y)
## Reference:     Box, G. P., W. G. Hunter, and J. S. Hunter (1978).
##               Statistics for Experimenters.  
##               New York, NY: Wiley, pp. 483-487.
# Formula               y = b1*(1-exp[-b2*x])  +  e
boxbodfrm<-   y ~ b1*(1-exp(-b2*x))
#      Starting values          Certified Values
#
#      Start 1  Start 2     Parameter     Standard Deviation
#  b1 =   1       100   2.1380940889E+02  1.2354515176E+01
#  b2 =   1       0.75  5.4723748542E-01  1.0455993237E-01

#Residual Sum of Squares:          1.1680088766E+03
#Residual Standard Deviation:      1.7088072423E+01
#Degrees of Freedom:                                4
#Number of Observations:                            6  

st1<-c(b1=1, b2=1)
st2<-c(b1=100, b2=0.75)
library(nlsr)
bbnlsr1<-nlsr::nlxb(formula=boxbodfrm, data=boddta, start=st1)
pshort(bbnlsr1)
bbnlsr2<-nlsr::nlxb(formula=boxbodfrm, data=boddta, start=st2)
pshort(bbnlsr2)
print(bbnlsr2)

# nls

bbnls1<-try(nls(formula=boxbodfrm, data=boddta, start=st1))
pnls(bbnls1)
bbnls2<-try(nls(formula=boxbodfrm, data=boddta, start=st2))
pnls(bbnls2)
bbnls1p<-try(nls(formula=boxbodfrm, data=boddta, start=st1, algorithm="port"))
pnls(bbnls1p)
bbnls2p<-try(nls(formula=boxbodfrm, data=boddta, start=st2, algorithm="port"))
pnls(bbnls2p)

library(minpack.lm)
bbnlsLM1<-nlsLM(formula=boxbodfrm, data=boddta, start=st1)
pnlslm(bbnlsLM1)
bbnlsLM2<-nlsLM(formula=boxbodfrm, data=boddta, start=st2)
pnlslm(bbnlsLM2)









