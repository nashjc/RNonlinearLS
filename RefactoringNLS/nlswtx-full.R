weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
wts <- 0.5^tt # simple weights
frmlogis <- weed ~ Asym/(1 + exp((xmid - tt)/scal))
Asym<-1; xmid<-1; scal<-1
strt <- c(Asym=Asym, xmid=xmid, scal=scal)
modlogis<-model2rjfun(frmlogis, strt)
Asym<-1; xmid<-1; scal<-1
nowt<-nls(weed ~ SSlogis(tt, Asym, xmid, scal)) # UNWEIGHTED
nowt
nlsmod0<-nowt$m$resid() # This has UNWEIGHTED residual and Jacobian
nlsmod0
nlsrmod<-modlogis(coef(nowt)) # generate model using nlsr function
nlsrmod
inpar<-getInitial(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf)
print(inpar)
nowtx<-nlxb(frmlogis, start=c(Asym=1, xmid=1, scal=1), trace=FALSE)
nowtx
# nowtSSx0<-nlxb(weed ~ SSlogis(tt, Asym=1, xmid=1, scal=1), start=c(Asym=1, xmid=1, scal=1), trace=TRUE)
# nowtSsx0
# nowtSSx<-nlxb(weed ~ SSlogis(tt, Asym=1, xmid=1, scal=1), start=c(Asym=100, xmid=2, scal=1), trace=TRUE, 
#               control=list(japprox="jand"))
# nowtSSx
nowtxi<-nlxb(frmlogis, start=inpar, trace=FALSE)
nowtxi
usewt <- nls(weed ~ SSlogis(tt, Asym, xmid, scal), weights=wts)
usewt
usewt$m$resid() 
nmod <- nlsModel(frmlogis, data=weeddf, start=c(Asym=1, xmid=1, scal=1), wts=wts)
nmod$resid()
usewtx<-nlxb(frmlogis,  start=c(Asym=1, xmid=1, scal=1), weights=wts, trace=TRUE)
usewtx
# NOTE jafwd and jaback fail!
usewtxn<-nlxb(frmlogis,  start=c(Asym=1, xmid=1, scal=1), weights=wts, trace=TRUE, control=list(japprox="jacentral"))
usewtxn
usewtxi<-nlxb(frmlogis, start=inpar, weights=wts, trace=TRUE)
usewtxi
# usewtx$resid
