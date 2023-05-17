# Hobbsmaskx.R
## masks in Hobbs with formula specification of problem
require(nlsr)
require(minpack.lm)
stu <- c(b1=200, b2=50, b3=0.3) # a default starting vector (named!)
sts <- c(b1=2, b2=5, b3=3) # a default scaled starting vector (named!)
traceval<-FALSE

# fix first parameter
anxbmsk1 <- try(nlxb(frmu, start=stu, data=weeddf, lower=c(200,0,0), 
			upper=c(200, 60, .3), trace=traceval))
print(anxbmsk1)
anlM1 <- try(nlsLM(frmu, start=stu, data=weeddf, lower=c(200,0,0), 
			upper=c(200, 60, .3), trace=traceval))
print(anlM1)
# nls gives warnings
anlsmsk1 <- try(nls(frmu, start=stu, data=weeddf, lower=c(200,0,0), 
		upper=c(200, 60, .3),  algorithm="port", trace=traceval))
print(anlsmsk1)


# Hobbs scaled problem with bounds, formula specification
anlshmsk1 <- nls(frms, start=sts, trace=traceval, data=weeddf, lower=c(2,0,0),
             upper=c(2,6,3), algorithm='port')
print(anlshmsk1)
cat("More precisely...crossprod(resid(anlshmsk1))=",crossprod(resid(anlshmsk1)),"\n")
# nlsLM does not always work with bounds
anlsLMmsks1 <- nlsLM(frms, start=sts, data=weeddf, lower=c(2,0,0),
                 upper=c(2,6,3))
print(anlsLMmsks1)

anlxmsks1 <- nlxb(frms, start=sts, data=weeddf, lower=c(2,0,0),
                  upper=c(2,6,3))
print(anlxmsks1)

# Test with all parameters masked
anlxmskall<- try(nlxb(frms, start=sts, data=weeddf, lower=sts, upper=sts))
print(anlxmskall)
