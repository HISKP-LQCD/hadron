#  source(file="fit_fp.R")
#
#


## data
mu <- c(0.0040 , 0.0064 , 0.0100 )
fpsr0 <- c(0.33939648 , 0.3669825 , 0.38985338) 
fpsr0err <- c(0.001451 , 0.0018203 , 0.00102620  ) 

## now try and fit
dummy <- data.frame(x=mu , y = fpsr0)

## chiral log function

fpsLOG = function(f,slope,mu,lamb) {

tmp <- f * ( 1 - ( slope * mu * log(slope * mu / lamb)   ) / (4.0 * pi * f^2 )^2 )

return (tmp )
}


fpslinear = function(f,slope,mu,lamb) {

tmp <- f  +  slope * mu 

return (tmp )
}


#
# from my fits
#


fff <- function(mu) {

return (fpsLOG(0.308944,120.834,mu,0.4) )
}


out <- fff(mu)

##
## simple linear fit
##
##myFIT <- nls(fpsr0 ~ a*mu + b  ,  start = list(a=1 , b =0 ) , weights = 1/fpsr0err^2  , trace = TRUE ) 

# cut off in the chiral log term
lamb <- 0.8


###myFIT <- nls(fpsr0 ~ fpsLOG(f,slope,mu,lamb) ,  start = list(f=0.308 , slope =8 ) , weights = 1/fpsr0err^2  , trace = TRUE ) 

myFIT <- nls(fpsr0 ~ fpslinear(f,slope,mu,lamb) ,  start = list(f=0.308 , slope =8 ) , weights = 1/fpsr0err^2  , trace = TRUE ) 



summary(myFIT) 

fitted(myFIT)
fP <- coef(myFIT)[1]
print(fP)

slopeP <- coef(myFIT)[2]
print(slopeP)

##q()

##plot(myFIT)

##pp <- predict(myFIT)
##plot(pp)

#
#  plot the fit function
#

xx_start <- 0.0 
xx_end   <- 0.015
tot <- 100 
delta <- ( xx_end - xx_start ) / tot 

xx <- xx_start 

for ( i in 1:tot ) {

  yy <- fpslinear(fP,slopeP,xx,lamb)

  ss <- c(xx,yy)
  print(ss)
  xx <- xx + delta 
}
