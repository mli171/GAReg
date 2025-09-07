##Author: Xueheng Shi

library(splines)
library(GA)

##Generate Data for Testing
n = 100
x = 1:n
beta0 = numeric(100)
beta0[1:40] = (1:40-20)^3
beta0[40:50] = -60*(40:50-50)^2 + 60*100+20^3
beta0[50:70] = -20*(50:70-50)^2 + 60*100+20^3
beta0[70:100] = -1/6*(70:100-110)^3 + -1/6*40^3 + 6000
beta0 = -beta0
beta0 = (beta0-min(beta0))*10/diff(range(beta0))
y = beta0 + rnorm(n)
plot(x,y)


##bs(y, knots = c(20,35,70))

p=5
#knot.bin = rep(0, n-10)
#knot.bin[c(23,56)]=1


fitness.knots = function(knot.bin){
  ##First and last p are not defined to be knots

  knot.bin.ext = c(rep(0,p), knot.bin, rep(0,p))
  knot.vec = sort( which(knot.bin.ext != 0) )

  bs.fit = lm(y ~ bs(x, degree =3, knots = knot.vec ))
  -AIC(bs.fit, k=log(n)) ##Obj function, -BIC to be maximized
}

BICGA = GA::ga(type="binary", fitness = fitness.knots,
               nBits = (n-2*p), maxiter = 20000, run = 5000,
               popSize = 200, monitor = T)

library(changepointGA)


fitness.knots = function(knot.bin){
  ##First and last p are not defined to be knots

  knot.bin.ext = c(rep(0,p), knot.bin, rep(0,p))
  knot.vec = sort( which(knot.bin.ext != 0) )

  bs.fit = lm(y ~ bs(x, degree =3, knots = knot.vec ))
  -AIC(bs.fit, k=log(n)) ##Obj function, -BIC to be maximized
}


knot.bin = c(2, 23, 56, n+1)

fitness.knots.1 = function(knot.bin, plen=0){

  n = tail(knot.bin, 1)
  m = knot.bin[1]
  knot.vec = knot.bin[2:(m+1)]

  bs.fit = lm(y ~ bs(x, degree =3, knots = knot.vec ))
  AIC(bs.fit, k=log(n)) ##Obj function, -BIC to be maximized
}

GA.res = changepointGA::cptga(ObjFunc=fitness.knots.1, N=n, minDist = 5)
summary(GA.res)

GA.res@overbestchrom # need add 5

##Check results
BICsol = BICGA@solution
##[1, ] ##can be 2 sols since X[1] is free, take 1st sol
#BICsol[1] = 0
mBIC = sum(BICsol)
mBIC
#5

##knot locations
tBIC = as.vector( which(BICsol == 1) )
tBIC+5



bs.fit = lm(y ~ bs(x, knots = sort(tBIC+5) ))
#summary(fm1)
#AIC(bs.fit, k=log(n))

plot(y, xlab = "x", ylab = "y", main="Test B-Splines Fitting")
ht = seq(1, 100, length.out = 100)
lines(ht, predict(bs.fit, data.frame(x = ht)))


