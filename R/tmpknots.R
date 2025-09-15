
knots.obj <- function(knot.bin, plen=0, y, x_base=NULL, p=3L){

  n <- as.integer(length(y))
  m <- as.integer(knot.bin[1])
  ones <- rep(1, n)
  if (!is.null(x_base)) x_base <- as.matrix(x_base)

  if(m == 0L){
    x <- cbind(ones, x_base)
  }else{
    knot.vec <- knot.bin[2:(m+1)]
    x <- splines::bs(1:n, degree=p, knots = knot.vec)
    x <- cbind(ones, x_base, x)
  }

  ## Fastter fit: a thin wrapper to the "innermost" C code performing the
  ##              QR decomposition
  fit <- .lm.fit(x, y)
  SSRes <- sum(fit$residuals^2)
  BIC.val <- n*log(SSRes/n) + (p + 1L + m)*log(n)

  return(BIC.val)
}


knots_obj_fixknots <- function(knot.bin, plen=0, y, x_base=NULL, p=3L, fixedm){

  n <- as.integer(length(y))
  ones <- rep(1, n)
  if (!is.null(x_base)) x_base <- as.matrix(x_base)

  knot.vec <- knot.bin[2:(fixedm+1)]
  x <- splines::bs(1:n, degree=p, knots = knot.vec)
  x <- cbind(ones, x_base, x)

  ## Fastter fit: a thin wrapper to the "innermost" C code performing the
  ##              QR decomposition
  fit <- .lm.fit(x, y)
  SSRes <- sum(fit$residuals^2)
  BIC.val <- n*log(SSRes/n) + (p + 1L + fixedm)*log(n)

  return(BIC.val)
}

Popinitial_fixknots <- function(popSize, prange=NULL, N, minDist, Pb, mmax, lmax, fixedm){

  pop <- matrix(0, nrow=lmax, ncol=popSize)

  for(j in 1:popSize){
    pop[,j] = selectTau_uniform_exact(N, fixedm, minDist, lmax)
  }

  return(pop)
}

crossover_fixknots <- function(mom, dad, prange=NULL, minDist, lmax, N){

  output <- rep(0, lmax)

  m.child <- as.integer(dad[1])

  child <- dad[2:(m.child+1)]
  mom_only <- setdiff(mom[2:(m.child+1)], dad[2:(m.child+1)])

  if (length(mom_only) > 0L) {
    mom_only <- sample(mom_only, length(mom_only))
    for (v in mom_only) {
      if (runif(1) >= 0.5) next
      diffs <- abs(child - v)
      order_j <- sample(order(diffs),
                        size=sample(1:length(diffs), size=1) # number of swap from mom
                        )
      placed <- FALSE
      for (j in order_j) {
        cand <- child
        cand[j] <- v
        cand <- sort(cand)
        # cat("\n cand=", cand)
        if (!all(diff(cand) >= minDist)) next
        ok <- TRUE
        if (ok) { child <- cand; placed <- TRUE; break }
        # cat("\n final child=", child)
      }
    }
  }

  output[1] <- m.child
  output[2:(m.child+1)] <- child
  output[m.child+2] <- N+1

  return(output)
}

crossover_fixknots_2 <- function(mom, dad, prange=NULL, minDist, lmax, N){

  up_tol <- 30
  output <- rep(0, lmax)

  m_child <- as.integer(dad[1])
  child <- rep(NA, m_child)

  mom_tau <- mom[2:(m_child+1)]
  dad_tau <- dad[2:(m_child+1)]
  co_tab <- as.vector(rbind(mom_tau, dad_tau))

  i <- 1
  ii <- 0
  while (i <= m_child) {
    if(i == 1){
      # first pick random
      tmppick <- sample(co_tab[1:2], size=1)
      child[1] <- tmppick
    }else{
      tmp_co_tab <- co_tab[((i-1)*2+1):(i*2)]
      tmp_diff <- which(tmp_co_tab - tmppick > minDist)
      if(length(tmp_diff) > 1){
        tmppick <- sample(tmp_co_tab,1)
        child[i] <- tmppick
      }else if(length(tmp_diff) == 1){
        tmppick <- tmp_co_tab[tmp_diff[1]]
        child[i] <- tmppick
      }else{
        i <- 1
        child <- rep(NA, m_child)
        ii <- ii + 1
        next
      }
    }
    if(i==m_child){
      if(all(child - dad_tau == 0) | all(child - mom_tau == 0)){
        if(ii > up_tol){break}
        i <- 1
        child <- rep(NA, m_child)
      }else{
        i <- i + 1
      }
    }else{
      i <- i + 1
    }
    ii <- ii + 1
  }

  output[1] <- m_child
  output[2:(m_child+1)] <- child
  output[m_child+2] <- N+1

  return(output)
}

# mom <- c(3, 20, 40, 70, 101)
# dad <- c(3, 22, 42, 76, 101)
# crossover_fixknots_2(dad, mom, prange=NULL, minDist=5, lmax=51, N=100)
#
# dad <- c(3, 20, 25, 45, 101)
# mom <- c(3, 8, 23, 55, 101)
# crossover_fixknots_2(dad, mom, prange=NULL, minDist=5, lmax=51, N=100)
# mom <- c(3, 20, 40, 70, 101)
# dad <- c(3, 22, 42, 76, 101)
# crossover_fixknots_2(dad, mom, prange=NULL, minDist=5, lmax=51, N=100)
# dad <- c(5, 5, 15, 26, 40, 55, 101)
# mom <- c(5, 7, 18, 27, 44, 60, 101)
# crossover_fixknots_2(dad, mom, prange=NULL, minDist=5, lmax=51, N=100)
# mom <- c(3, 28, 45, 72, 101)
# dad <- c(3, 28, 41, 72, 101)
# crossover_fixknots_2(dad, mom, prange=NULL, minDist=5, lmax=51, N=100)

selectTau_uniform_exact <- function(N, m, minDist, lmax){

  output <- rep(0, lmax)

  if (N < (m + 1L) * minDist + 2L){stop("Infeasible: need N >= (m+1)*minDist + 2.")}
  L0 <- 1L + minDist # 6
  U0 <- N - m * minDist - 1L # (N - minDist - 1) - (m-1)*minDist
  S  <- (U0 - L0 + 1L) + m - 1L
  if (S < m){stop("Infeasible (S < m).")}
  picks <- sort(sample.int(S, m, replace = FALSE))
  h <- (L0 - 1L) + picks - (seq_len(m) - 1L)
  i <- h + (seq_len(m) - 1L) * minDist
  output[1] <- m
  output[2:(m+1)] <- i
  output[m+2] <- N+1

  return(output)
}

# N =100
# m = 5
# minDist=5
# lmax=51
#
# selectTau_uniform_exact(N, m, minDist, lmax)

mutation_fixknots <- function(child, p.range = NULL, minDist, Pb, lmax, mmax, N) {

  m <- child[1]
  childMut <- selectTau_uniform_exact(N, m, minDist, lmax)

  return(childMut)
}
















library(splines)

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
true.knots = c(40,50,70)
plot(x,y)



fitness.knots <- function(knot.bin, y){

  p = 3L

  knot.bin.ext = c(rep(0,p), knot.bin, rep(0,p))
  knot.vec = sort( which(knot.bin.ext != 0) )

  bs.fit = lm(y ~ bs(x, degree =3, knots = knot.vec ))
  -AIC(bs.fit, k=log(n)) ##Obj function, -BIC to be maximized

}

p=3
BICGA = GA::ga(type="binary", fitness = fitness.knots,
               nBits = (n-2*p), maxiter = 20000, run = 500,
               popSize = 200, monitor = T, y)
##Check results
BICsol = BICGA@solution
##[1, ] ##can be 2 sols since X[1] is free, take 1st sol
#BICsol[1] = 0
mBIC = sum(BICsol)
mBIC

##knot locations
BICsol = BICsol[-1]
tBIC = as.vector( which(BICsol == 1) )
tBIC_est = tBIC

bs.fit <- lm(y ~ bs(1:n, knots = sort(tBIC_est)))
plot(y, xlab = "x", ylab = "y",
     main=paste0("B-Splines Fitting with GA_AIC=", round(BIC(bs.fit),3)))
ht = seq(1, 100, length.out = 100)
lines(ht, predict(bs.fit, data.frame(x = ht)))
abline(v=tBIC_est, col="red")





library(changepointGA)

GA.res1 <- cptga(ObjFunc=knots.obj, N=n, minDist = 5, seed=1234, y=y)
summary(GA.res1)
tBIC1 <- GA.res1@overbestchrom[2:(1+GA.res1@overbestchrom[1])]

GA.res2 <- cptga(ObjFunc=knots_obj_fixknots,
                 N=n,
                 minDist=5,
                 popInitialize="Popinitial_fixknots",
                 selection="selection_linearrank",
                 crossover="crossover_fixknots",
                 mutation="mutation_fixknots",
                 y=y,
                 fixedm=3)
summary(GA.res2)
tBIC2 <- GA.res2@overbestchrom[2:(1+GA.res2@overbestchrom[1])]


GA.res3 <- cptga(ObjFunc=knots_obj_fixknots,
                 N=n,
                 minDist=5,
                 popInitialize="Popinitial_fixknots",
                 selection="selection_linearrank",
                 crossover="crossover_fixknots_2",
                 mutation="mutation_fixknots",
                 y=y,
                 fixedm=3)
summary(GA.res3)
tBIC3 <- GA.res3@overbestchrom[2:(1+GA.res3@overbestchrom[1])]



bs.fit <- lm(y ~ bs(1:n, knots = true.knots ))
plot(y, xlab = "x", ylab = "y",
     main=paste0("B-Splines Fitting with TRUE BIC=", round(BIC(bs.fit),3)))
ht = seq(1, 100, length.out = 100)
lines(ht, predict(bs.fit, data.frame(x = ht)))
abline(v=true.knots, col="blue")

bs.fit <- lm(y ~ bs(1:n, knots = sort(tBIC1)))
plot(y, xlab = "x", ylab = "y",
     main=paste0("B-Splines Fitting with GA BIC=", round(BIC(bs.fit),3)))
ht = seq(1, 100, length.out = 100)
lines(ht, predict(bs.fit, data.frame(x = ht)))
abline(v=tBIC1, col="red")

bs.fit <- lm(y ~ bs(1:n, knots = sort(tBIC2)))
plot(y, xlab = "x", ylab = "y",
     main=paste0("B-Splines Fitting with cptGA BIC + co1=", round(BIC(bs.fit),3)))
ht = seq(1, 100, length.out = 100)
lines(ht, predict(bs.fit, data.frame(x = ht)))
abline(v=tBIC2, col="red")

bs.fit <- lm(y ~ bs(1:n, knots = sort(tBIC3)))
plot(y, xlab = "x", ylab = "y",
     main=paste0("B-Splines Fitting with cptGA BIC + co2=", round(BIC(bs.fit),3)))
ht = seq(1, 100, length.out = 100)
lines(ht, predict(bs.fit, data.frame(x = ht)))
abline(v=tBIC3, col="red")
