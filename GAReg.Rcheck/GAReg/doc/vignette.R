## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

set.seed(12345)

## ----setup--------------------------------------------------------------------
library(GAReg)

## -----------------------------------------------------------------------------
sim_subset_data <- function(n = 60, p = 50, s0 = 25, sigma = 1.5,
                            magnitudes_range = c(0.5, 2),
                            rho = NULL,
                            seed = NULL) {
  stopifnot(n > 0, p > 0, s0 >= 0, s0 <= p, sigma >= 0)
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)

  # Active set and coefficients
  true_idx <- if (s0 > 0) sort(sample.int(p, s0)) else integer(0)
  signs <- if (s0 > 0) sample(c(-1, 1), s0, replace = TRUE) else numeric(0)
  magnitudes <- if (s0 > 0) runif(s0, magnitudes_range[1], magnitudes_range[2]) else numeric(0)

  beta_true <- numeric(p)
  if (s0 > 0) beta_true[true_idx] <- magnitudes * signs

  if (is.null(rho)) {
    e <- rnorm(n, sd = sigma)
  } else {
    sd_innov <- sigma * sqrt(1 - rho^2)
    burn_in <- 100
    z <- rnorm(n + burn_in, sd = sd_innov)
    e_full <- numeric(n + burn_in)
    for (t in 2:(n + burn_in)) e_full[t] <- rho * e_full[t - 1] + z[t]
    e <- e_full[(burn_in + 1):(burn_in + n)]
  }

  y <- as.numeric(X %*% beta_true + e)

  DF <- data.frame(y = y, as.data.frame(X))
  colnames(DF)[-1] <- paste0("X", seq_len(p))

  list(
    X = X,
    y = y,
    beta_true = beta_true,
    true_idx = true_idx,
    DF = DF,
    rho = if (is.null(rho)) NULL else rho,
    args = list(
      n = n, p = p, s0 = s0, sigma = sigma,
      magnitudes_range = magnitudes_range,
      rho = rho, seed = seed
    )
  )
}

## -----------------------------------------------------------------------------
sim <- sim_subset_data(n = 100, p = 50, s0 = 25, sigma = 1.5, rho = NULL, seed = 123)
y <- sim$y
X <- sim$X
ga <- gareg_subset(
  y = y, X = X, gaMethod = "GA", monitor = FALSE,
  gacontrol = list(
    popSize = 120,
    maxiter = 20000,
    run = 4000,
    pmutation = 0.02
  )
)
summary(ga)

res <- FDRCalc(truelabel = sim$true_idx, predlabel = ga@bestsol, N = 50)
# FALSE Discover Rate
res$fdr
# TRUE Positive Rate
res$tpr

## -----------------------------------------------------------------------------
BIC.cpt <- function(chromosome, Xt) {
  m <- chromosome[1]
  tau <- chromosome[2:(2 + m - 1)]
  N <- length(Xt)

  if (m == 0) {
    mu.hat <- mean(Xt)
    sigma.hatsq <- sum((Xt - mu.hat)^2) / N
    BIC.obj <- 0.5 * N * log(sigma.hatsq) + 2 * log(N)
  } else {
    tau.ext <- c(1, tau, N + 1)
    seg.len <- diff(tau.ext)
    ff <- rep(0:m, times = seg.len)
    Xseg <- split(Xt, ff)
    mu.seg <- unlist(lapply(Xseg, mean), use.names = F)
    mu.hat <- rep(mu.seg, seg.len)
    sigma.hatsq <- sum((Xt - mu.hat)^2) / N
    BIC.obj <- 0.5 * N * log(sigma.hatsq) + (2 * m + 2) * log(N)
  }
  return(BIC.obj)
}

# IID data
set.seed(1234)
n <- 200
et <- rnorm(n)
Xt <- et + rep(c(0, 2, 0, 2), each = n / 4)

library(changepointGA)
GA.res <- cptga(
  ObjFunc = BIC.cpt, N = n, popSize = 500,
  pcrossover = 0.95, pmutation = 0.3, pchangepoint = 10 / n,
  Xt = Xt
)
summary(GA.res)

## -----------------------------------------------------------------------------
library(MASS)
library(splines)

data(mcycle)
head(mcycle)

## -----------------------------------------------------------------------------
g1 <- gareg_knots(
  y = mcycle$accel, x = mcycle$times,
  minDist = 5,
  gaMethod = "cptga",
  cptgactrl = cptgaControl(popSize = 200, pcrossover = 0.9, pmutation = 0.3),
  ic_method = "BIC"
)
summary(g1)

# knots location
g1@bestsol

## -----------------------------------------------------------------------------
g2 <- gareg_knots(
  y = mcycle$accel, x = mcycle$times,
  minDist = 5,
  gaMethod = "cptgaisl",
  cptgactrl = cptgaControl(
    numIslands = 5, popSize = 200, maxMig = 250,
    pcrossover = 0.9, pmutation = 0.3
  ),
  ic_method = "BIC"
)
summary(g2)

## -----------------------------------------------------------------------------
g3 <- gareg_knots(
  y = mcycle$accel, x = mcycle$times,
  fixedknots = 3,
  minDist = 5,
  gaMethod = "cptga",
  cptgactrl = cptgaControl(popSize = 200, pcrossover = 0.9, pmutation = 0.3),
  ic_method = "BIC"
)
summary(g3)

## -----------------------------------------------------------------------------
g4 <- gareg_knots(
  y = mcycle$accel, x = mcycle$times,
  fixedknots = 4,
  minDist = 5,
  gaMethod = "cptgaisl",
  cptgactrl = cptgaControl(
    numIslands = 5, popSize = 200, maxMig = 250,
    pcrossover = 0.9, pmutation = 0.3
  ),
  ic_method = "BIC"
)
summary(g4)

## ----fig.width=7, fig.height=5, out.width="90%", fig.align='center'-----------
y <- mcycle$accel
x <- mcycle$times
x_unique <- unique(x)

tBIC.vary.ga <- g1@bestsol
tBIC.vary.gaisl <- g2@bestsol
tBIC.fix.3.ga <- g3@bestsol
tBIC.fix.4.gaisl <- g4@bestsol

bsfit.vary.ga <- lm(y ~ bs(x, knots = x_unique[g1@bestsol], Boundary.knots = range(x)))
bsfit.vary.gaisl <- lm(y ~ bs(x, knots = x_unique[g2@bestsol], Boundary.knots = range(x)))
bsfit.fix.3.ga <- lm(y ~ bs(x, knots = x_unique[g3@bestsol], Boundary.knots = range(x)))
bsfit.fix.4.gaisl <- lm(y ~ bs(x, knots = x_unique[g4@bestsol], Boundary.knots = range(x)))

plot(x, y, xlab = "Time (ms)", ylab = "Acceleration (g)")
ht <- seq(min(x), max(x), length.out = 200)
lines(ht, predict(bsfit.vary.ga, data.frame(x = ht)), col = "blue", lty = 5, lwd = 2)
lines(ht, predict(bsfit.vary.gaisl, data.frame(x = ht)), col = "orange", lty = 4, lwd = 2)
lines(ht, predict(bsfit.fix.3.ga, data.frame(x = ht)), col = "purple", lty = 3, lwd = 2)
lines(ht, predict(bsfit.fix.4.gaisl, data.frame(x = ht)), col = "#D55E00", lty = 2, lwd = 2)
legend("bottomright",
  legend = c(
    "Varying knots GA",
    "Varying knots island model GA",
    "Fixed 3 knots GA",
    "Fixed 4 knots island model GA"
  ),
  lty = 5:2, lwd = 2,
  col = c("blue", "orange", "purple", "#D55E00"), bty = "n"
)

## ----vary-ppolys--------------------------------------------------------------
g_pp3 <- gareg_knots(
  y = y, x = x,
  minDist = 5,
  gaMethod = "cptga",
  ObjFunc = NULL, # use default varyknotsIC
  type = "ppolys",
  degree = 3, # degree-3 piecewise polynomial
  intercept = TRUE,
  ic_method = "BIC"
)
summary(g_pp3)

## ----vary-ns------------------------------------------------------------------
g_ns <- gareg_knots(
  y = y, x = x,
  minDist = 5,
  gaMethod = "cptga",
  type = "ns", # natural cubic (degree ignored)
  degree = 3, # ignored for "ns"
  intercept = TRUE,
  ic_method = "BIC"
)
summary(g_ns)

## ----vary-bs------------------------------------------------------------------
g_bs1 <- gareg_knots(
  y = y, x = x,
  minDist = 5,
  gaMethod = "cptga",
  type = "bs",
  degree = 1, # linear B-splines
  intercept = TRUE,
  ic_method = "BIC"
)
summary(g_bs1)

## -----------------------------------------------------------------------------
sessionInfo()

