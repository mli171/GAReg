pkgname <- "GAReg"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "GAReg-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('GAReg')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("FDRCalc")
### * FDRCalc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: FDRCalc
### Title: False Discovery Rate (FDR) and True Positive Rate (TPR) from
###   index labels
### Aliases: FDRCalc

### ** Examples

# Simple example
N <- 10
true <- c(2, 4, 7)
pred <- c(4, 5, 7, 7) # duplicates are ignored
FDRCalc(true, pred, N)

# Empty predictions
FDRCalc(true, integer(0), N)

# All correct predictions
FDRCalc(true, true, N)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("FDRCalc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("crossover_fixknots")
### * crossover_fixknots

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: crossover_fixknots
### Title: Crossover Operator (Fixed-m) with Feasibility-First Restarts
### Aliases: crossover_fixknots

### ** Examples

## Not run: 
##D N <- 120
##D lmax <- 30
##D minDist <- 5
##D m <- 3
##D mom <- c(m, c(20, 50, 90), rep(0, lmax - 1 - m))
##D mom[m + 2] <- N + 1
##D dad <- c(m, c(18, 55, 85), rep(0, lmax - 1 - m))
##D dad[m + 2] <- N + 1
##D child <- crossover_fixknots(mom, dad, minDist = minDist, lmax = lmax, N = N)
##D child
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("crossover_fixknots", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fixknotsIC")
### * fixknotsIC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fixknotsIC
### Title: Information criterion for a fixed–knot spline regression
### Aliases: fixknotsIC

### ** Examples

library(MASS)
y <- mcycle$accel
x <- mcycle$times
x_unique <- sort(unique(x))
# chromosome encoding 5 interior knot indices with sentinel:
chrom <- c(5, 24, 30, 46, 49, 69, length(x_unique) + 1)
fixknotsIC(chrom,
  y = y, x = x, x_unique = x_unique,
  fixedknots = 5, ic_method = "BIC"
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fixknotsIC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gareg_knots")
### * gareg_knots

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gareg_knots
### Title: Genetic-Algorithm–based Optimal Knot Selection
### Aliases: gareg_knots

### ** Examples

## Not run: 
##D set.seed(1)
##D N <- 120
##D y <- c(rnorm(40, 0), rnorm(40, 3), rnorm(40, 0))
##D x <- seq_len(N)
##D 
##D # 1) Varying-knots with single-pop GA
##D g1 <- gareg_knots(
##D   y, x,
##D   minDist = 5,
##D   gaMethod = "cptga",
##D   cptgactrl = cptgaControl(popSize = 150, pcrossover = 0.9, maxgen = 500)
##D )
##D summary(g1)
##D 
##D # 2) Fixed knots (operators auto-injected unless overridden)
##D g2 <- gareg_knots(
##D   y, x,
##D   fixedknots = 5,
##D   minDist = 5
##D )
##D summary(g2)
##D 
##D # 3) Island GA with island-specific controls
##D g3 <- gareg_knots(
##D   y, x,
##D   gaMethod = "cptgaisl",
##D   minDist = 6,
##D   cptgactrl = cptgaControl(
##D     engine = "cptgaisl",
##D     numIslands = 8, maxMig = 250,
##D     popSize = 120, pcrossover = 0.9
##D   )
##D )
##D summary(g3)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gareg_knots", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gareg_subset")
### * gareg_subset

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gareg_subset
### Title: Genetic-Algorithm Best Subset Selection (GA / GAISL)
### Aliases: gareg_subset

### ** Examples

## Not run: 
##D if (requireNamespace("GA", quietly = TRUE)) {
##D   set.seed(1)
##D   n <- 100
##D   p <- 12
##D   X <- matrix(rnorm(n * p), n, p)
##D   y <- 1 + X[, 1] - 0.7 * X[, 4] + rnorm(n, sd = 0.5)
##D 
##D   # Default: subsetBIC (Gaussian – negative BIC), engine = GA::ga
##D   fit1 <- gareg_subset(y, X,
##D     gaMethod = "ga",
##D     gacontrol = list(popSize = 60, maxiter = 80, run = 40)
##D   )
##D   summary(fit1)
##D 
##D   # Island model: GA::gaisl
##D   fit2 <- gareg_subset(y, X,
##D     gaMethod = "gaisl",
##D     gacontrol = list(popSize = 40, maxiter = 60, islands = 4)
##D   )
##D   summary(fit2)
##D 
##D   # Logistic objective (subsetBIC handles GLM via ...):
##D   ybin <- rbinom(n, 1, plogis(0.3 + X[, 1] - 0.5 * X[, 2]))
##D   fit3 <- gareg_subset(ybin, X,
##D     gaMethod = "ga",
##D     family = stats::binomial(), # <- passed to subsetBIC via ...
##D     gacontrol = list(popSize = 60, maxiter = 80)
##D   )
##D   summary(fit3)
##D }
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gareg_subset", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("splineX")
### * splineX

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: splineX
### Title: Build spline design matrices (piecewise polynomials, natural
###   cubic, B-spline)
### Aliases: splineX

### ** Examples

set.seed(1)
x <- sort(rnorm(100))
k <- quantile(x, probs = c(.25, .5, .75))

# 1) Piecewise polynomials (degree 3)
X_pp <- splineX(x, knots = k, degree = 3, type = "ppolys", intercept = TRUE)
dim(X_pp) # n x ((3+1) + 3) = n x 7

# 2) Natural cubic spline (cubic, degree ignored)
X_ns <- splineX(x, knots = k, type = "ns", intercept = TRUE)

# 3) B-spline basis (degree 3)
X_bs <- splineX(x, knots = k, degree = 3, type = "bs", intercept = TRUE)

# Fit without a duplicated intercept:
# fit <- lm(y ~ 0 + X_pp)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("splineX", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("varyknotsIC")
### * varyknotsIC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: varyknotsIC
### Title: Information criterion for spline regression with a variable
###   number of knots
### Aliases: varyknotsIC

### ** Examples

## Example with 'mcycle' data (MASS)
# y <- mcycle$accel; x <- mcycle$times
# x_unique <- sort(unique(x))
# chrom <- c(5, 24, 30, 46, 49, 69, length(x_unique) + 1)
# varyknotsIC(chrom, y=y, x=x, x_unique=x_unique,
#             type="ppolys", degree=3, ic_method="BIC")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("varyknotsIC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
