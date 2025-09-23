#' Genetic-Algorithm Best Subset Selection (GA / GAISL)
#'
#' @description
#' Runs a GA-based search over variable subsets using a user-specified
#' objective (default: \code{\link{subsetBIC}}) and returns a compact
#' \code{"gareg"} S4 result with \code{method = "subset"}.
#' The engine can be \code{\link[GA]{ga}} (single population) or
#' \code{\link[GA]{gaisl}} (islands), selected via \code{gaMethod}.
#'
#' @param y Numeric response vector (length \code{n}).
#' @param X Numeric matrix of candidate predictors (\code{n} rows by \code{p} columns).
#' @param ObjFunc Objective function or its name. Defaults to \code{\link{subsetBIC}}.
#'   The objective must accept as its first argument a binary chromosome
#'   (0/1 mask of length \code{p}) and may accept additional arguments passed via \code{...}.
#'   By convention, \code{subsetBIC} returns \emph{negative} BIC, so the GA maximizes fitness.
#' @param gaMethod GA backend to call: \code{"ga"} or \code{"gaisl"} (functions from package \pkg{GA}),
#'   or a GA-compatible function with the same interface as \code{\link[GA]{ga}}.
#' @param gacontrol Optional named list of GA engine controls (e.g., \code{popSize}, \code{maxiter},
#'   \code{run}, \code{pcrossover}, \code{pmutation}, \code{elitism}, \code{seed}, \code{parallel},
#'   \code{keepBest}, \code{monitor}, ...). These are passed to the GA engine, not to the objective.
#' @param monitoring Logical; if \code{TRUE}, prints a short message and (if not supplied in
#'   \code{gacontrol}) sets \code{monitor = GA::gaMonitor} for live progress.
#' @param seed Optional RNG seed (convenience alias for \code{gacontrol$seed}).
#' @param ... Additional arguments forwarded to \code{ObjFunc} (not to the GA engine). For
#'   \code{\link{subsetBIC}} these typically include \code{family}, \code{weights}, \code{offset},
#'   and \code{control}.
#'
#' @details
#' The fitness passed to \pkg{GA} is \code{ObjFunc} itself. Because the engine expects
#' a function with signature \code{f(chrom, ...)}, your \code{ObjFunc} must interpret
#' \code{chrom} as a 0/1 mask over the columns of \code{X}. The function then computes a score
#' (e.g., negative BIC) using \code{y}, \code{X}, and any extra arguments supplied via \code{...}.
#'
#' With the default \code{\link{subsetBIC}}, the returned value is \code{-BIC}, so we set
#' \code{max = TRUE} in the GA call to maximize fitness. If you switch to an objective that
#' returns a quantity to \emph{minimize}, either negate it in your objective or change
#' the engine setting to \code{max = FALSE}.
#'
#' Engine controls belong in \code{gacontrol}; objective-specific options belong in \code{...}.
#' This separation prevents accidental name collisions between GA engine parameters and
#' objective arguments.
#'
#' @return An object of S4 class \code{"gareg"} (with \code{method = "subset"}) containing:
#' \itemize{
#'   \item \code{call} – the matched call.
#'   \item \code{N} – number of observations.
#'   \item \code{objFunc} – the objective function used.
#'   \item \code{gaMethod} – \code{"ga"} or \code{"gaisl"}.
#'   \item \code{gaFit} – the GA fit object returned by \pkg{GA} (if your class allows it).
#'   \item \code{featureNames} – column names of \code{X} (or empty).
#'   \item \code{bestFitness} – best fitness value (\code{GA::ga}@fitnessValue).
#'   \item \code{bestChrom} – \code{c(m, idx)}: number of selected variables and their indices.
#'   \item \code{bestnumbsol} – \code{m}, number of selected variables.
#'   \item \code{bestsol} – vector of selected column indices in \code{X}.
#' }
#'
#' @seealso
#' \code{\link{subsetBIC}},
#' \code{\link[GA]{ga}},
#' \code{\link[GA]{gaisl}}
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("GA", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 100; p <- 12
#'   X <- matrix(rnorm(n*p), n, p)
#'   y <- 1 + X[,1] - 0.7*X[,4] + rnorm(n, sd = 0.5)
#'
#'   # Default: subsetBIC (Gaussian – negative BIC), engine = GA::ga
#'   fit1 <- gareg_subset(y, X, gaMethod = "ga",
#'                        gacontrol = list(popSize = 60, maxiter = 80, run = 40))
#'   summary(fit1)
#'
#'   # Island model: GA::gaisl
#'   fit2 <- gareg_subset(y, X, gaMethod = "gaisl",
#'                        gacontrol = list(popSize = 40, maxiter = 60, islands = 4))
#'   summary(fit2)
#'
#'   # Logistic objective (subsetBIC handles GLM via ...):
#'   ybin <- rbinom(n, 1, plogis(0.3 + X[,1] - 0.5*X[,2]))
#'   fit3 <- gareg_subset(ybin, X, gaMethod = "ga",
#'                        family = stats::binomial(),           # <- passed to subsetBIC via ...
#'                        gacontrol = list(popSize = 60, maxiter = 80))
#'   summary(fit3)
#' }
#' }
#'
#' @export
gareg_subset = function(y,
                        X,
                        ObjFunc = NULL,
                        gaMethod = "ga",
                        gacontrol = NULL,
                        monitoring = FALSE,
                        seed = NULL,
                        ...){

  if (!requireNamespace("GA", quietly = TRUE))
    stop("Package 'GA' is required: install.packages('GA')")

  call <- match.call()
  X    <- as.matrix(X)
  nobs <- NROW(X)
  p    <- NCOL(X)

  ga_fun <- if (is.function(gaMethod)) {
    gaMethod
  } else {
    switch(tolower(as.character(gaMethod)),
           "ga"    = GA::ga,
           "gaisl" = GA::gaisl,
           stop("gaMethod must be 'ga' or 'gaisl' (from the GA package), or a GA-compatible function.")
    )
  }
  ga_name <- if (identical(ga_fun, GA::gaisl)) "gaisl" else "ga"

  if (is.null(ObjFunc)) {
    ObjFunc <- subsetBIC
  } else if (is.character(ObjFunc)) {
    ObjFunc <- get(ObjFunc, mode = "function")
  }
  if (!is.function(ObjFunc)) stop("ObjFunc must be a function or a function name.", call. = FALSE)

  engine_args <- as.list(gacontrol %||% list())
  if (!is.null(seed)) engine_args$seed <- seed
  if (isTRUE(monitoring) && is.null(engine_args$monitor)) engine_args$monitor <- GA::gaMonitor

  if (isTRUE(monitoring)) message("Running best subset via GA (engine = ", ga_name, ")")

  base_args <- list(
    fitness = ObjFunc,   # your subsetBIC; expects first arg = chromosome (0/1 mask)
    type    = "binary",
    nBits   = p,
    monitor = monitoring
  )

  GA.res <- do.call(ga_fun,
                    c(base_args,
                      engine_args,
                      list(y = y, X = X),
                      list(...)))

  sol <- GA.res@solution
  if (is.matrix(sol)) sol <- sol[1L, ]
  idx      <- which(sol != 0)
  mhat     <- length(idx)
  bestChrom <- c(mhat, idx)

  `%||%` <- function(a, b) if (!is.null(a)) a else b
  feat_names <- colnames(X) %||% character(p)

  # Build gareg S4 (gaFit kept NULL unless your class accepts GA objects)
  object <- methods::new("gareg",
                         call         = call,
                         method       = "subset",
                         N            = nobs,
                         objFunc      = ObjFunc,
                         gaMethod     = ga_name,
                         gaFit        = GA.res,
                         featureNames = feat_names,
                         bestFitness  = as.numeric(GA.res@fitnessValue),
                         bestChrom    = as.numeric(bestChrom),
                         bestnumbsol  = as.numeric(mhat),
                         bestsol      = as.numeric(idx)
  )

  return(object)

}

#' Unified BIC-style Objective for Subset Selection (GLM & Gaussian)
#'
#' @description
#' Computes a BIC-like criterion for a chromosome that encodes a variable
#' subset. The same expression
#' \deqn{\mathrm{BIC} = n \log(\mathrm{rss\_like}/n) + k \log n}
#' is used for all families, where:
#' \itemize{
#'   \item For Gaussian with identity link, \code{rss_like} is the residual sum of squares (RSS),
#'         computed via a fast \code{.lm.fit}.
#'   \item For other GLM families, \code{rss_like} is the residual \emph{deviance}
#'         from \code{glm.fit}.
#' }
#' The effective parameter count \eqn{k} includes the intercept.
#'
#' @details
#' The chromosome \code{subset_bin} is a \emph{binary} vector (0/1 by column),
#' indicating which predictors from \code{X} are included. The design matrix
#' always includes an intercept. Rank-deficient selections return \code{Inf}
#' (which the GA maximizer treats as a very poor score). The value returned is
#' \strong{-BIC} so that GA engines can \emph{maximize} it.
#'
#' @param subset_bin Integer/numeric 0–1 vector (length \code{ncol(X)}); 1 means
#'   the corresponding column of \code{X} is included in the model.
#' @param y Numeric response vector of length \code{n}.
#' @param X Numeric matrix of candidate predictors; columns correspond to variables.
#' @param family A GLM family object (default \code{stats::gaussian()}).
#' @param weights Optional prior weights (passed to \code{glm.fit}).
#' @param offset Optional offset (passed to \code{glm.fit}).
#' @param control GLM fit controls; default \code{stats::glm.control()}.
#'
#' @return A single numeric value: \strong{-BIC}. Larger is better for GA maximizers.
#'   Returns \code{Inf} for rank-deficient designs.
#'
#' @seealso \code{\link[stats]{glm.fit}}, \code{\link[stats]{glm.control}},
#'   \code{\link[stats]{.lm.fit}}
#'
#' @export
#' @importFrom stats glm.fit glm.control .lm.fit
subsetBIC <- function(subset_bin,
                      y,
                      X,
                      family = stats::gaussian(),
                      weights = NULL,
                      offset  = NULL,
                      control = stats::glm.control()) {

  n   <- as.integer(length(y))
  X   <- as.matrix(X)

  idX <- which(as.integer(abs(subset_bin) !=0 ) == 1L)
  m   <- length(idX)

  ## Construct the Design Matrix
  if (m == 0L) {
    Xmat <- matrix(1, nrow = n, ncol = 1L,
                   dimnames = list(NULL, "(Intercept)"))
  } else {
    Xsel <- X[, idX, drop = FALSE]
    Xmat <- cbind(`(Intercept)` = 1, Xsel)
    if (qr(Xmat)$rank < ncol(Xmat)) BIC_val <- Inf
  }

  k_eff <- NCOL(Xmat) # include intercept

  if (identical(family$family, "gaussian") && identical(family$link, "identity")) {
    ## Gaussian-identity: use RSS
    fit <- stats::.lm.fit(Xmat, y)
    if (!is.null(fit$rank) && fit$rank < k_eff) BIC_val <- Inf
    rss_like <- sum(fit$residuals^2)
  }else{
    ## GLM: use deviance
    fit <- stats::glm.fit(x = Xmat, y = y,
                          family  = family,
                          weights = weights,
                          offset  = offset,
                          control = control)
    if (!is.null(fit$rank) && fit$rank < k_eff) BIC_val <- Inf
    rss_like <- fit$deviance
    if (!is.finite(rss_like) || rss_like <= 0) {
      # fall back to a small positive value to avoid -Inf BIC
      rss_like <- .Machine$double.eps
    }
  }

  BIC_val <- n * log(rss_like / n) + k_eff * log(n)

  return(-BIC_val) # GA::ga or GA::gaisl maximizing

}
