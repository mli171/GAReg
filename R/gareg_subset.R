#' Genetic-Algorithm based Best Subset Selection
#'
#' @description
#' Runs a GA search over variable subsets using a BIC-like objective
#' (default: \code{subsetBIC}) and returns a compact \code{"gareg"} S4 result.
#'
#' @param y Numeric response vector (length \code{N}).
#' @param X Numeric matrix of candidate predictors.
#' @param ObjFunc Objective function or its name. Defaults to \code{subsetBIC}.
#'   A custom function must accept the chromosome and needed data via named args
#'   (at least \code{y}, \code{X}).
#' @param gaMethod GA backend to call: function or name. Supports
#'   \code{"cptga"} (single population) and \code{"cptgaisl"} (islands).
#' @param cptgactrl Control list from \code{cptgaControl()} (or a named list of
#'   overrides). If \code{gaMethod = "cptgaisl"}, island knobs like
#'   \code{numIslands}, \code{maxMig} are recognized.
#' @param monitoring Logical; print short progress messages (also forwarded
#'   into the backend control).
#' @param seed Optional RNG seed; also stored into the backend control.
#' @param ... Additional arguments passed to the GA backend and then to the
#'   objective (e.g., \code{family}, \code{weights}, \code{offset}).
#'
#' @return An object of class \code{"gareg"} with \code{method="subset"}.
#'   Use \code{summary()} to see GA settings and the best subset.
#' @export
gareg_subset = function(y,
                        X,
                        ObjFunc=NULL,
                        minDist=1,
                        gaMethod="cptga",
                        cptgactrl=NULL,
                        monitoring=FALSE,
                        seed=NULL,
                        ...){

  call <- match.call()
  dots <- list(...)
  n <- NCOL(X)

  ga_name <- if (is.function(gaMethod)) deparse(substitute(gaMethod)) else as.character(gaMethod)
  ga_fun  <- if (is.function(gaMethod)) gaMethod else get(gaMethod, mode = "function")
  engine  <- if (tolower(ga_name) == "cptgaisl") "cptgaisl" else "cptga"

  cptgactrl <- if (is.null(cptgactrl)) {
    cptgaControl(engine = engine)
  } else if (inherits(cptgactrl, "cptgaControl")) {
    cptgaControl(.list = unclass(cptgactrl), engine = engine)
  } else if (is.list(cptgactrl)) {
    cptgaControl(.list = cptgactrl, engine = engine)
  } else {
    stop("'cptgactrl' must be cptgaControl() or a named list.")
  }

  if (!missing(monitoring) && !is.null(monitoring)) cptgactrl$monitoring <- monitoring
  if (!missing(seed)       && !is.null(seed))       cptgactrl$seed       <- seed
  if (!is.null(cptgactrl$seed)) set.seed(cptgactrl$seed)

  if (is.null(ObjFunc)) {
    if (monitoring) cat("\nDefault Objective Function (subsetBIC) in use ...")
    ObjFunc <- subsetBIC
  } else {
    if (monitoring) cat("\nSelf-defined Objective Function in use ...")
    if (is.character(ObjFunc)) ObjFunc <- get(ObjFunc, mode = "function")
    if (!is.function(ObjFunc)) stop("ObjFunc must be a function or a function name.")
  }

  core_args <- list(
    ObjFunc     = ObjFunc,
    N           = n,
    y           = y,
    X           = X
  )

  # Combine: control < core < dots
  ga_args <- utils::modifyList(cptgactrl, core_args, keep.null = TRUE)
  if (length(dots)) ga_args <- utils::modifyList(ga_args, dots, keep.null = TRUE)

  # If backend exposes 'option', request subset mode (only if supported)
  fm <- try(names(formals(ga_fun)), silent = TRUE)
  if (!inherits(fm, "try-error")) {
    if ("option" %in% fm) ga_args$option <- "subset"
    if (!("..." %in% fm)) {
      ga_args <- ga_args[intersect(names(ga_args), fm)]
    }
  }

  GA.res <- do.call(ga_fun, ga_args)

  object <- new("gareg",
                call        = call,
                method      = "subset",
                N           = n,
                objFunc     = ObjFunc,
                gaMethod    = ga_name,
                gaFit       = GA.res,
                minDist     = minDist,
                bestFitness = GA.res@overbestfit,
                bestChrom   = GA.res@overbestchrom)

  mhat <- object@bestnumbsol <- as.integer(object@bestChrom[1])
  if (is.na(mhat) || mhat <= 0) {
    object@bestsol <- numeric()
  } else {
    object@bestsol <- as.integer(object@bestChrom[2:(1 + mhat)])
  }

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
#' The chromosome \code{subset_bin} is assumed to store:
#' \itemize{
#'   \item \code{subset_bin[1]} = \eqn{m}, the number of selected predictors,
#'   \item \code{subset_bin[2:(1+m)]} = column indices in \code{X} to include.
#' }
#' The design matrix always includes an intercept. Rank-deficient selections are
#' rejected by returning \code{Inf} so that the GA avoids them. For non-Gaussian
#' families, if the fitted deviance is non-finite or non-positive, a very small
#' positive fallback is used to keep the criterion defined.
#'
#' @param subset_bin Integer vector encoding the subset: first element \eqn{m},
#'   followed by \eqn{m} column indices into \code{X}.
#' @param plen Unused placeholder (kept for compatibility with \code{changepointGA} interfaces).
#' @param y Numeric response vector of length \code{n}.
#' @param X Numeric matrix of candidate predictors; columns are referenced by
#'   \code{subset_bin}.
#' @param family A GLM family object (default \code{stats::gaussian()}).
#' @param weights Optional prior weights (passed to \code{glm.fit}).
#' @param offset Optional offset (passed to \code{glm.fit}).
#' @param control GLM fit controls; default \code{stats::glm.control()}.
#'
#' @return A single numeric value: the BIC-like score
#'   \code{n * log(rss_like / n) + k * log(n)}. Smaller is better.
#'   Returns \code{Inf} for rank-deficient designs.
#'
#' @seealso \code{\link[stats]{glm.fit}}, \code{\link[stats]{glm.control}},
#'   \code{\link[stats]{.lm.fit}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 100; p <- 10
#' X <- matrix(rnorm(n*p), n, p)
#' y <- 1 + X[,1] - 0.5*X[,3] + rnorm(n)
#'
#' # choose variables 1 and 3 (m = 2)
#' chrom <- c(2L, 1L, 3L)
#'
#' # Gaussian (fast RSS path)
#' subsetBIC(chrom, y = y, X = X, family = stats::gaussian())
#'
#' # Logistic (use deviance)
#' ybin <- rbinom(n, size = 1, prob = plogis(0.5 + X[,1] - X[,2]))
#' subsetBIC(chrom, y = ybin, X = X, family = stats::binomial())
#' }
#'
#' @export
#' @importFrom stats glm.fit glm.control .lm.fit
subsetBIC <- function(subset_bin,
                      plen=0,
                      y,
                      X,
                      family = stats::gaussian(),
                      weights = NULL,
                      offset  = NULL,
                      control = stats::glm.control()) {

  n   <- as.integer(length(y))
  X   <- as.matrix(X)
  m   <- as.integer(subset_bin[1])

  ## Construct the Design Matrix
  if (m == 0L) {
    Xmat <- matrix(1, nrow = n, ncol = 1L,
                   dimnames = list(NULL, "(Intercept)"))
  } else {
    idX  <- subset_bin[2:(1 + m)]
    Xsel <- X[, idX, drop = FALSE]
    Xmat <- cbind(`(Intercept)` = 1, Xsel)
    if (qr(Xmat)$rank < ncol(Xmat)) BIC_val <- Inf
  }

  k_eff <- NCOL(Xmat) # include intercept

  if (identical(family$family, "gaussian") && identical(family$link, "identity")) {
    ## Gaussian-identity: use RSS
    fit <- stats::.lm.fit(Xmat, y)
    if (!is.null(fit$rank) && fit$rank < p) BIC_val <- Inf
    rss_like <- sum(fit$residuals^2)
  }else{
    ## GLM: use deviance
    fit <- stats::glm.fit(x = Xmat, y = y,
                          family  = family,
                          weights = weights,
                          offset  = offset,
                          control = control)
    if (!is.null(fit$rank) && fit$rank < p) BIC_val <- Inf
    rss_like <- fit$deviance
    if (!is.finite(rss_like) || rss_like <= 0) {
      # fall back to a small positive value to avoid -Inf BIC
      rss_like <- .Machine$double.eps
    }
  }

  BIC_val <- n * log(rss_like / n) + k_eff * log(n)

  return(BIC_val)

}
