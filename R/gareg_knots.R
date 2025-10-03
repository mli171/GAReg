#' Genetic-Algorithm based Optimal Knots Selection
#'
#' @description
#' Runs a GA-based search for changepoints/knots and returns a compact
#' \code{"gareg"} S4 result that stores the backend GA fit
#' (\code{"cptga"} or \code{"cptgaisl"}) plus the essential run settings.
#'
#' @param y Numeric vector of responses (length \code{N}).
#' @param x Optional index or time vector aligned with \code{y}. Currently
#'   kept for API symmetry; not directly used by the backend.
#' @param ObjFunc Objective function or its name. If \code{NULL}, a default
#'   is chosen:
#'   \itemize{
#'     \item \code{fixknotsBIC} when \code{fixedknots} is supplied;
#'     \item \code{varyknotsBIC} otherwise.
#'   }
#'   A custom function must accept the chromosome and needed data via
#'   named arguments (see the defaults for a template).
#' @param fixedknots \code{NULL} (varying-knots search) or an integer vector
#'   of fixed changepoint locations. If non-\code{NULL}, the method is
#'   \code{"fixknots"} and specialized operators are injected unless
#'   overridden in \code{cptgactrl}.
#' @param minDist Integer minimum distance between adjacent changepoints.
#'   If omitted (\code{missing()} or \code{NULL}), the value in \code{cptgactrl}
#'   is used. If supplied here, it overrides the control value.
#' @param polydegree Spline degree for the default BIC objectives.
#' @param gaMethod GA backend to call: function or name. Supports
#'   \code{"cptga"} (single population) and \code{"cptgaisl"} (islands).
#' @param cptgactrl Control list built with \code{\link{cptgaControl}()}
#'   (or a named list of overrides). When \code{gaMethod = "cptgaisl"},
#'   island-specific knobs like \code{numIslands} and \code{maxMig} are
#'   recognized.
#' @param monitoring Logical; print short progress messages (also forwarded
#'   into the backend control).
#' @param seed Optional RNG seed; also stored into the backend control.
#' @param ... Additional arguments passed to the GA backend. If the backend
#'   does not accept \code{...}, unknown arguments are silently dropped
#'   (the call is filtered against the backend formals).
#'
#' @details
#' \strong{Engine selection and controls.}
#' The function detects the engine from \code{gaMethod} and constructs a
#' matching control via \code{\link{cptgaControl}()}:
#' \itemize{
#'   \item \code{"cptga"} uses \code{.cptga.default}.
#'   \item \code{"cptgaisl"} uses \code{.cptgaisl.default} (supports
#'         \code{numIslands}, \code{maxMig}, etc.).
#'   \item see other details in \link[changepointGA:cptga]{cptga} and
#'         \link[changepointGA:cptgaisl]{cptgaisl}.
#' }
#' Top-level \code{monitoring}, \code{seed}, and \code{minDist} given to
#' \code{gareg_knots()} take precedence over the control list.
#'
#' \strong{Fix-knots operators.}
#' When \code{fixedknots} is provided and the control does not already
#' override them, the following operators are injected:
#' \code{Popinitial_fixknots}, \code{crossover_fixknots},
#' \code{mutation_fixknots}.
#'
#' @return An object of class \code{"gareg"} with (key slots):
#' \itemize{
#'   \item \code{call}, \code{method} (\code{"varyknots"} or \code{"fixknots"}),
#'         \code{regMethod}, \code{N}.
#'   \item \code{objFunc}, \code{gaMethod}, \code{gaFit} (the GA result of
#'         class \code{"cptga"} or \code{"cptgaisl"}), \code{ctrl}.
#'   \item \code{fixedknots}, \code{minDist}, \code{polydegree}.
#' }
#' Use \code{summary(g)} to print GA settings and the best solution (extracted
#' from \code{g@gaFit}); \code{show(g)} prints a compact header.
#'
#' @section Argument precedence:
#' Values are combined as \emph{control < core < \code{...}}. That is,
#' \code{cptgactrl} provides defaults, then core arguments from
#' \code{gareg_knots()} override those, and finally any matching names in
#' \code{...} override both.
#'
#' @seealso \code{\link{cptgaControl}}, \code{changepointGA::cptga},
#'   \code{changepointGA::cptgaisl}, \code{\link{fixknotsIC}},
#'   \code{\link{varyknotsIC}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' N <- 120
#' y <- c(rnorm(40, 0), rnorm(40, 3), rnorm(40, 0))
#' x <- seq_len(N)
#'
#' # 1) Varying-knots with single-pop GA
#' g1 <- gareg_knots(
#'   y, x,
#'   minDist = 5,
#'   gaMethod = "cptga",
#'   cptgactrl = cptgaControl(popSize = 150, pcrossover = 0.9, maxgen = 500)
#' )
#' summary(g1)
#'
#' # 2) Fixed knots (operators auto-injected unless overridden)
#' g2 <- gareg_knots(
#'   y, x,
#'   fixedknots = 5,
#'   minDist = 5
#' )
#' summary(g2)
#'
#' # 3) Island GA with island-specific controls
#' g3 <- gareg_knots(
#'   y, x,
#'   gaMethod = "cptgaisl",
#'   minDist = 6,
#'   cptgactrl = cptgaControl(engine = "cptgaisl",
#'                            numIslands = 8, maxMig = 250,
#'                            popSize = 120, pcrossover = 0.9)
#' )
#' summary(g3)
#' }
#'
#' @export
gareg_knots = function(y,
                       x,
                       ObjFunc=NULL,
                       fixedknots=NULL,
                       minDist=3L,
                       polydegree=3L,
                       gaMethod="cptga",
                       cptgactrl=NULL,
                       monitoring=FALSE,
                       seed=NULL,
                       ...) {

  call <- match.call()
  dots <- list(...)

  x_unique <- sort(unique(x))
  n <- length(x_unique)

  gareg_method <- if (is.null(fixedknots)) "varyknots" else "fixknots"
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

  if (missing(minDist) || is.null(minDist)) {
    if (!is.null(cptgactrl$minDist)) minDist <- cptgactrl$minDist
  } else {
    cptgactrl$minDist <- NULL
  }
  if (!is.null(cptgactrl$seed)) set.seed(cptgactrl$seed)

  if (is.null(ObjFunc)) {
    if (!is.null(fixedknots)) {
      if (monitoring) cat("\nDefault Objective Function (fixknotsBIC) in use ...")
      ObjFunc <- fixknotsIC
    } else {
      if (monitoring) cat("\nDefault Objective Function (varyknotsBIC) in use ...")
      ObjFunc <- varyknotsIC
    }
  } else {
    if (monitoring) cat("\nSelf-defined Objective Function in use ...")
    if (is.character(ObjFunc)) ObjFunc <- get(ObjFunc, mode = "function")
    if (!is.function(ObjFunc)) stop("ObjFunc must be a function or a function name.")
  }

  if (!is.null(fixedknots)) {
    if (is.null(cptgactrl$popInitialize) || identical(cptgactrl$popInitialize, .cptga.default$popInitialize))
      cptgactrl$popInitialize <- "Popinitial_fixknots"
    if (is.null(cptgactrl$crossover) || identical(cptgactrl$crossover, .cptga.default$crossover))
      cptgactrl$crossover <- "crossover_fixknots"
    if (is.null(cptgactrl$mutation) || identical(cptgactrl$mutation, .cptga.default$mutation))
      cptgactrl$mutation <- "mutation_fixknots"
  }

  core_args <- list(
    ObjFunc     = ObjFunc,
    N           = n,
    minDist     = minDist,
    y           = y,
    x           = x,
    x_unique    = x_unique,
    polydegree  = polydegree
  )
  if (!is.null(fixedknots)) core_args$fixedknots <- fixedknots

  ga_args <- utils::modifyList(cptgactrl, core_args, keep.null = TRUE)
  if (length(dots)) ga_args <- utils::modifyList(ga_args, dots, keep.null = TRUE)

  ga_formals  <- try(names(formals(ga_fun)),  silent = TRUE);
  if (inherits(ga_formals,  "try-error") || is.null(ga_formals))  ga_formals  <- character(0)
  obj_formals <- try(names(formals(ObjFunc)), silent = TRUE);
  if (inherits(obj_formals, "try-error") || is.null(obj_formals)) obj_formals <- character(0)

  if ("..." %in% ga_formals) {
    keep_names <- union(ga_formals, "ObjFunc")
    keep_names <- union(keep_names, intersect(names(ga_args), obj_formals))
    ga_args <- ga_args[intersect(names(ga_args), keep_names)]
  } else {
    ga_args <- ga_args[intersect(names(ga_args), ga_formals)]
  }

  GA.res <- do.call(ga_fun, ga_args)

  object <- new("gareg",
                call        = call,
                method      = gareg_method,
                N           = n,
                objFunc     = ObjFunc,
                gaMethod    = ga_name,
                gaFit       = GA.res,
                fixedknots  = fixedknots,
                minDist     = minDist,
                polydegree  = polydegree,
                bestFitness = GA.res@overbestfit,
                bestChrom   = GA.res@overbestchrom)

  mhat <- object@bestnumbsol <- object@bestChrom[1]
  if(mhat == 0){
    object@bestsol <- NULL
  }else{
    object@bestsol <- object@bestChrom[2:(1+mhat)]
  }

  return(object)
}


#' Information criterion for B-spline regression with a **variable** number of knots
#'
#' Evaluates an information criterion (BIC, AIC, or AICc) for a regression of
#' \code{y} on a B-spline basis of \code{x} where the **number and locations of
#' interior knots are encoded in the chromosome**. Designed for use as a GA
#' objective/fitness function.
#'
#' @param knot_bin Integer vector (chromosome). The first element is
#'   \eqn{m}, the number of interior knots (\eqn{m \ge 0}). The next
#'   \eqn{m} elements are **indices into** \code{x_unique} selecting interior
#'   knot locations: \code{knot_bin[2L:(1L + m)]}. Indices must be interior,
#'   i.e., in \code{2:(length(x_unique) - 1)}. Any non-finite, duplicate, or
#'   out-of-range index causes the function to return \code{Inf}.
#' @param plen Unused placeholder kept for API compatibility; ignored.
#' @param y Numeric response vector of length \eqn{n}.
#' @param x Numeric predictor (same length as \code{y}) on which the spline
#'   basis is constructed.
#' @param x_unique Optional numeric vector of unique candidate knot locations.
#'   If missing or \code{NULL}, defaults to \code{sort(unique(x))}. Must have
#'   at least three values (two boundaries + one interior) to allow any knots.
#' @param x_base Optional matrix (or vector) of additional covariates to include
#'   linearly alongside the spline basis; coerced to a matrix if supplied.
#' @param polydegree Integer polynomial degree of the B-spline segments
#'   (default \code{3L}).
#' @param ic_method Which information criterion to return: \code{"BIC"},
#'   \code{"AIC"}, or \code{"AICc"}.
#'
#' @details
#' The function:
#' \enumerate{
#' \item Reads \eqn{m <- knot_bin[1]} and the \eqn{m} interior indices
#'   \code{idx <- knot_bin[2L:(1L + m)]}.
#' \item Validates that \code{idx} are finite, interior, and that the resulting
#'   knot locations \code{x_unique[idx]} lie strictly inside \code{range(x)}.
#'   Violations return \code{Inf} (useful for GA pruning).
#' \item Builds the design matrix:
#'   \itemize{
#'     \item If \eqn{m = 0}: \code{X <- cbind(1, x_base)} (no spline terms).
#'     \item Else: \code{X <- cbind(1, x_base, splines::bs(
#'             x, degree = polydegree, knots = x_unique[idx],
#'             Boundary.knots = range(x)))}.
#'   }
#' \item Fits OLS via \code{stats::.lm.fit} and computes the selected criterion
#'   with parameter count \code{p <- NCOL(X)}.
#' }
#'
#' The criteria are computed as:
#' \deqn{\mathrm{BIC} = n \log(\mathrm{SSRes}/n) + p \log n,}
#' \deqn{\mathrm{AIC} = n \log(\mathrm{SSRes}/n) + 2p,}
#' \deqn{\mathrm{AICc} = n \log(\mathrm{SSRes}/n) + 2p +
#'       \frac{2p(p+1)}{n-p-1},}
#' where \eqn{\mathrm{SSRes}} is the residual sum of squares and \eqn{p} is the
#' number of columns in \code{X}.
#'
#' @return A single numeric value: the requested information criterion (lower is
#'   better). Returns \code{Inf} for invalid chromosomes/inputs.
#'
#' @note This function allows \eqn{m=0} (no spline terms) so that the GA can
#'   compare against a pure-linear baseline (intercept + \code{x_base}).
#'   Spacing constraints (e.g., minimum distance between indices) should be
#'   enforced by the GA operators or an external penalty.
#'
#' @seealso [fixknotsIC()], [splines::bs()]
#'
#' @examples
#' ## Example with 'mcycle' data (MASS)
#' # y <- mcycle$accel; x <- mcycle$times
#' # x_unique <- sort(unique(x))
#' # # chromosome encoding m=4 and four interior indices:
#' # chrom <- c(4, 25, 40, 52, 77)
#' # varyknotsIC(chrom, y = y, x = x, x_unique = x_unique, ic_method = "BIC")
#' @export
varyknotsIC <- function(knot_bin,
                        plen = 0,
                        y,
                        x,
                        x_unique,
                        x_base = NULL,
                        polydegree = 3L,
                        ic_method = "BIC") {

  ic_method <- match.arg(ic_method)
  stopifnot(is.numeric(y), is.numeric(x), length(y) == length(x))
  if (!is.null(x_base)) x_base <- as.matrix(x_base)

  if (missing(x_unique) || is.null(x_unique)) {
    x_unique <- sort(unique(x))
  } else {
    x_unique <- sort(unique(x_unique))
  }
  if (length(x_unique) < 3L) {
    return(Inf)
  } # need interior points

  n <- length(y)
  m <- as.integer(knot_bin[1L])

  xr <- range(x, na.rm = TRUE)

  # Design matrix
  ones <- rep(1, n)
  if (m == 0L) {
    X <- cbind(ones, x_base)
  } else {
    idx <- as.integer(knot_bin[2L:(1L + m)])
    knotvec <- sort(x_unique[idx])

    L <- 2L
    U <- length(x_unique) - 1L
    if (U < L) {
      return(Inf)
    }
    if (any(!is.finite(idx)) || any(idx < L) || any(idx > U)) {
      return(Inf)
    }

    if (length(knotvec) != m) {
      return(Inf)
    }
    if (anyNA(knotvec)) {
      return(Inf)
    }
    if (any(!(knotvec > xr[1] & knotvec < xr[2]))) {
      return(Inf)
    }

    x_bs <- splines::bs(x,
                        degree = polydegree,
                        knots = knotvec,
                        Boundary.knots = xr
    )
    X <- cbind(ones, x_base, x_bs)
  }

  fit <- stats::.lm.fit(X, y)
  SSRes <- sum(fit$residuals^2)
  p <- NCOL(X)

  val <- switch(ic_method,
                BIC  = n * log(SSRes / n) + p * log(n),
                AIC  = n * log(SSRes / n) + p * 2,
                AICc = if (n - p - 1 > 0) n * log(SSRes / n) + p * 2 + (2 * p * (p + 1)) / (n - p - 1) else Inf
  )

  return(val)
}

#' Information criterion for a fixed–knot B-spline regression
#'
#' Computes an information criterion (BIC, AIC, or AICc) for a regression of
#' \code{y} on a B-spline basis of \code{x} when the number of interior knots is
#' fixed. This is designed to be used as a fitness/objective function inside a
#' GA search where the chromosome encodes the **indices** of the interior knots.
#'
#' @param knot_bin Integer vector (chromosome). The first element may be ignored
#'   in a fixed-\eqn{m} setup; the next \eqn{m} elements are **indices into**
#'   \code{x_unique} selecting the interior knots:
#'   \code{knot_bin[2L:(1L + m)]}. Indices must refer to interior candidates
#'   \code{2:(length(x_unique)-1)}. Non-finite, duplicate, out-of-range, or
#'   boundary indices cause the function to return \code{Inf}.
#' @param plen Unused placeholder kept for API compatibility with other
#'   objective functions. Ignored.
#' @param y Numeric response vector of length \eqn{n}.
#' @param x Numeric predictor (same length as \code{y}) on which the spline
#'   basis is built.
#' @param x_unique Optional numeric vector of unique candidate knot locations.
#'   If \code{NULL} or missing, it defaults to \code{sort(unique(x))}. Must
#'   contain at least \eqn{m + 2} values (interior + two boundaries).
#' @param x_base Optional matrix (or vector) of additional covariates to include
#'   linearly alongside the spline basis. If supplied, it is coerced to a
#'   matrix and column-bound to the design.
#' @param fixedknots Integer \eqn{m}: the number of **interior** knots to use.
#'   Internally this determines how many indices are read from \code{knot_bin}.
#' @param polydegree Integer polynomial degree of the B-spline segments
#'   (default \code{3L}).
#' @param ic_method Character; which information criterion to return:
#'   \code{"BIC"}, \code{"AIC"}, or \code{"AICc"}.
#'
#' @details
#' The function:
#' \enumerate{
#' \item Derives \eqn{m <- as.integer(fixedknots)} and reads indices
#'   \code{idx <- knot_bin[2L:(1L + m)]}.
#' \item Validates that \code{idx} are finite, interior, unique, and sorted; the
#'   corresponding knot locations \code{x_unique[idx]} must lie strictly inside
#'   \code{range(x)}. Any violation returns \code{Inf} (useful for GA pruning).
#' \item Builds the design matrix
#'   \code{cbind(1, x_base, splines::bs(x, degree = polydegree,
#'   knots = x_unique[idx], Boundary.knots = range(x)))}.
#' \item Fits OLS via \code{stats::.lm.fit} and computes the selected criterion
#'   with parameter count \code{p = NCOL(X)} (intercept + \code{x_base} columns
#'   + spline basis columns).
#' }
#'
#' @return A single numeric value: the requested information criterion. Lower is
#'   better. Returns \code{Inf} for invalid chromosomes/inputs.
#'
#' @note Earlier drafts mistakenly referenced an undefined/global \code{m}.
#' Always set \code{m <- as.integer(fixedknots)} inside the function to ensure
#' that exactly \code{fixedknots} interior knots are evaluated.
#'
#' @seealso [varyknotsIC()], [splines::bs()]
#'
#' @examples
#' library(MASS)
#' y <- mcycle$accel; x <- mcycle$times
#' x_unique <- sort(unique(x))
#' # chromosome encoding 5 interior knot indices:
#' chrom <- c(5, 24, 30, 46, 49, 69, 95)
#' fixknotsIC(chrom, y = y, x = x, x_unique = x_unique,
#'            fixedknots = 5, ic_method = "BIC")
#' @export
fixknotsIC <- function(knot_bin,
                       plen=0,
                       y,
                       x,
                       x_unique,
                       x_base=NULL,
                       fixedknots,
                       polydegree=3L,
                       ic_method = "BIC"){

  ic_method <- match.arg(ic_method)
  stopifnot(is.numeric(y), is.numeric(x), length(y) == length(x))
  if (!is.null(x_base)) x_base <- as.matrix(x_base)

  if (missing(x_unique) || is.null(x_unique)) {
    x_unique <- sort(unique(x))
  } else {
    x_unique <- sort(unique(x_unique))
  }
  if (length(x_unique) < fixedknots) {
    stop("unique value less than fixed knots")
  } # need interior points

  n <- length(y)
  xr <- range(x, na.rm = TRUE)

  # Design matrix
  idx <- as.integer(knot_bin[2L:(1L + fixedknots)])
  knotvec <- sort(x_unique[idx])

  L <- 2L
  U <- length(x_unique) - 1L
  if (U < L) {
    return(Inf)
  }
  if (any(!is.finite(idx)) || any(idx < L) || any(idx > U)) {
    return(Inf)
  }
  if (anyNA(knotvec)) {
    return(Inf)
  }
  if (any(!(knotvec > xr[1] & knotvec < xr[2]))) {
    return(Inf)
  }

  x_bs <- splines::bs(x,
                      degree = polydegree,
                      knots = knotvec,
                      Boundary.knots = xr
  )
  X <- cbind(rep(1, n), x_base, x_bs)

  fit <- stats::.lm.fit(X, y)
  SSRes <- sum(fit$residuals^2)
  p <- NCOL(X)

  val <- switch(ic_method,
                BIC  = n * log(SSRes / n) + p * log(n),
                AIC  = n * log(SSRes / n) + p * 2,
                AICc = if (n - p - 1 > 0) n * log(SSRes / n) + p * 2 + (2 * p * (p + 1)) / (n - p - 1) else Inf
  )

  return(val)
}


.is_scalar <- function(x) is.atomic(x) && length(x) == 1L
.is_intish <- function(x) is.numeric(x) && length(x)==1L && is.finite(x) && abs(x - round(x)) < 1e-9
.chk_prob  <- function(x, nm) if(!(is.numeric(x)&&length(x)==1L&&is.finite(x)&&x>=0&&x<=1))
  stop(sprintf("`%s` must be a single number in [0,1].", nm), call.=FALSE)

.validate_ctrl <- function(x, engine = c("cptga","cptgaisl")) {
  engine <- match.arg(engine)
  if(!.is_intish(x$popSize) || x$popSize < 1) stop("`popSize` must be integer >= 1.")
  .chk_prob(x$pcrossover, "pcrossover")
  .chk_prob(x$pmutation,  "pmutation")
  .chk_prob(x$pchangepoint, "pchangepoint")
  if(!.is_intish(x$minDist) || x$minDist < 0) stop("`minDist` must be integer >= 0.")
  if(!is.null(x$mmax)  && (!.is_intish(x$mmax)  || x$mmax  < 0)) stop("`mmax` must be NULL or integer >= 0.")
  if(!is.null(x$lmax)  && (!.is_intish(x$lmax)  || x$lmax  < 0)) stop("`lmax` must be NULL or integer >= 0.")
  if(!.is_intish(x$maxgen)  || x$maxgen  < 1)  stop("`maxgen` must be integer >= 1.")
  if(!.is_intish(x$maxconv) || x$maxconv < 1)  stop("`maxconv` must be integer >= 1.")
  if(!(x$option %in% c("cp","both"))) stop("`option` must be 'cp' or 'both'.")
  if(!is.logical(x$monitoring) || length(x$monitoring)!=1L) stop("`monitoring` must be logical(1).")
  if(!is.logical(x$parallel)   || length(x$parallel)!=1L)   stop("`parallel` must be logical(1).")
  if(!is.null(x$nCore) && (!.is_intish(x$nCore) || x$nCore < 1)) stop("`nCore` must be NULL or integer >= 1.")
  if(!(is.numeric(x$tol) && length(x$tol)==1L && is.finite(x$tol) && x$tol > 0)) stop("`tol` must be > 0.")
  if(!is.null(x$seed) && !.is_intish(x$seed)) stop("`seed` must be NULL or an integer.")
  if(!(is.character(x$popInitialize) || is.function(x$popInitialize))) stop("`popInitialize` must be name or function.")
  if(!(is.character(x$selection)     || is.function(x$selection)))     stop("`selection` must be name or function.")
  if(!(is.character(x$crossover)     || is.function(x$crossover)))     stop("`crossover` must be name or function.")
  if(!(is.character(x$mutation)      || is.function(x$mutation)))      stop("`mutation` must be name or function.")
  if (engine == "cptgaisl") {
    if(!.is_intish(x$numIslands) || x$numIslands < 1) stop("`numIslands` must be integer >= 1.")
    if(!.is_intish(x$maxMig)     || x$maxMig < 0)     stop("`maxMig` must be integer >= 0.")
  }
  x
}

#' Build Control List for \code{cptga}/\code{cptgaisl}
#'
#' @description
#' Convenience constructor for GA control parameters used by
#' \code{changepointGA::cptga} and \code{changepointGA::cptgaisl}. It merges
#' named overrides into engine-specific defaults
#' (\link{.cptga.default} or \link{.cptgaisl.default}), with light validation.
#'
#' @param ... Named overrides for control fields (e.g., \code{popSize},
#'   \code{pcrossover}, \code{minDist}, \code{numIslands}).
#' @param .list Optional named list of overrides (merged with \code{...}).
#' @param .persist Logical; if \code{TRUE}, persist updated defaults back into
#'   the target environment (not usually recommended in user code).
#' @param .env Environment where defaults live (defaults to \code{parent.frame()}).
#' @param .validate Logical; validate values/ranges (default \code{TRUE}).
#' @param engine Character; one of \code{"cptga"} or \code{"cptgaisl"} to select
#'   the default set and validation rules.
#'
#' @details
#' Unknown names are rejected. When both \code{...} and \code{.list} are present,
#' they are combined, with later entries overwriting earlier ones.
#'
#' @return A list of class \code{"cptgaControl"}.
#'
#' @seealso \link{gareg_knots}, \link{.cptga.default}, \link{.cptgaisl.default}
#' @export
cptgaControl <- function(..., .list = NULL, .persist = FALSE,
                         .env = asNamespace("GAReg"), .validate = TRUE,
                         engine = NULL) {

  overrides <- list(...)
  if (!is.null(.list)) {
    if (!is.list(.list)) stop("`.list` must be a named list.")
    overrides <- c(overrides, list(.list))
  }
  if (length(overrides) == 1L && is.list(overrides[[1]]) && !length(names(overrides)))
    overrides <- overrides[[1]]

  if (!is.null(engine)) {
    engine <- match.arg(engine, c("cptga","cptgaisl"))
  } else {
    nm <- names(overrides)
    engine <- if (length(nm) && any(c("numIslands","maxMig") %in% nm)) "cptgaisl" else "cptga"
  }

  defaults_name <- if (engine == "cptga") ".cptga.default" else ".cptgaisl.default"
  if (!exists(defaults_name, envir = .env, inherits = TRUE))
    stop(sprintf("`%s` not found in target environment.", defaults_name))
  current <- get(defaults_name, envir = .env, inherits = TRUE)

  if (!length(overrides)) return(structure(current, class = c("cptgaControl","list")))

  if (is.null(names(overrides)) || any(names(overrides) == ""))
    stop("All control overrides must be *named*.")

  unknown <- setdiff(names(overrides), names(current))
  if (length(unknown))
    stop("Unknown control parameter(s): ", paste(unknown, collapse = ", "))

  updated <- current
  updated[names(overrides)] <- overrides
  if (.validate) updated <- .validate_ctrl(updated, engine = engine)

  if (.persist) {
    assign(defaults_name, updated, envir = .env)
    invisible(structure(updated, class = c("cptgaControl","list")))
  } else {
    structure(updated, class = c("cptgaControl","list"))
  }
}

#' Fixed-Knots Population Initializer
#'
#' @description
#' Initializes a population matrix for the fixed-knots GA. Each column is a
#' feasible chromosome sampled by \link{selectTau_uniform_exact}.
#'
#' @param popSize Integer; number of individuals (columns).
#' @param prange Optional hyperparameter range (unused here).
#' @param N Series length.
#' @param minDist Integer minimum spacing between adjacent changepoints.
#' @param Pb Unused placeholder (kept for compatibility).
#' @param mmax,lmax Integers; maximum number of knots and chromosome length.
#' @param fixedknots Integer; number of knots to place.
#'
#' @return Integer matrix of size \code{lmax x popSize}; each column is a
#'   chromosome \code{c(m, tau_1, ..., tau_m, N+1, ...)}.
#'
#' @seealso \link{selectTau_uniform_exact}, \link{gareg_knots}
#' @export
Popinitial_fixknots <- function(popSize, prange=NULL, N, minDist, Pb, mmax, lmax, fixedknots){

  pop <- matrix(0, nrow=lmax, ncol=popSize)

  for(j in 1:popSize){
    pop[,j] = selectTau_uniform_exact(N, fixedknots, minDist, lmax)
  }

  return(pop)
}

#' Crossover Operator (Fixed-\eqn{m}) with Feasibility-First Restarts
#'
#' @description
#' Produces a child chromosome from two fixed-\eqn{m} parents (same number of
#' knots) by alternately sampling candidate knot locations from the parents and
#' enforcing the spacing constraint \code{diff(child) > minDist}. If a conflict
#' is encountered, the routine restarts the construction up to a small cap.
#'
#' @details
#' Let \code{mom} and \code{dad} be chromosomes of the form
#' \code{c(m, tau_1, ..., tau_m, ...)}. This operator:
#' \enumerate{
#'   \item Initializes an empty child of size \eqn{m}.
#'   \item Picks the first knot at random from \code{mom} or \code{dad}.
#'   \item For each subsequent position \eqn{i=2,\dots,m}, considers the
#'         pair \code{(mom[i], dad[i])} and chooses the first value that
#'         maintains the spacing constraint relative to the previously chosen
#'         knot (\code{> minDist}); if both work, one is chosen at random.
#'   \item If no feasible choice exists at some step, the construction restarts
#'         from the first position (up to a small cap governed internally by
#'         \code{up_tol}).
#' }
#' The result is written back as a full-length chromosome with the sentinel
#' \code{N+1} in position \code{m+2}, and zeros elsewhere.
#'
#' @param mom,dad Integer vectors encoding parent chromosomes:
#'   first entry \eqn{m} (number of changepoints), followed by \eqn{m} ordered
#'   knot locations.
#' @param prange Unused placeholder (kept for compatibility with other GA
#'   operators). Default \code{NULL}.
#' @param minDist Integer; minimum spacing between adjacent knots in the child.
#' @param lmax Integer; chromosome length (number of rows in the population
#'   matrix).
#' @param N Integer; series length. Used to place the sentinel \code{N+1} at
#'   position \code{m+2}.
#'
#' @return
#' An integer vector of length \code{lmax} encoding the child chromosome:
#' \code{c(m, child_knots, N+1, 0, 0, ...)}.
#'
#' @seealso
#' \link{crossover_fixknots}, \link{mutation_fixknots},
#' \link{selectTau_uniform_exact}, \link{Popinitial_fixknots},
#' \link{gareg_knots}
#'
#' @examples
#' \dontrun{
#' N <- 120; lmax <- 30; minDist <- 5
#' m <- 3
#' mom <- c(m, c(20, 50, 90), rep(0, lmax - 1 - m)); mom[m+2] <- N + 1
#' dad <- c(m, c(18, 55, 85), rep(0, lmax - 1 - m)); dad[m+2] <- N + 1
#' child <- crossover_fixknots(mom, dad, minDist = minDist, lmax = lmax, N = N)
#' child
#' }
#'
#' @export
crossover_fixknots <- function(mom, dad, prange=NULL, minDist, lmax, N){

  up_tol <- 30L
  max_restarts <- 50L

  output <- rep.int(0L, lmax)

  m_child <- as.integer(dad[1L])
  child <- rep.int(NA_integer_, m_child)

  if (all(mom[2:(m_child+1)] == dad[2:(m_child+1)])) {
    mom <- selectTau_uniform_exact(N, m_child, minDist, lmax)
  }

  mom_tau <- as.integer(mom[2L:(m_child + 1L)])
  dad_tau <- as.integer(dad[2L:(m_child + 1L)])
  co_tab  <- as.integer(as.vector(rbind(mom_tau, dad_tau)))

  i <- 1L
  ii <- 0L
  restarts <- 0L
  tmppick <- NA_integer_
  while (i <= m_child) {
    if (i == 1L) {
      # first pick random from first pair of mom and dad
      tmppick   <- sample(co_tab[1:2], size = 1L)
      child[1L] <- tmppick
      i <- i + 1L
    } else {
      # pick from the i^th pair of mom and dad (many cases)
      tmp_co_tab <- co_tab[((i - 1L) * 2L + 1L):(i * 2L)]
      ok <- which(tmp_co_tab - tmppick > minDist)
      if (length(ok) >= 1L) {
        tmppick <- if (length(ok) == 2L) sample(tmp_co_tab[ok], 1L) else tmp_co_tab[ok]
        child[i] <- tmppick

        if (i == m_child) {
          if (all(child == dad_tau) || all(child == mom_tau)) {
            if (ii > up_tol) break
            i <- 1L; child[] <- NA_integer_; restarts <- restarts + 1L
          } else {
            i <- i + 1L
          }
        } else {
          i <- i + 1L
        }
      } else {
        i <- 1L; child[] <- NA_integer_; restarts <- restarts + 1L
      }
    }

    ii <- ii + 1L
    if (restarts > max_restarts) {
      # select a new one to guaranteed progress
      child <- selectTau_uniform_exact(N, m_child, minDist, lmax)[2L:(m_child + 1L)]
      break
    }
  }

  output[1] <- m_child
  output[2:(m_child+1)] <- child
  output[m_child+2] <- N+1

  return(output)
}

#' Exact Uniform Sampler of Feasible Changepoints
#'
#' @description
#' Samples \eqn{m} ordered changepoint indices uniformly from all feasible
#' configurations on \code{1:N} subject to a minimum spacing \code{minDist}.
#' Encodes the result as a chromosome for downstream GA operators.
#'
#' @param N Integer series length.
#' @param m Integer number of changepoints to place.
#' @param minDist Integer minimum spacing between adjacent changepoints.
#' @param lmax Integer chromosome length.
#'
#' @return Integer vector length \code{lmax}:
#'   \code{c(m, tau_1, ..., tau_m, N+1, 0, 0, ...)}.
#'
#' @seealso \link{Popinitial_fixknots}, \link{mutation_fixknots}
#' @export
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

#' Mutation Operator (Fixed-Knots)
#'
#' @description
#' Replaces a child with a fresh feasible sample having the same \eqn{m},
#' drawn by \link{selectTau_uniform_exact}.
#'
#' @param child Current chromosome (its first entry defines \eqn{m}).
#' @param p.range,Pb Unused placeholders (kept for compatibility).
#' @param minDist Integer minimum spacing.
#' @param lmax,mmax Integers; chromosome length and maximum \eqn{m} (unused).
#' @param N Integer series length.
#'
#' @return New feasible chromosome with the same \eqn{m}.
#'
#' @seealso \link{crossover_fixknots}
#' @export
mutation_fixknots <- function(child, p.range = NULL, minDist, Pb, lmax, mmax, N) {

  m <- child[1]
  childMut <- selectTau_uniform_exact(N, m, minDist, lmax)

  return(childMut)
}

#' Default Controls for \code{cptga}
#'
#' @description
#' Engine defaults used by \link{cptgaControl} when \code{engine = "cptga"}.
#' Not exported; shown here for reference.
#'
#' @format A named list with fields like \code{popSize}, \code{pcrossover},
#' \code{pmutation}, \code{pchangepoint}, \code{minDist}, \code{maxgen},
#' \code{option}, \code{selection}, \code{crossover}, \code{mutation}, etc.
#'
#' @seealso \link{cptgaControl}, \link{.cptgaisl.default}
#' @keywords internal
.cptga.default <- list(
  popSize = 200,
  pcrossover = 0.95,
  pmutation = 0.3,
  pchangepoint = 0.01,
  minDist = 1,
  mmax = NULL,
  lmax = NULL,
  maxgen = 50000,
  maxconv = 5000,
  option = "cp",
  monitoring = FALSE,
  parallel = FALSE,
  nCore = NULL,
  tol = 1e-05,
  seed = NULL,
  popInitialize = "random_population",
  suggestions = NULL,
  selection = "selection_linearrank",
  crossover = "uniformcrossover",
  mutation = "mutation")

#' Default Controls for \code{cptgaisl} (Island GA)
#'
#' @description
#' Engine defaults used by \link{cptgaControl} when \code{engine = "cptgaisl"}.
#' Includes island-specific fields (e.g., \code{numIslands}, \code{maxMig}).
#'
#' @format A named list with fields like \code{popSize}, \code{numIslands},
#' \code{pcrossover}, \code{pmutation}, \code{maxMig}, \code{maxgen}, etc.
#'
#' @seealso \link{cptgaControl}, \link{.cptga.default}
#' @keywords internal
.cptgaisl.default <- list(
  popSize = 200,
  numIslands = 5,
  pcrossover = 0.95,
  pmutation  = 0.3,
  pchangepoint = 0.01,
  minDist = 1,
  mmax = NULL,
  lmax = NULL,
  maxMig = 1000,
  maxgen = 50,
  maxconv = 100,
  option = "cp",
  monitoring = FALSE,
  parallel = FALSE,
  nCore = NULL,
  tol = 1e-05,
  seed = NULL,
  popInitialize = "random_population",
  suggestions = NULL,
  selection = "selection_linearrank",
  crossover = "uniformcrossover",
  mutation = "mutation"
)
