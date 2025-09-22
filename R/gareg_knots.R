#' Genetic-Algorithm based Optimal Knots Selection
#'
#' @description
#' Runs a GA-based search for changepoints/knots and returns a compact
#' \code{"gareg"} S4 result that stores the backend GA fit
#' (\code{"cptga"} or \code{"cptgaisl"}) plus the essential run settings.
#'
#' @param y Numeric vector of responses (length \code{N}).
#' @param t Optional index or time vector aligned with \code{y}. Currently
#'   kept for API symmetry; not directly used by the backend.
#' @param X Optional matrix of additional covariates (ignored by the default
#'   objectives, but available to custom objectives).
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
#' }
#' Top-level \code{monitoring}, \code{seed}, and \code{minDist} given to
#' \code{gareg_knots()} take precedence over the control list.
#'
#' \strong{Fix-knots operators.}
#' When \code{fixedknots} is provided and the control does not already
#' override them, the following operators are injected:
#' \code{Popinitial_fixknots}, \code{crossover_fixknots_2},
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
#'   \code{changepointGA::cptgaisl}, \code{\link{fixknotsBIC}},
#'   \code{\link{varyknotsBIC}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' N <- 120
#' y <- c(rnorm(40, 0), rnorm(40, 3), rnorm(40, 0))
#' t <- seq_len(N)
#'
#' # 1) Varying-knots with single-pop GA
#' g1 <- gareg_knots(
#'   y, t,
#'   minDist = 5,
#'   gaMethod = "cptga",
#'   cptgactrl = cptgaControl(popSize = 150, pcrossover = 0.9, maxgen = 500)
#' )
#' summary(g1)
#'
#' # 2) Fixed knots (operators auto-injected unless overridden)
#' g2 <- gareg_knots(
#'   y, t,
#'   fixedknots = 5,
#'   minDist = 5
#' )
#' summary(g2)
#'
#' # 3) Island GA with island-specific controls
#' g3 <- gareg_knots(
#'   y, t,
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
                       t,
                       X=NULL,
                       ObjFunc=NULL,
                       fixedknots=NULL,
                       minDist=3,
                       polydegree=3,
                       gaMethod="cptga",
                       cptgactrl=NULL,
                       monitoring=FALSE,
                       seed=NULL,
                       ...) {

  call <- match.call()
  dots <- list(...)
  n <- length(y)

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
      ObjFunc <- fixknotsBIC
    } else {
      if (monitoring) cat("\nDefault Objective Function (varyknotsBIC) in use ...")
      ObjFunc <- varyknotsBIC
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
      cptgactrl$crossover <- "crossover_fixknots_2"
    if (is.null(cptgactrl$mutation) || identical(cptgactrl$mutation, .cptga.default$mutation))
      cptgactrl$mutation <- "mutation_fixknots"
  }

  core_args <- list(
    ObjFunc     = ObjFunc,
    N           = n,
    minDist     = minDist,
    y           = y,
    polydegree  = polydegree
  )
  if (!is.null(fixedknots)) core_args$fixedknots <- fixedknots

  ga_args <- utils::modifyList(cptgactrl, core_args, keep.null = TRUE)
  if (length(dots)) ga_args <- utils::modifyList(ga_args, dots, keep.null = TRUE)

  fm <- try(names(formals(ga_fun)), silent = TRUE)
  if (!inherits(fm, "try-error") && !("..." %in% fm)) {
    ga_args <- ga_args[intersect(names(ga_args), fm)]
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


#' BIC Objective for Varying-Knots GA
#'
#' @description
#' Computes a BIC-style objective for a chromosome encoding \emph{varying}
#' changepoints/knots. The chromosome stores the number of knots \eqn{m} in
#' \code{knot_bin[1]} and the ordered knot locations in \code{knot_bin[2:(m+1)]}.
#' A spline basis of degree \code{polydegree} is built on \code{1:n}, with an
#' intercept and optional \code{x_base}, and the BIC is returned.
#'
#' @param knot_bin Numeric vector; GA chromosome with \eqn{m} and ordered knots.
#' @param plen Unused placeholder (kept for compatibility).
#' @param y Numeric response vector.
#' @param x_base Optional matrix of additional covariates (will be column-bound).
#' @param polydegree Integer spline degree (default \code{3L}).
#'
#' @return Numeric scalar: BIC value (lower is better).
#'
#' @seealso \link{fixknotsBIC}, \link{gareg_knots}
#' @export
varyknotsBIC <- function(knot_bin,
                         plen=0,
                         y,
                         x_base=NULL,
                         polydegree=3L){

  n <- as.integer(length(y))
  m <- as.integer(knot_bin[1])
  ones <- rep(1, n)
  if (!is.null(x_base)) x_base <- as.matrix(x_base)

  if(m == 0L){
    x <- cbind(ones, x_base)
  }else{
    knot.vec <- knot_bin[2:(m+1)]
    x <- splines::bs(1:n, degree=polydegree, knots = knot.vec)
    x <- cbind(ones, x_base, x)
  }

  ## Fastter fit: a thin wrapper to the "innermost" C code performing the
  ##              QR decomposition
  fit <- stats::.lm.fit(x, y)
  SSRes <- sum(fit$residuals^2)
  BIC_val <- n*log(SSRes/n) + (polydegree + 1L + m)*log(n)

  return(BIC_val)
}

#' BIC Objective for Fixed-Knots GA
#'
#' @description
#' Computes a BIC-style objective for a chromosome when the number of knots
#' is \emph{fixed} at \code{fixedknots}. The chromosome places those many
#' ordered knot locations in \code{knot_bin[2:(fixedknots+1)]}. A spline basis
#' of degree \code{polydegree} is built on \code{1:n}, with an intercept and
#' optional \code{x_base}, and the BIC is returned.
#'
#' @param knot_bin Numeric vector; GA chromosome with ordered knots.
#' @param plen Unused placeholder (kept for \code{changepointGA} compatibility).
#' @param y Numeric response vector.
#' @param x_base Optional matrix of additional covariates (will be column-bound).
#' @param fixedknots Integer; number of knots encoded in the chromosome.
#' @param polydegree Integer spline degree (default \code{3L}).
#'
#' @return Numeric scalar: BIC value (lower is better).
#'
#' @seealso \link{varyknotsBIC}, \link{gareg_knots}
#' @export
fixknotsBIC <- function(knot_bin, plen=0, y, x_base=NULL, fixedknots, polydegree=3L){

  n <- as.integer(length(y))
  ones <- rep(1, n)
  if (!is.null(x_base)) x_base <- as.matrix(x_base)

  knot.vec <- knot_bin[2:(fixedknots+1)]
  x <- splines::bs(1:n, degree=polydegree, knots = knot.vec)
  x <- cbind(ones, x_base, x)

  ## Fastter fit: a thin wrapper to the "innermost" C code performing the
  ##              QR decomposition
  fit <- stats::.lm.fit(x, y)
  SSRes <- sum(fit$residuals^2)
  BIC_val <- n*log(SSRes/n) + (polydegree + 1L + fixedknots)*log(n)

  return(BIC_val)
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
                         .env = parent.frame(), .validate = TRUE,
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

#' Crossover Operator (Fixed-Knots, Local Swap)
#'
#' @description
#' Attempts to inject mom-only knot locations into the child's vector while
#' preserving feasibility (\code{diff(child) >= minDist}) and using only parent
#' knots. Operates on fixed-\eqn{m} chromosomes.
#'
#' @param mom,dad Parent chromosomes (integer vectors).
#' @param prange Optional hyperparameter range (unused here).
#' @param minDist Integer minimum spacing constraint.
#' @param lmax Chromosome length.
#' @param N Series length (used to place \code{N+1} sentinel).
#'
#' @return Child chromosome vector (length \code{lmax}).
#'
#' @seealso \link{crossover_fixknots_2}, \link{mutation_fixknots}
#' @export
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
#' @seealso \link{crossover_fixknots}, \link{crossover_fixknots_2}
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
#' @noRd
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
#' @noRd
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
