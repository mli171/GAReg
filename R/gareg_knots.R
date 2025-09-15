#' GA-based spline knot selection
#'
#' @description
#' Selects interior spline knots on a time/grid variable `t` via a Genetic Algorithm.
#' Supports fixed `M` or a set of candidate M's (e.g., 0:6). Fits a GLM:
#'   y ~ Intercept + X + Spline(t; knots)
#'
#' @param y Numeric vector (response).
#' @param t Numeric vector (time/grid), same length as `y`.
#' @param X Optional matrix/data.frame of additional covariates (defaults to \code{NULL}).
#' @param family A \code{stats::family} object (default \code{gaussian()}).
#' @param method fixed knots or varying knots.
#' @param M Integer or integer vector. If length = 1, we apply the fixed knots selection.
#'          If length > 1, we fit each M and pick the best by `criterion`. Default NULL.
#' @param natural Use natural cubic spline basis via \code{splines::ns}. If FALSE, use \code{splines::bs}. Default TRUE.
#' @param degree Polynomial degree for \code{bs()} (ignored for \code{ns()}, which is cubic). Default 3.
#' @param minspace Minimum spacing between interior knots. Default 1.
#' @param popSize,maxiter,run GA hyperparameters passed to \code{GA::ga}.
#' @param seed RNG seed for reproducibility.
#' @param parallel Passed to \code{GA::ga} (TRUE recommended).
#'
#'
#' @return A list with:
#' \itemize{
#'   \item \code{best$M} chosen number of knots
#'   \item \code{best$knots} selected interior knot locations
#'   \item \code{best$model} fitted \code{glm} object at the optimum
#'   \item \code{best$score} minimized criterion
#'   \item \code{ga} the GA result for the winning M
#'   \item \code{scores_by_M} data.frame of scores for each tried M
#' }
#' @examples
#' \dontrun{


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

gareg_knots()


#' }
#' @export

gareg_knots = function(y,
                       t,
                       X=NULL,
                       family,
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

  ga_name <- if (is.function(gaMethod)) deparse(substitute(gaMethod)) else as.character(gaMethod)
  ga_fun  <- if (is.function(gaMethod)) gaMethod else get(gaMethod, mode = "function")
  engine  <- if (tolower(ga_name) == "cptgaisl") "cptgaisl" else "cptga"

  cptgactrl <- if (is.null(cptgactrl)) {
    cptgaControl(engine = engine)
  } else if (inherits(cptgactrl, "cptgaControl")) {
    cptgactrl
  } else if (is.list(cptgactrl)) {
    cptgaControl(.list = cptgactrl, engine = engine)
  } else stop("'cptgactrl' must be cptgaControl() or a named list.")

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

  regMethod_str <- if (missing(family)) "unspecified" else as.character(substitute(family))

  object <- new("gareg",
                call        = call,
                regMethod   = regMethod_str,
                N           = n,
                objFunc     = ObjFunc,
                gaMethod    = ga_name,
                gaFit       = GA.res,
                fixedknots  = if (is.null(fixedknots)) numeric() else as.numeric(fixedknots),
                minDist     = as.numeric(minDist),
                polydegree  = as.numeric(polydegree))

  return(object)

}



varyknotsBIC <- function(knot.bin,
                         plen=0,
                         y,
                         x_base=NULL,
                         polydegree=3L){

  n <- as.integer(length(y))
  m <- as.integer(knot.bin[1])
  ones <- rep(1, n)
  if (!is.null(x_base)) x_base <- as.matrix(x_base)

  if(m == 0L){
    x <- cbind(ones, x_base)
  }else{
    knot.vec <- knot.bin[2:(m+1)]
    x <- splines::bs(1:n, degree=polydegree, knots = knot.vec)
    x <- cbind(ones, x_base, x)
  }

  ## Fastter fit: a thin wrapper to the "innermost" C code performing the
  ##              QR decomposition
  fit <- .lm.fit(x, y)
  SSRes <- sum(fit$residuals^2)
  BIC.val <- n*log(SSRes/n) + (polydegree + 1L + m)*log(n)

  return(BIC.val)
}


fixknotsBIC <- function(knot.bin, plen=0, y, x_base=NULL, fixedknots, polydegree=3L){

  n <- as.integer(length(y))
  ones <- rep(1, n)
  if (!is.null(x_base)) x_base <- as.matrix(x_base)

  knot.vec <- knot.bin[2:(fixedknots+1)]
  x <- splines::bs(1:n, degree=polydegree, knots = knot.vec)
  x <- cbind(ones, x_base, x)

  ## Fastter fit: a thin wrapper to the "innermost" C code performing the
  ##              QR decomposition
  fit <- .lm.fit(x, y)
  SSRes <- sum(fit$residuals^2)
  BIC.val <- n*log(SSRes/n) + (polydegree + 1L + fixedknots)*log(n)

  return(BIC.val)
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

cptgaControl <- function(..., .list = NULL, .persist = FALSE,
                         .env = parent.frame(), .validate = TRUE,
                         engine = c("cptga","cptgaisl")) {
  engine <- match.arg(engine)
  defaults_name <- if (engine=="cptga") ".cptga.default" else ".cptgaisl.default"
  if (!exists(defaults_name, envir = .env, inherits = TRUE))
    stop(sprintf("`%s` not found in target environment.", defaults_name))
  current <- get(defaults_name, envir = .env, inherits = TRUE)

  overrides <- list(...)
  if (!is.null(.list)) {
    if (!is.list(.list)) stop("`.list` must be a named list.")
    overrides <- c(overrides, list(.list))
  }
  if (length(overrides) == 1L && is.list(overrides[[1]]) && !length(names(overrides)))
    overrides <- overrides[[1]]

  if (!length(overrides)) return(structure(current, class=c("cptgaControl","list")))

  if (is.null(names(overrides)) || any(names(overrides)==""))
    stop("All control overrides must be *named*.")

  unknown <- setdiff(names(overrides), names(current))
  if (length(unknown))
    stop("Unknown control parameter(s): ", paste(unknown, collapse = ", "))

  updated <- current
  updated[names(overrides)] <- overrides
  if (.validate) updated <- .validate_ctrl(updated, engine = engine)

  if (.persist) {
    assign(defaults_name, updated, envir = .env)
    invisible(structure(updated, class=c("cptgaControl","list")))
  } else {
    structure(updated, class=c("cptgaControl","list"))
  }
}

Popinitial_fixknots <- function(popSize, prange=NULL, N, minDist, Pb, mmax, lmax, fixedknots){

  pop <- matrix(0, nrow=lmax, ncol=popSize)

  for(j in 1:popSize){
    pop[,j] = selectTau_uniform_exact(N, fixedknots, minDist, lmax)
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

mutation_fixknots <- function(child, p.range = NULL, minDist, Pb, lmax, mmax, N) {

  m <- child[1]
  childMut <- selectTau_uniform_exact(N, m, minDist, lmax)

  return(childMut)
}


.cptga.default <- list(popSize = 200,
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



# if (!is.function(gaMethod)) {gaMethod <- get(gaMethod)}
#
# if(!is.null(fixedknots)){
#   # fixedknots
#   if(!is.null(ObjFunc)){
#     if(monitoring){cat("\n Default Objective Function in Use ...")}
#     ObjFunc <- fixknotsBIC
#     object@objFunc <- ObjFunc
#   }else{
#     if(monitoring){cat("\n Self-defined Objective Function in Use ...")}
#     ObjFunc=get(ObjFunc)
#     if (!is.function(ObjFunc)) {
#       stop("A fitness function must be provided")
#     }
#     object@objFunc <- ObjFunc
#   }
#   GA.res <- gaMethod(ObjFunc=ObjFunc,
#                      N=n,
#                      minDist=minDist,
#                      popInitialize="Popinitial_fixknots",
#                      selection="selection_linearrank",
#                      crossover="crossover_fixknots_2",
#                      mutation="mutation_fixknots",
#                      y=y,
#                      degree=degree,
#                      fixedknots=fixedknots)
#   object@gaFit <- GA.res
# }else{
#   # varyknots
#   if(!is.null(ObjFunc)){
#     if(monitoring){cat("\n Default Objective Function in Use ...")}
#     ObjFunc=varyknotsBIC
#     object@objFunc <- ObjFunc
#   }else{
#     if(monitoring){cat("\n Self-defined Objective Function in Use ...")}
#     ObjFunc=get(ObjFunc)
#     if (!is.function(ObjFunc)) {
#       stop("A fitness function must be provided")
#     }
#     object@objFunc <- ObjFunc
#   }
#   GA.res <- gaMethod(ObjFunc=ObjFunc,
#                      N=n,
#                      minDist=minDist,
#                      y=y,
#                      polydegree=polydegree)
#   object@gaFit <- GA.res
# }
# knotsEst <- GA.res@overbestchrom[2:(1+GA.res@overbestchrom[1])]







