
gareg_subset = function(y,
                        t,
                        X=NULL,
                        ObjFunc=NULL,
                        minDist=1,
                        gaMethod="cptga",
                        cptgactrl=NULL,
                        monitoring=FALSE,
                        seed=NULL,
                        ...){

  call <- match.call()
  dots <- list(...)
  n <- length(y)

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
    if (monitoring) cat("\nDefault Objective Function (fixknotsBIC) in use ...")
      ObjFunc <-
    } else {
    if (monitoring) cat("\nSelf-defined Objective Function in use ...")
    if (is.character(ObjFunc)) ObjFunc <- get(ObjFunc, mode = "function")
    if (!is.function(ObjFunc)) stop("ObjFunc must be a function or a function name.")
  }






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

  mhat <- object@bestnumbsol <- object@bestChrom[1]
  if(mhat == 0){
    object@bestsol <- NULL
  }else{
    object@bestsol <- object@bestChrom[2:(1+mhat)]
  }

  return(object)

}

#' BIC Objective for Best Subset GA
#'
#' @description
#' Computes a BIC-style objective for a chromosome encoding
#'
#'
#'
#'
#'
#'
#' @export
subsetBIC <- function(subset_bin,
                      plen=0,
                      y,
                      X) {

  n   <- as.integer(length(y))
  X   <- as.matrix(X)
  m   <- as.integer(subset_bin[1])
  K   <- subset_bin[2+m] - 1

  if (m == 0L) {
    rss <- sum((y - mean(y))^2)
    k_eff <- 1L
  }else{
    idX <- subset_bin[2:(1+m)]
    Xsel <- X[, idX, drop = FALSE]
    fit <- stats::.lm.fit(cbind(1, Xsel), y)
    if (!is.null(fit$rank) && fit$rank < (m + 1L)) return(Inf)
    rss <- sum(fit$residuals^2)
    k_eff <- m + 1L
  }
  BIC_val <- n * log(rss/n) + k_eff * log(n)

  return(BIC_val)

}
