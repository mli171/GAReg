#' S4 Class for Genetic Algorithm-Based Regression
#'
#' @name gareg-class
#' @title S4 Class Definition for `gareg`
#' #'
#' @slot call The matched call that created the object.
#'
#' @slot N The sample size of the time series.
#' @importClassesFrom changepointGA cptga cptgaisl
NULL

setClassUnion("numericOrNULL", members = c("numeric", "NULL"))
setClassUnion("numericOrChara", members = c("numeric", "character"))
setClassUnion("listOrNULL", members = c("list", "NULL"))
setClassUnion("functionOrNULL", c("function","NULL"))
setClassUnion("gaBackendORNULL", c("cptga","cptgaisl","ga","gaisl","NULL"))

#' GAReg result container
#'
#' @description
#' S4 container for GA-based regression/changepoint tasks. Holds the GA
#' backend fit and a normalized summary of the best solution.
#'
#' @slot call language. The original call.
#' @slot method character. One of "varyknots", "fixknots", "subset".
#' @slot N numeric. Length of the response vector.
#' @slot objFunc functionOrNULL. Objective function used.
#' @slot gaMethod character. GA engine name ("cptga","cptgaisl","ga","gaisl").
#' @slot gaFit Backend GA fit object (union of classes from GA and changepointGA).
#' @slot ctrl listOrNULL. Control list used to run the GA.
#' @slot fixedknots numericOrNULL. Fixed knot locations (NULL or numeric()).
#' @slot minDist numeric. Minimum distance between adjacent changepoints.
#' @slot polydegree numericOrNULL. Spline degree for default objectives.
#' @slot subsetSpec listOrNULL. Constraints for subset selection (unused for knots).
#' @slot featureNames character. Candidate feature names (subset tasks).
#' @slot bestFitness numeric. Best fitness value found.
#' @slot bestChrom numeric. Raw best chromosome returned by the backend.
#' @slot bestnumbsol numeric. Count of selected elements (e.g., m for knots).
#' @slot bestsol numericOrCharacter. For knots: the m knot locations; for subset: mask/indices/names.
#' @seealso \link{gareg_knots}, \link{cptgaControl}
#' @exportClass gareg
setClass(
  Class = "gareg",
  slots = c(
    call        = "language",
    method      = "character",
    N           = "numeric",
    objFunc     = "functionOrNULL",
    gaMethod    = "character",
    gaFit       = "gaBackendORNULL",
    ctrl        = "listOrNULL",
    # knots
    fixedknots  = "numericOrNULL",
    minDist     = "numeric",
    polydegree  = "numericOrNULL",
    # best subset
    subsetSpec  = "listOrNULL",
    featureNames= "character",
    # general results
    bestFitness = "numeric",
    bestChrom   = "numeric",
    bestnumbsol = "numeric",
    bestsol     = "numericOrChara"
  ),
  prototype = list(
    method       = character(),
    fixedknots   = NA_real_,
    minDist      = numeric(),
    polydegree   = NA_real_,
    subsetSpec   = NULL,
    featureNames = character(),
    ctrl         = NULL,
    bestFitness  = numeric(),
    bestChrom    = numeric(),
    bestnumbsol  = numeric(),
    bestsol      = numeric()
  ),
  package = "GAReg"
)

setMethod("print", "gareg", function(x, ...) str(x))

# .header_from_method <- function(m) switch(
#   m,
#   subset    = "# Best Subset Variable Selection via GA                  #",
#   varyknots = "# Varying Knots Detection via changepointGA              #",
#   fixknots  = "# Fixed Knots Detection via changepointGA                #",
#   "# GAReg Result                                                 #"
# )

.s <- function(x, nm, default = NA) {
  if (nm %in% methods::slotNames(x)) methods::slot(x, nm) else default
}

setMethod("show", "gareg", function(object) {
  hdr <- switch(
    object@method,
    subset    = "# Best Subset Variable Selection via GA                  #",
    varyknots = "# Varying Knots Detection via changepointGA              #",
    fixknots  = "# Fixed Knots Detection via changepointGA                #",
    "# GAReg Result                                                 #"
  )

  cat("###############################################\n")
  # cat(.header_from_method(object@method), "\n", sep = "")
  cat(hdr, "\n", sep = "")
  cat("###############################################\n")
  cat("Call: "); print(object@call)
  cat("   gaMethod: ", object@gaMethod, "\n", sep = "")
  cat("N: ", object@N, "\n", sep = "")
  if (!is.null(object@gaFit)) cat("\nUse summary() for GA settings and best solution.\n")
})


#' @keywords internal
#' @noRd
print_summary_gareg <- function(x, digits = getOption("digits"), max_display = 5, ...) {
  gf <- x@gaFit
  cat("##########################################################\n")
  hdr <- switch(
    x@method,
    subset    = "# Best Subset Variable Selection via GA                  #",
    varyknots = "# Varying Knots Detection via changepointGA              #",
    fixknots  = "# Fixed Knots Detection via changepointGA                #",
    "# GAReg Result                                                 #"
  )
  cat(hdr, "\n", sep = "")
  cat("##########################################################\n")

  cat("   Settings: \n")
  cat("   Population size         = ", .s(gf, "popSize", "NA"), "\n", sep = "")
  cat("   Number of generations   = ", .s(gf, "count", "NA"), "\n", sep = "")
  cat("   GA convergence          = ", .s(gf, "convg", FALSE), "\n", sep = "")
  cat("   Crossover probability   = ", format(.s(gf, "pcrossover", NA), digits=digits), "\n", sep = "")
  cat("   Mutation probability    = ", format(.s(gf, "pmutation",  NA), digits=digits), "\n", sep = "")
  cat("   Changepoint probability = ", format(.s(gf, "pchangepoint",NA), digits=digits), "\n", sep = "")
  cat("   Parallel Usage          = ", .s(gf, "parallel", FALSE), "\n", sep = "")
  if (isTRUE(.s(gf, "parallel", FALSE))) {
    cat("   Number of thread        = ", .s(gf, "nCore", "NA"), "\n", sep = "")
  }
  sugg <- .s(gf, "suggestions", NULL)
  if (!is.null(sugg)) {
    cat("   Suggestions:\n")
    for (i in seq_along(sugg)) {
      cat("     [", i, "]: ", paste(sugg[[i]], collapse=" "), "\n", sep = "")
      if (i >= max_display) { cat("     ...\n"); break }
    }
  }

  cat("\n##### Best ##### \n")
  m  <- x@bestnumbsol
  bs <- x@bestsol

  if (length(bs) == 0L || all(is.na(bs))) {
    cat("   <no best solution stored>\n")

  } else if (x@method %in% c("varyknots","fixknots","cptdetect")) {
    ## Knots-style: bestnumbsol = m, bestsol = tau vector
    m_int <- as.integer(if (length(m)) m[1] else NA_integer_)
    tau   <- as.integer(bs)
    if (!is.na(m_int) && length(tau) >= m_int) {
      tau <- tau[seq_len(m_int)]
    }
    cat("   Fitness   =", format(x@bestFitness, digits = digits), "\n")
    cat("   m         =", if (is.na(m_int)) 0L else m_int, "\n")
    cat("   knots     =", if (length(tau)) paste(tau, collapse = " ") else "<none>", "\n")

  } else if (identical(x@method, "subset")) {
    ## Subset-style: bestsol is either a 0/1 mask, indices, or names
    fn <- x@featureNames

    if (is.numeric(bs) && length(fn) > 0L && length(bs) == length(fn) && all(bs %in% c(0,1))) {
      # 0/1 mask
      idx <- which(bs != 0)
      k   <- length(idx)
      nm  <- if (k) fn[idx] else character(0)
      cat("   Fitness =", format(x@bestFitness, digits = digits), "\n")
      cat("   k =", k, "\n")
      cat("   subset =", if (k) paste(nm, collapse = ", ") else "<none>", "\n")

    } else if (is.numeric(bs)) {
      # integer indices
      idx <- sort(unique(as.integer(bs[is.finite(bs) & bs > 0])))
      k   <- length(idx)
      nm  <- if (length(fn) >= max(c(idx, 0))) fn[idx] else as.character(idx)
      cat("   Fitness =", format(x@bestFitness, digits = digits), "\n")
      cat("   k =", k, "\n")
      cat("   Subset Id =", format(x@bestsol, digits = digits), "\n")
      cat("   Best subset =", if (k) paste(nm, collapse = ", ") else "<none>", "\n")

    } else if (is.character(bs)) {
      # feature names directly
      nm <- unique(bs[nzchar(bs)])
      k  <- length(nm)
      cat("   Fitness =", format(x@bestFitness, digits = digits), "\n")
      cat("   k =", k, "\n")
      cat("   subset =", if (k) paste(nm, collapse = ", ") else "<none>", "\n")

    } else {
      cat("   Fitness =", format(x@bestFitness, digits = digits), "\n")
      cat("   bestsol =", paste(bs, collapse = " "), "\n")
    }

  } else {
    ## Fallback: just show the vector
    cat("   Fitness =", format(x@bestFitness, digits = digits), "\n")
    cat("   bestsol =", paste(bs, collapse = " "), "\n")
  }

  invisible(x)
}

#' @export
setMethod("summary", "gareg", function(object, ...) print_summary_gareg(object, ...))
