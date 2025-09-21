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
setClassUnion("listOrNULL", members = c("list", "NULL"))
setClassUnion("functionOrNULL", c("function","NULL"))
setClassUnion("cptgaORcptgaislORNULL", c("cptga", "cptgaisl", "NULL"))

#' @rdname gareg-class
#' @export
setClass(
  Class = "gareg",
  slots = c(
    call        = "language",
    method      = "character",
    regMethod   = "character",
    family      = "character",
    N           = "numeric",
    objFunc     = "functionOrNULL",
    gaMethod    = "character",
    gaFit       = "cptgaORcptgaislORNULL",
    ctrl        = "listOrNULL",
    # knots
    fixedknots  = "numericOrNULL",
    minDist     = "numeric",
    polydegree  = "numericOrNULL",
    # best subset
    subsetSpec  = "listOrNULL",
    featureNames= "character",
    bestFitness = "numeric",
    bestChrom   = "numeric",
    bestsol     = "numeric"
  ),
  prototype = list(
    method       = "varyknots",
    regMethod    = "unspecified",
    family       = "unspecified",
    fixedknots   = numeric(),
    minDist      = NA_real_,
    polydegree   = NA_real_,
    subsetSpec   = NULL,
    featureNames = character(),
    ctrl         = NULL,
    bestFitness  = NA_real_,
    bestChrom    = numeric(),
    bestsol      = NULL
  ),
  package = "GAReg"
)

setMethod("print", "gareg", function(x, ...) str(x))

.header_from_method <- function(m) switch(
  m,
  subset    = "# Best Subset Variable Selection via GA       #",
  varyknots = "# Varying Knots Detection via GA              #",
  fixknots  = "# Fixed Knots Detection via GA                #",
  "# GAReg Result                                      #"
)

.s <- function(x, nm, default = NA) {
  if (is.null(x)) return(default)
  if (methods::hasSlot(x, nm)) slot(x, nm) else default
}

setMethod("show", "gareg", function(object) {
  cat("###############################################\n")
  cat(.header_from_method(object@method), "\n", sep = "")
  cat("###############################################\n")
  cat("Call: "); print(object@call)
  cat("regMethod: ", object@regMethod,
      "   family: ", object@family,
      "   gaMethod: ", object@gaMethod, "\n", sep = "")
  cat("N: ", object@N, "\n", sep = "")
  if (!is.null(object@gaFit)) cat("\nUse summary() for GA settings and best solution.\n")
})


#' @export
print.summary.gareg <- function(x, digits = getOption("digits"), max_display = 5, ...) {
  gf <- x@gaFit
  cat("###############################################\n")
  cat(.header_from_method(x@method), "\n", sep = "")
  cat("###############################################\n")

  cat("   Settings: \n")
  cat("   Population size         = ", .s(gf, "popSize", "NA"), "\n", sep = "")
  cat("   Number of generations   = ", .s(gf, "count", "NA"), "\n", sep = "")
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

  cat("\n##### Best (normalized) ##### \n")
  if (is.null(x@best)) {
    cat("   <no normalized best stored>\n")
  } else if (identical(x@best$type, "changepoint")) {
    m <- x@best$m; tau <- x@best$knots
    cat("   Fitness =", format(x@bestFitness, digits=digits), "\n")
    cat("   m =", m, "\n")
    cat("   knots =", if (length(tau)) paste(tau, collapse=" ") else "<none>", "\n")
  } else if (identical(x@best$type, "subset")) {
    vv <- x@best$vars; k <- x@best$k
    cat("   Fitness =", format(x@bestFitness, digits=digits), "\n")
    cat("   k =", k, "\n")
    cat("   subset =", if (length(vv)) paste(vv, collapse=", ") else "<none>", "\n")
  } else {
    cat("   type =", x@best$type, " (custom)\n")
    cat("   Fitness =", format(x@bestFitness, digits=digits), "\n")
  }

  invisible(x)
}

setMethod("summary", "gareg", function(object, ...) print.summary.gareg(object, ...))
