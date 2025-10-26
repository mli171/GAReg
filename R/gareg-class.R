#' S4 Class for Genetic Algorithm-Based Regression
#'
#' @name gareg-class
#' @title S4 Class Definition for `gareg`
#'
#' @slot call The matched call that created the object.
#' @slot N The effective size of the x grid used for knot search
#'   (i.e., `length(x_unique)`), typically the number of unique `x`.
#' @importClassesFrom changepointGA cptga cptgaisl
#' @importFrom methods show slotNames slot
NULL
setClassUnion("numericOrNULL",   members = c("numeric", "NULL"))
setClassUnion("numericOrChara",  members = c("numeric", "character"))
setClassUnion("listOrNULL",      members = c("list", "NULL"))
setClassUnion("functionOrNULL",  members = c("function","NULL"))
setClassUnion("gaBackendORNULL", members = c("cptga","cptgaisl","ga","gaisl","NULL"))


#' GAReg result container
#'
#' @description
#' S4 container for GA-based regression/changepoint tasks. Holds the GA
#' backend fit and a normalized summary of the best solution.
#'
#' @slot call language. The original call.
#' @slot method character. One of "varyknots", "fixknots", "subset".
#' @slot N numeric. Length of `x_unique` used by the GA (also `sentinel-1`).
#' @slot objFunc functionOrNULL. Objective function used.
#' @slot gaMethod character. GA engine name ("cptga","cptgaisl","ga","gaisl").
#' @slot gaFit Backend GA fit object (union of classes from GA and changepointGA).
#' @slot ctrl listOrNULL. Control list used to run the GA (if stored by caller).
#' @slot fixedknots numericOrNULL. Fixed number of interior knots (`m`) for fixed-knots mode, or NULL.
#' @slot minDist numeric. Minimum distance between adjacent changepoints.
#' @slot polydegree numericOrNULL. Spline degree for default objectives.
#' @slot type character. One of `c("ppolys", "ns", "bs")` indicating piecewise
#'       polynomials, natural cubic, or B-spline.
#' @slot intercept logical. Whether the spline basis included an intercept column.
#' @slot subsetSpec listOrNULL. Constraints for subset selection (unused for knots).
#' @slot featureNames character. Candidate feature names (subset tasks).
#' @slot bestFitness numeric. Best fitness value found.
#' @slot bestChrom numeric. Raw best chromosome returned by the backend
#'       (may include a sentinel equal to `N+1` and optional padding).
#' @slot bestnumbsol numeric. Count of selected elements (e.g., `m` for knots).
#' @slot bestsol numericOrChara. For knots: the `m` interior indices (pre-sentinel);
#'       for subset: mask/indices/names.
#' @seealso [gareg_knots], [cptgaControl]
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
    type        = "character",
    intercept   = "logical",
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
    type         = character(),
    intercept    = logical(),
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


.s <- function(x, nm, default = NA) {
  if (nm %in% methods::slotNames(x)) methods::slot(x, nm) else default
}

# -------------------------------
# show() and summary() methods
# -------------------------------

#' @export
setMethod("show", "gareg", function(object) {
  hdr <- switch(
    object@method,
    subset    = "# Best Subset Variable Selection via GA                  #",
    varyknots = "# Varying-Knots Spline Regression via changepointGA     #",
    fixknots  = "# Fixed-Knots Spline Regression via changepointGA       #",
    "# GAReg Result                                                 #"
  )

  cat("##########################################################\n")
  cat(hdr, "\n", sep = "")
  cat("##########################################################\n")
  cat("Call: "); print(object@call)
  cat("   gaMethod: ", object@gaMethod, "\n", sep = "")
  cat("   N (|x_unique|): ", object@N, "\n", sep = "")
  if (object@method %in% c("varyknots","fixknots")) {
    cat("   Spline type / degree / intercept: ",
        paste0(object@type, " / ", object@polydegree, " / ", object@intercept), "\n", sep = "")
    cat("   minDist: ", object@minDist, "\n", sep = "")
  }
  if (!is.null(object@gaFit)) cat("\nUse summary() for GA settings and best solution.\n")
})

# internal pretty-printer used by summary()
#' @keywords internal
#' @noRd
print_summary_gareg <- function(x, digits = getOption("digits"), max_display = 5, ...) {
  gf <- x@gaFit
  cat("##########################################################\n")
  hdr <- switch(
    x@method,
    subset    = "# Best Subset Variable Selection via GA                  #",
    varyknots = "# Varying-Knots Spline Regression via changepointGA     #",
    fixknots  = "# Fixed-Knots Spline Regression via changepointGA       #",
    "# GAReg Result                                                 #"
  )
  cat(hdr, "\n", sep = "")
  cat("##########################################################\n")

  # --- Settings (best we can, across engines) ---
  cat("   Settings: \n")
  cat("   Population size         = ", .s(gf, "popSize", "NA"), "\n", sep = "")
  cat("   Number of generations   = ", .s(gf, "count", "NA"), "\n", sep = "")
  cat("   GA convergence          = ", .s(gf, "convg", FALSE), "\n", sep = "")
  cat("   Crossover probability   = ", format(.s(gf, "pcrossover", NA), digits=digits), "\n", sep = "")
  cat("   Mutation probability    = ", format(.s(gf, "pmutation",  NA), digits=digits), "\n", sep = "")
  cat("   Changepoint probability = ", format(.s(gf, "pchangepoint",NA), digits=digits), "\n", sep = "")
  cat("   Parallel Usage          = ", .s(gf, "parallel", FALSE), "\n", sep = "")
  if (isTRUE(.s(gf, "parallel", FALSE))) {
    cat("   Number of threads       = ", .s(gf, "nCore", "NA"), "\n", sep = "")
  }
  sugg <- .s(gf, "suggestions", NULL)
  if (!is.null(sugg)) {
    cat("   Suggestions:\n")
    for (i in seq_along(sugg)) {
      cat("     [", i, "]: ", paste(sugg[[i]], collapse=" "), "\n", sep = "")
      if (i >= max_display) { cat("     ...\n"); break }
    }
  }

  # --- Best solution ---
  cat("\n##### Best ##### \n")
  m  <- x@bestnumbsol
  bs <- x@bestsol

  if (length(bs) == 0L || all(is.na(bs))) {
    # try to decode from the raw chromosome if available (sentinel-aware)
    chrom <- x@bestChrom
    if (length(chrom)) {
      m_int <- as.integer(chrom[1L])
      Nu    <- as.integer(x@N)
      tail  <- as.integer(chrom[-1L])
      endp  <- match(Nu + 1L, tail)           # sentinel = N+1 (length(x_unique)+1)
      idx   <- if (!is.na(endp)) tail[seq_len(endp - 1L)] else integer()
      idx   <- idx[idx != 0L]
      if (!is.na(m_int) && length(idx) >= m_int) {
        idx <- idx[seq_len(m_int)]
      }
      bs <- idx
      m  <- m_int
    }
  }

  if (length(bs) == 0L || all(is.na(bs))) {
    cat("   <no best solution stored>\n")

  } else if (x@method %in% c("varyknots","fixknots","cptdetect")) {
    ## Knots-style: bestnumbsol = m, bestsol = indices into x_unique
    m_int <- as.integer(if (length(m)) m[1] else NA_integer_)
    tau   <- as.integer(bs)
    if (!is.na(m_int) && length(tau) >= m_int) {
      tau <- tau[seq_len(m_int)]
    }
    cat("   Fitness   = ", format(x@bestFitness, digits = digits), "\n", sep = "")
    cat("   m         = ", if (is.na(m_int)) 0L else m_int, "\n", sep = "")
    cat("   indices   = ", if (length(tau)) paste(tau, collapse = " ") else "<none>", "\n", sep = "")
    cat("   sentinel  = N+1 = ", as.integer(x@N) + 1L, "\n", sep = "")

  } else if (identical(x@method, "subset")) {
    ## Subset-style: bestsol is either a 0/1 mask, indices, or names
    fn <- x@featureNames

    if (is.numeric(bs) && length(fn) > 0L && length(bs) == length(fn) && all(bs %in% c(0,1))) {
      # 0/1 mask
      idx <- which(bs != 0)
      k   <- length(idx)
      nm  <- if (k) fn[idx] else character(0)
      cat("   Fitness = ", format(x@bestFitness, digits = digits), "\n", sep = "")
      cat("   k       = ", k, "\n", sep = "")
      cat("   subset  = ", if (k) paste(nm, collapse = ", ") else "<none>", "\n", sep = "")

    } else if (is.numeric(bs)) {
      # integer indices
      idx <- sort(unique(as.integer(bs[is.finite(bs) & bs > 0])))
      k   <- length(idx)
      nm  <- if (length(fn) >= max(c(idx, 0))) fn[idx] else as.character(idx)
      cat("   Fitness    = ", format(x@bestFitness, digits = digits), "\n", sep = "")
      cat("   k          = ", k, "\n", sep = "")
      cat("   Subset Id  = ", paste(format(x@bestsol, digits = digits), collapse = " "), "\n", sep = "")
      cat("   Best subset= ", if (k) paste(nm, collapse = ", ") else "<none>", "\n", sep = "")

    } else if (is.character(bs)) {
      # feature names directly
      nm <- unique(bs[nzchar(bs)])
      k  <- length(nm)
      cat("   Fitness = ", format(x@bestFitness, digits = digits), "\n", sep = "")
      cat("   k       = ", k, "\n", sep = "")
      cat("   subset  = ", if (k) paste(nm, collapse = ", ") else "<none>", "\n", sep = "")

    } else {
      cat("   Fitness = ", format(x@bestFitness, digits = digits), "\n", sep = "")
      cat("   bestsol = ", paste(bs, collapse = " "), "\n", sep = "")
    }
  } else {
    ## Fallback: just show the vector
    cat("   Fitness = ", format(x@bestFitness, digits = digits), "\n", sep = "")
    cat("   bestsol = ", paste(bs, collapse = " "), "\n", sep = "")
  }

  invisible(x)
}

#' @export
setMethod("summary", "gareg", function(object, ...) print_summary_gareg(object, ...))
