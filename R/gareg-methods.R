#' Methods for displaying and summarizing `gareg` objects
#'
#' @title Show and summary methods for \code{gareg}
#' @description
#' \itemize{
#'   \item \code{show(object)}: Compact header with call, engine, and N.
#'   \item \code{summary(object, ...)}: GA settings (when available) and best solution.
#' }
#' @param object A \code{"gareg"} object.
#' @param ... Currently unused.
#' @return \code{show}: invisible NULL. \code{summary}: invisibly returns \code{object}.
#' @seealso \link{gareg-class}
#' @name gareg-methods
NULL

#' @keywords internal
#' @noRd
if (!isGeneric("summary")) {
  setGeneric("summary", function(object, ...) standardGeneric("summary"))
}

#' @rdname gareg-methods
#' @export
setMethod("show", "gareg", function(object) {
  hdr <- switch(
    object@method,
    subset    = "# Best Subset Variable Selection via GA                  #",
    varyknots = "# Varying Knots Detection via changepointGA              #",
    fixknots  = "# Fixed Knots Detection via changepointGA                #",
    "# GAReg Result                                                 #"
  )
  cat("##########################################################\n")
  cat(hdr, "\n", sep = "")
  cat("##########################################################\n")
  cat("Call: "); print(object@call)
  cat("   gaMethod: ", object@gaMethod, "\n", sep = "")
  cat("N: ", object@N, "\n", sep = "")
  if (!is.null(object@gaFit)) cat("\nUse summary() for GA settings and best solution.\n")
  invisible(NULL)
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
    m_int <- as.integer(if (length(m)) m[1] else NA_integer_)
    tau   <- as.integer(bs)
    if (!is.na(m_int) && length(tau) >= m_int) tau <- tau[seq_len(m_int)]
    cat("   Fitness   =", format(x@bestFitness, digits = digits), "\n")
    cat("   m         =", if (is.na(m_int)) 0L else m_int, "\n")
    cat("   knots     =", if (length(tau)) paste(tau, collapse = " ") else "<none>", "\n")

  } else if (identical(x@method, "subset")) {
    fn <- x@featureNames
    if (is.numeric(bs) && length(fn) > 0L && length(bs) == length(fn) && all(bs %in% c(0,1))) {
      idx <- which(bs != 0); k <- length(idx); nm <- if (k) fn[idx] else character(0)
      cat("   Fitness =", format(x@bestFitness, digits = digits), "\n")
      cat("   k =", k, "\n")
      cat("   subset =", if (k) paste(nm, collapse = ", ") else "<none>", "\n")
    } else if (is.numeric(bs)) {
      idx <- sort(unique(as.integer(bs[is.finite(bs) & bs > 0])))
      k <- length(idx); nm <- if (length(fn) >= max(c(idx, 0))) fn[idx] else as.character(idx)
      cat("   Fitness =", format(x@bestFitness, digits = digits), "\n")
      cat("   k =", k, "\n")
      cat("   Subset Id =", format(x@bestsol, digits = digits), "\n")
      cat("   Best subset =", if (k) paste(nm, collapse = ", ") else "<none>", "\n")
    } else if (is.character(bs)) {
      nm <- unique(bs[nzchar(bs)]); k <- length(nm)
      cat("   Fitness =", format(x@bestFitness, digits = digits), "\n")
      cat("   k =", k, "\n")
      cat("   subset =", if (k) paste(nm, collapse = ", ") else "<none>", "\n")
    } else {
      cat("   Fitness =", format(x@bestFitness, digits = digits), "\n")
      cat("   bestsol =", paste(bs, collapse = " "), "\n")
    }
  } else {
    cat("   Fitness =", format(x@bestFitness, digits = digits), "\n")
    cat("   bestsol =", paste(bs, collapse = " "), "\n")
  }

  invisible(x)
}

#' @rdname gareg-methods
#' @export
setMethod("summary", "gareg", function(object, ...) {
  print_summary_gareg(object, ...)
  invisible(object)
})
