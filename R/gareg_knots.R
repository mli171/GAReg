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
#' @param M Integer or integer vector. If length > 1, we fit each M and pick the best by `criterion`. Default 0:6.
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
#'   set.seed(1)
#'   n <- 150; t <- sort(runif(n, 0, 10))
#'   f <- function(x) 1 + 0.8*x - 0.5*pmax(0, x-3) + 0.6*pmax(0, x-7)
#'   y <- f(t) + rnorm(n, 0, 0.4)
#'   out <- gareg_knots(y, t, X = NULL, family = gaussian(), M = 0:5)
#'   out$best$M
#'   out$best$knots
#' }
#' @export
gareg_knots = function(y, t, X=NULL, family, method="vary", minspace=1, seed=NULL){



}
