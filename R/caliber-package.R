#' caliber: monotone calibration of a digital measure to a multi-item PRO
#'
#' One latent trait \eqn{\eta \sim N(0, 1)} drives an ordinal multi-item PRO
#' (logistic graded response model, items scored 0 to 3) and a continuous
#' digital or performance measure \eqn{y}:
#' \deqn{\mathrm{logit}\, P(C_j \ge k \mid \eta) = a_j \eta + d_{jk}, \qquad
#'   y = \mu(\eta) + e, \quad e \sim N(0, \sigma(\eta)^2), \quad
#'   \log \sigma(\eta) = g_0 + g_1 \eta,}
#' with \eqn{\mu(\eta) = \sum_m \beta_m B_m(\eta)} a non-increasing cubic
#' B-spline (24 coefficients, 20 interior knots on \eqn{[-4, 4]}). Higher
#' \eqn{\eta} means more of the PRO construct (worse sleep, more mobility
#' difficulty); the digital measure falls as \eqn{\eta} rises.
#'
#' The curve is fitted by monotonic likelihood boosting of B-splines
#' (Leitenstorfer and Tutz 2007) with the AIC stop, in one of three ways that
#' differ only in how \eqn{\eta} is handled (see [fit_ms()]):
#' \describe{
#'   \item{MS-Joint}{concurrent calibration: GRM items and curve estimated
#'     together, \eqn{\eta} integrated out over \eqn{p(\eta \mid C, y)};}
#'   \item{MS-EAP}{EAP plug-in: GRM on the PRO alone, EAP scores treated as
#'     known;}
#'   \item{MS-Fixed}{fixed-item calibration: GRM on the PRO alone, items then
#'     fixed, \eqn{\eta} integrated out (the recommended primary analysis).}
#' }
#'
#' Throughout the code the digital measure is called `sqi` (the sleep-quality
#' index of the worked example); any continuous measure that falls with the
#' latent trait can be used (gait speed in the NHANES example).
#'
#' @keywords internal
#' @importFrom stats plogis dnorm rnorm runif optim optimize optimHess uniroot
#' @importFrom stats isoreg quantile sd median na.omit setNames complete.cases
#' @importFrom splines splineDesign
#' @importFrom statmod gauss.quad.prob
#' @importFrom grDevices hcl.colors
#' @importFrom utils head tail
#' @import ggplot2
#' @importFrom patchwork plot_layout plot_annotation
"_PACKAGE"

# ---- Constants --------------------------------------------------------------

#' The three estimation methods
#'
#' `ms_methods` are the short labels, `ms_method_labels` the long names and
#' `ms_method_colours` the colours used in all figures.
#' @format Named character vectors.
#' @export
ms_methods <- c("MS-Joint", "MS-EAP", "MS-Fixed")

#' @rdname ms_methods
#' @export
ms_method_labels <- c("MS-Joint" = "Concurrent calibration",
                      "MS-EAP"   = "EAP plug-in",
                      "MS-Fixed" = "Fixed-item calibration")

#' @rdname ms_methods
#' @export
ms_method_colours <- c("MS-Joint" = "#2a78d6", "MS-EAP" = "#eb6834",
                       "MS-Fixed" = "#1baf7a")

ms_mode <- c("MS-Joint" = "joint", "MS-EAP" = "eap", "MS-Fixed" = "fixed")

# eta is integrated out over Q = 61 Gauss-Hermite nodes for N(0, 1)
quad    <- statmod::gauss.quad.prob(61, dist = "normal")
nodes   <- quad$nodes
log_wts <- log(quad$weights)

#' Quadrature nodes and log weights for N(0, 1)
#'
#' The 61 Gauss-Hermite nodes over which the latent trait is integrated out.
#' @return list(nodes, log_wts)
#' @export
ms_quadrature <- function() list(nodes = nodes, log_wts = log_wts)

#' Names of the thresholds of an item set
#'
#' Thresholds are stored item-major within category: C1_1, C2_1, ..., C7_1,
#' C1_2, ..., C7_3, matching `c(b)` for an items x 3 matrix `b`.
#' @param items item names
#' @return character vector of length 3 * length(items)
#' @export
threshold_names <- function(items)
  paste0(rep(items, 3), "_", rep(1:3, each = length(items)))

row_max <- function(L) L[cbind(seq_len(nrow(L)), max.col(L, ties.method = "first"))]
row_log_sum_exp <- function(L) { m <- row_max(L); m + log(rowSums(exp(L - m))) }
row_posterior   <- function(L) { p <- exp(L - row_max(L)); p / rowSums(p) }

utils::globalVariables(c("x", "xmin", "xmax", "y", "cat", "alpha", "SQI", "score",
                         "item", "source"))
