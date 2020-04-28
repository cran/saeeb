#' Small Area Estimation for Count Data
#'
#' This package provides functions for small area estimation using Empirical Bayes (EB) Poisson-Gamma model. This model only accomodates count data type and gives option whether to use covariates in the estimation or not. Each function returns EB estimators and mean squared error (MSE) estimators for each area. The EB estimators are obtained using the model proposed by Wakefield (2006) and refined by Kismiantini (2007) and the MSE estimators are obtained using Jackknife method by Jiang et. al. (2002).
#'
#' @section Functions:
#' \describe{
#'  \item{\code{\link{ebcov}}}{Gives the EB Poisson-Gamma with covariates and the Jackknife MSE estimators.}
#'  \item{\code{\link{ebnocov}}}{Gives the EB Poisson-Gamma without covariates and the Jackknife MSE estimators.}
#' }
#'
#' @docType package
#' @name saeeb
#' @author Rizki Ananda Fauziah, Ika Yuni Wulansari
#' @references Clayton, David & Kaldor, John. (1987). Empirical Bayes Estimates of Age-Standardized Relative Risks for Use in Disease Mapping. Biometrics, 43, 671-681. doi:10.2307/2532003.
#' @references Jiang, J., Lahiri, P., & Wan, S. M. (2002). A Unified Jackknife Theory for Empirical Best Prediction with M-Estimation. The Annals of Statistics, 30, 6, 1782-1810. doi:10.1214/aos/1043351257.
#' @references Kismiantini. (2007). Pendugaan Statistik Area Kecil Berbasis Model Poisson-Gamma [Tesis]. Bogor: Institut Pertanian Bogor.
#' @references Rao, J. N. K. & Molina, Isabel. (2015). Small Area Estimation (2nd ed.). New Jersey: John Wiley & Sons, Inc.
#' @references Wakefield, Jon. (2006). Disease Mapping and Spatial Regression with Count Data. Biostatistics, 8, 2, 158–183. doi:10.1093/biostatistics/kxl008.
NULL
