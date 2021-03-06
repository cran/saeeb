#' Lip Cancer in Scotland
#'
#' This dataset sets out observed and "expected" cases of lip cancer registered during the 6 years from 1975 to 1980 in each of the 56 counties of Scotland. These are the districts prior to the 1995 reorganization of local government. The dataset includes district names and identifying numbers and for district \code{i} with \code{i = 1, ..., 56}: the number of observed cases \emph{Yi}; the number of expected cases \emph{Ei}; and the value of a single covariate (percent of population employed in agriculture, fishing and forestry).
#'
#' @format A data frame with 56 rows and 5 variables:
#' \describe{
#'   \item{ID}{The district identifying number.}
#'   \item{district.name}{The district name.}
#'   \item{Y}{The number of observed lip cancer cases.}
#'   \item{E}{The number of expected lip cancer cases.}
#'   \item{AFF}{The percentage of population employed in agriculture, fishing and forestry.}
#' }
"lip"
