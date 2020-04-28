#' @title EB Poisson-Gamma without Covariates
#'
#' @description This function gives the area level EB and MSE estimator based on Clayton & Kaldor (1987).
#'
#' @param y a response variable contains the numbers of cases in each area.
#' @param data a mandatory data frame containing the \code{y} and \code{e} variables.
#' @param e a variable that contains the expected numbers of cases in each area.
#'
#' @details This function only accomodates variables with count data type.
#'
#' @return The function returns a list with the following objects:
#' \describe{
#'  \item{EB}{data frame with number of rows equal to number of areas containing the EB estimator. For domains with zero sample size, the EB estimators are based on the synthetic regression.}
#'  \item{Parameters}{ }
#'  \itemize{
#'    \item alpha: the scale estimator in Gamma distribution
#'    \item v: the shape parameter estimator in Gamma distribution
#'  }
#'  \item{MSE.EB}{ }
#'  \itemize{
#'    \item method: Jackknife
#'    \item mse: the mean squared error estimator of the EB estimators
#'  }
#'  \item{direct}{ }
#'  \itemize{
#'    \item est: direct estimators for the response variable
#'    \item mse: the mean squared error estimator of the direct estimators
#'  }
#' }
#'
#' @examples
#' #Load dataset
#' data(lip)
#'
#' #Save output as an object
#' results <- ebnocov(Y, lip, E)
#' results
#'
#' @export

ebnocov <- function(y, data, e){
  result <- list(EB = NA, Parameter = list(alpha = NA, v = NA), MSE.EB = list(method = "Jackknife", mse = NA), direct = list(est = NA, mse = NA))
  response <- deparse(substitute(y))
  namevar <- deparse(substitute(e))
  y <- data[ ,response]
  e <- data[ ,namevar]
  if (any(is.na(y)))
    stop("Argument y=", response, " contains NA values.")
  if (any(is.na(e)))
    stop("Argument e=", namevar, " contains NA values.")
  m <- length(y)
  direct.est <- y/e
  mse.dir <- direct.est/e
  result$direct$est <- direct.est
  result$direct$mse <- mse.dir
  e. <- mean(e)
  theta_e._cap <- sum((e/e.)*direct.est)/m
  var_e <- sum((e/e.)*(direct.est-theta_e._cap)^2)/m
  alpha_cap <- theta_e._cap/(var_e-(theta_e._cap/e.))
  v_cap <- theta_e._cap^2/(var_e-(theta_e._cap/e.))
  gamma_i <- e/(e+alpha_cap)
  EB <- gamma_i*direct.est+(1-gamma_i)*theta_e._cap
  g1 <- (y+v_cap)/((e+alpha_cap)^2)
  result$Parameter$alpha <- alpha_cap
  result$Parameter$v <- v_cap
  result$EB <- EB
  jackknife <- function(y, e, j){
    direct.est_jk <- 0
    e.jk <- 0
    theta_e._cap_jk <- 0
    direct.est_jk <- y[-j]/e[-j]
    e.jk <- mean(e[-j])
    theta_e._cap_jk <- sum((e[-j]/e.jk)*direct.est_jk)/m
    var_e <- sum((e[-j]/e.jk)*(direct.est_jk-theta_e._cap_jk)^2)/m
    alpha_cap <- theta_e._cap_jk/(var_e-(theta_e._cap_jk/e.jk))
    v_cap <- theta_e._cap_jk^2/(var_e-(theta_e._cap_jk/e.jk))
    gamma_i <- e/(e+alpha_cap)
    EB <- gamma_i*direct.est+(1-gamma_i)*theta_e._cap
    g1 <- (y+v_cap)/((e+alpha_cap)^2)
    result <- list(gamma = gamma_i, EB = EB, g1 = g1)
    return(result)
  }
  jk <- lapply(1:m, function(j) jackknife(y, e, j))
  M1 <- sapply(1:m, function(i){
    m1 <- g1[i]-(m-1)/m*sum(sapply(1:m, function(j){
      return(jk[[j]]$g1[i]-g1[i])
    }))
    return(m1)
  })
  M2 <- sapply(1:m, function(i){
    m2 <- ((m-1)/m)*sum(sapply(1:m, function(j){
      return((jk[[j]]$EB[i]-EB[i])^2)
    }))
    return(m2)
  })
  mse <- M1+M2
  result$MSE.EB$mse <- mse
  return(result)
}
