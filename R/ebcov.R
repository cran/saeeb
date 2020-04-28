#' @title EB Poisson-Gamma with Covariates
#'
#' @description This function gives the area level EB and MSE estimator based on Wakefield (2006) model and the refinement model by Kismiantini (2007).
#'
#' @param formula an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
#' @param data a mandatory data frame containing the variables in \code{formula} and \code{e}.
#' @param e a variable that contains the expected numbers of cases in each area.
#'
#' @details A typical model has the form response ~ terms where the response is a vector with numeric data type and terms is a set(s) of auxiliary variables.
#' @details A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link[stats]{formula}} for more details of allowed formulae.
#' @details The \code{formula} only accomodates variables with count data type and will be modeled using binomial negatif linear regression.
#'
#' @return The function returns a list with the following objects:
#' \describe{
#'  \item{EB}{data frame with number of rows equal to number of areas containing the EB estimator. For domains with zero sample size, the EB estimators are based on the synthetic regression.}
#'  \item{Parameter}{ }
#'  \itemize{
#'    \item alpha: the scale parameter estimator in Gamma distribution
#'    \item v: the shape parameter estimator in Gamma distribution
#'  }
#'  \item{fit}{ }
#'  \itemize{
#'    \item Estimate: maximum likelihood estimator of the model parameters
#'    \item SE: asymptotic estimate of the standard error of the the parameters
#'    \item Z: the Z statistic of the asymptotic hypothesis test that the population value for the parameter is 0
#'    \item LCL: lower 95\% confidence interval for the parameter estimators
#'    \item UCL: upper 95\% confidence interval for the parameter estimators
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
#' results <- ebcov(Y ~ AFF, lip, E)
#' results
#'
#' @seealso \code{\link[COUNT]{ml.nb2}}, \code{\link[MASS]{glm.nb}},
#'
#' @export

ebcov <- function(formula, data, e){
  result <- list(EB = NA, Parameter = list(alpha = NA, v = NA), fit = list(Estimate = NA, SE = NA, Z = NA, LCL = NA, UCL = NA), MSE.EB = list(method = "Jackknife", mse = NA), direct = list(est = NA, mse = NA))
  namevar <- deparse(substitute(e))
  if (!missing(data)) {
    formuladata <- model.frame(formula, na.action = na.omit, data)
    X <- model.matrix(formula, data)
    e <- data[ ,namevar]
    e <- as.matrix(e)
  }
  else {
    formuladata <- model.frame(formula, na.action = na.omit)
    X <- model.matrix(formula)
  }
  X <- X[,-1]
  X <- as.matrix(X)
  aux <- ncol(X)+1
  y <- formuladata[, 1]
  if (attr(attributes(formuladata)$terms, "response") == 1)
    textformula <- paste(formula[2], formula[1], formula[3])
  else textformula <- paste(formula[1], formula[2])
  if (length(na.action(formuladata)) > 0)
    stop("Argument formula=", textformula, " contains NA values.")
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
  result$Parameter$alpha <- alpha_cap
  result$Parameter$v <- v_cap
  regresi <- ml.nb2(formula, data, offset = log(e))
  result$fit$Estimate <- regresi$Estimate
  result$fit$SE <- regresi$SE
  result$fit$Z <- regresi$Z
  result$fit$LCL <- regresi$LCL
  result$fit$UCL <- regresi$UCL
  alpha <- matrix(c(rep(alpha_cap,m)), ncol = 1)
  X <- cbind(alpha, X)
  Bduga <- matrix(c(regresi$Estimate[1:aux]), nrow = aux)
  alpha_cov <- 1/regresi$Estimate[aux+1]
  mu_cap <- exp(X %*% Bduga)
  theta_i_cov <- y/(e*mu_cap)
  gamma_i <- e*mu_cap/(alpha_cov+e*mu_cap)
  EB <- gamma_i*theta_i_cov+(1-gamma_i)*mu_cap
  g1 <- (y+alpha_cov)/((e*mu_cap+alpha_cov)^2)
  result$EB <- EB
  jackknife <- function(y, e, j){
    X <- X[,-1]
    direct.est <- y[-j]/e[-j]
    e. <- mean(e[-j])
    theta_e._cap <- sum((e[-j]/e.)*direct.est)/m
    var_e <- sum((e[-j]/e.)*(direct.est-theta_e._cap)^2)/m
    alpha_cap <- theta_e._cap/(var_e-(theta_e._cap/e.))
    alpha <- matrix(c(rep(alpha_cap,m)), ncol = 1)
    X <- cbind(alpha, X)
    mu_cap <- exp(X %*% Bduga)
    theta_i_cov <- y/(e*mu_cap)
    gamma_i <- e*mu_cap/(alpha_cov+e*mu_cap)
    EB <- gamma_i*theta_i_cov+(1-gamma_i)*mu_cap
    g1 <- (y+alpha_cov)/((e*mu_cap+alpha_cov)^2)
    result <- list(EB = EB, g1 = g1)
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
