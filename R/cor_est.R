kendall_tau_multi <- function(p1, p2, rho){
  C_rho <- function(p1u, p2v){
    if(p1u == 0 | p2v == 0){
      return(0)
    }
    if(p1u == 1){
      return(p2v)
    } else {
      if(p2v == 1){
        return(p1u)
      }
      else {
        return(pbv::pbvnorm(qnorm(p1u), qnorm(p2v), rho))
      }
    }
  }

  tau <- 0
  k1 <- length(p1) - 1
  k2 <- length(p2) - 1
  C_rho_mat <- array(dim = c(k1+1, k2 + 1))
  for(u in 1:(k1 + 1)){
    for(v in 1:(k2 + 1)){
      C_rho_mat[u, v] <- C_rho(p1[u], p2[v])
    }
  }

  for(u in 1:k1){
    for(v in 1:k2){
      C1 <- C_rho_mat[u+1, v+1]
      C2 <- C_rho_mat[u+1, v]
      C3 <- C_rho_mat[u,v+1]
      C4 <- C_rho_mat[u, v]
      tau <- tau + (C1 + C4)^2 - (C2 + C3)^2
    }
  }
  return(tau - 1)
}


proxy_rho <- function(p1, p2, tau){
  k1 <- length(p1) - 1
  k2 <- length(p2) - 1

  d1 <- p1[3:(k1+1)] - p1[1:(k1 - 1)]
  d2 <- p2[3:(k2+1)] - p2[1:(k2 - 1)]
  A <- dnorm(qnorm(p1[2:k1]))
  B <- dnorm(qnorm(p2[2:k2]))

  a <- 0
  for (u in 1:(k1 - 1)){
    for (v in 1: (k2 - 1)){
      a <- a + d1[u] * d2[v] * A[u] * B[v]
    }
  }

  return(tau/(2*a))
}

dcoca_corr_rho_solving <- function(mapped_data, probs, nearPDproj = F, ncores = 2){
  #browser()
  res <- list()
  # probs <- list()
  p <- ncol(mapped_data)
  n <- nrow(mapped_data)

  # change this chunk for a hierachical version marginal distribtribution
  # for (j in 1:p){
  #   probs[[j]] <- c(0, cumsum(table(mapped_data[,j])/n))
  # }
  ## end
  # data_c <- gdtrans(mapped_data)$transformed_data
  # data_normal <- qnorm(data_c, mean = 0, sd = 1)
  # res$mapped_data <- mapped_data
  # res$data_c <- data_c
  # res$probs <- probs
  kendalltau <- fast_kendall(mapped_data)
  res$kendalltau <- kendalltau
  sigmahat_tau <- bigstatsr::FBM(p, p, init = 1)
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  foreach(i = 2:p, .combine = 'c', .export = c('kendall_tau_multi', 'proxy_rho')) %:%
    foreach(j = 1:(i-1), .combine = 'c', .export = c('kendall_tau_multi', 'proxy_rho'))%dopar%{
      # browser()
      tau_rho_equi <- function(x) {

        kendall_tau_multi(p1 = probs[[i]], p2 = probs[[j]], x) - kendalltau[i,j]
      }
      # tau_rho_fixed_point <- function(x) {
      #
      #   - kendall_tau_multi(p1 = probs[[i]], p2 = probs[[j]], x) - kendalltau[i,j] + x
      # }
      margin_epsilon <- 0.001
      epsilon_tau <- 1e-8
      if (abs(kendalltau[i,j]) < epsilon_tau){
        sigmahat_tau[i,j] <- 0}
      else{
        if (sign(tau_rho_equi(-1+margin_epsilon)) * sign(tau_rho_equi(1 - margin_epsilon)) < 0){
          init_rho <- proxy_rho(probs[[i]], probs[[j]], kendalltau[i,j])
          if (abs(init_rho) > (1 - margin_epsilon)){
            init_rho <- sign(init_rho)*0.5
          }
          finit <- tau_rho_equi(init_rho)
          if (abs(finit) < 1e-6){
            sigmahat_tau[i,j] <- init_rho
          } else {
            sign_init <- sign(finit)
            if (sign(init_rho) > 0 & sign_init > 0) interval_rho <- c(0, init_rho)
            if (sign(init_rho) > 0 & sign_init < 0) interval_rho <- c(init_rho, 1 - margin_epsilon)
            if (sign(init_rho) < 0 & sign_init > 0) interval_rho <- c(-1 + margin_epsilon, init_rho)
            if (sign(init_rho) < 0 & sign_init < 0) interval_rho <- c(init_rho, 0)
            sigmahat_tau[i,j] <- uniroot(tau_rho_equi, interval = interval_rho)$root
          }}
        else {
          if (sign(tau_rho_equi(-1+margin_epsilon)) > 0){
            sigmahat_tau[i,j] <- -1
          }
          if (sign(tau_rho_equi(1 - margin_epsilon)) < 0){
            sigmahat_tau[i,j] <- 1
          }
        }}

      sigmahat_tau[j, i] <- sigmahat_tau[i,j]
    }
  tau_corr <- sigmahat_tau[]
  if (nearPDproj){
    tau_corr <- Matrix::nearPD(tau_corr)$mat
  }
  res$tau_corr <- as.matrix(tau_corr)
  parallel::stopCluster(cl)

  return(res)
}





#' Fit discrete multivariate data using Gaussian copula
#'
#' @param mapped_data A discrete dataset whose discrete values are all numeric or have been mapped to numbers.
#' @param Fhat A list of estimated marginal cdfs of the dataset. Each item in the list is a \code{stepfun} object and the support of the cdfs (\code{knots(Fhats[[j]])}) should be the superset of the unique discrete values of the corresponding feature (\code{unique(mapped_data[,j])}). If not specified, \code{ecdf} will be used to estimate the marginal cdfs.
#' @param nearPDproj Logical. If TRUE, the estimated correlation matrix using Kendall's tau will be projected to the nearest positive definite matrix.
#' @param ncores Number of cores to use when parallelized computing of the correlation matrix.
#'
#' @return An object of class \code{DGCFit} with components
#' \item{mapped_data}{The input \code{mapped_data}.}
#' \item{tau}{The Kendall's tau matrix of the input discrete numerical dataset.}
#' \item{corr}{The estimated correlation matrix using the Kendall's tau matrix.}
#' \item{cdfs}{The list of estimated marginal cdfs.}
#'
#' @export
DGCFit <- function(mapped_data, Fhats, nearPDproj = T, ncores = 2,  ...){
  if (missing(Fhats)){
    Fhats <- Fhat(mapped_data)
  }

  p <- ncol(mapped_data)
  n <- nrow(mapped_data)

  probs <- list()

  for (j in 1:p){
    probs[[j]] <- c(0, Fhats[[j]](knots(Fhats[[j]])))
  }

  cor_est_res <- dcoca_corr_rho_solving(mapped_data, probs, nearPDproj = nearPDproj, ncores = ncores)

  res <- list()

  res$mapped_data <- mapped_data
  res$tau <- cor_est_res$kendalltau
  res$corr <- cor_est_res$tau_corr
  res$cdfs <- Fhats

  class(res) <- "DGCFit"

  return(res)

}
