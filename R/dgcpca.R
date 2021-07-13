truncnormalsim_trans <- function(lbounds, ubounds, inity, rotmat, lambda, sig, iter = 1){
  d <- sqrt(lambda - sig^2)
  iteration <- 0
  inity_backup <- inity
  while (iteration < iter){
    # browser()
    x_means <- sweep(inity %*% rotmat, 2, FUN = "*", STATS = d/(sig *sqrt(d^2 + sig^2)))
    xstar <- apply(x_means, MARGIN = c(1,2), function(x) rnorm(n = 1, mean = x))
    xstar <- sweep(xstar, 2, FUN = "*", STATS = d*sig/sqrt(sig^2 + d^2)) %*% t(rotmat)
    bounds <- array(c(lbounds, ubounds), dim = c(dim(lbounds), 2))
    bounds <- sweep(bounds, MARGIN = c(1,2), FUN = "-", STAT = xstar)/sig
    qustar <- apply(bounds, MARGIN = c(1,2), function(x) TruncatedNormal::rtnorm(n = 1, mu = 0, sd = 1, lb = x[1], ub = x[2]))
    ystar <- sig * qustar + xstar
    iteration <- iteration + 1
    inity <- ystar
  }
  # dim_ystar <- dim(ystar)
  # browser()
  ystar[ystar == Inf] <- inity_backup[ystar == Inf]
  ystar[ystar == -Inf] <- inity_backup[ystar == -Inf]

  return(ystar)
}

firstnpc <- function(npc, Rtilda, mapped_data, Fhats){
  # browser()
  p <- ncol(Rtilda)
  trace_R <- sum(diag(Rtilda))

  if (npc > p){
    stop("number of PC should be less than or equal to the dimension")
  }
  if (npc < p){
    firstneigs <- RSpectra::eigs_sym(Rtilda, k = npc, which = "LA")
    lambdas <- firstneigs$values
    rotmat <- firstneigs$vectors
    sigma2 <- (trace_R - sum(lambdas))/(p - npc)
    sig <- sqrt(sigma2)
    trans_res <- gdtrans(mapped_data, Fhats)
    latent_norm_ysim <- truncnormalsim_trans(lbounds = qnorm(trans_res$lower_bounds), ubounds = qnorm(trans_res$upper_bounds), inity = qnorm(trans_res$transformed_data), rotmat = rotmat, lambda = lambdas, sig = sig)
    scores <- latent_norm_ysim %*% rotmat
  } else {
    firstneigs <- eigen(Rtilda)
    lambdas <- firstneigs$values[-npc]
    rotmat <- firstneigs$vectors[,-npc]
    sigma2 <- firstneigs$values[npc]
    sig <- sqrt(sigma2)
    trans_res <- gdtrans(mapped_data)
    latent_norm_ysim <- truncnormalsim_trans(lbounds = qnorm(trans_res$lower_bounds), ubounds = qnorm(trans_res$upper_bounds), inity = qnorm(trans_res$transformed_data), rotmat = rotmat, lambda = lambdas, sig = sig)
    scores <- latent_norm_ysim %*% firstneigs$vectors
    rotmat <- firstneigs$vectors
  }

  return(list(lambdas = lambdas, rotmat = rotmat, scores = scores, latent_norm = latent_norm_ysim))

}

#' Perform PCA for discrete dataset using Gaussian copula
#'
#' @param x An \code{DGCFit} object or a discrete dataset.
#' @param npc Number of first PCs needed.
#'
#' @return An object of class \code{DGCPCA} with components:
#' \describe{
#' \item{fit}{An object of \code{DGCFit}. See \code{\link{DGCFit}}.}
#' \item{npc}{Number of first PCs needed.}
#' \item{lambdas}{Eigenvalues for the first \code{npc} PCs.}
#' \item{rotmat}{Rotation matrix.}
#' \item{scores}{The first \code{npc} surrogate PC scores.}
#' \item{latent_norm}{The simulated unobserved multivariate normal random vectors used to generate the surrogate PC scores.}
#' }
#' @export
DGCPCA <- function(x, npc, ...){
  UseMethod("DGCPCA")
}


#' @describeIn DGCPCA S3 method for \code{DGCFit} object
#' @export
DGCPCA.DGCFit <- function(x,npc = 5, ...){
  res <- list()
  res$fit <- x
  res$npc <- npc
  res <- c(res, firstnpc(npc = npc, Rtilda = x$corr, mapped_data = x$mapped_data, Fhats = x$cdfs))
  class(res) <- "DGCPCA"
  return(res)
}

#' @describeIn DGCPCA S3 method for default
#'
#' @param x A discrete dataset. Assumed to be all numerical-discrete-valued if \code{maps = FALSE}.
#' @param npc Number of first PCs needed.
#' @param maps Logical or a list of vectors of length \code{ncol(x)}. If FALSE, \code{x} must be a mapped dataset (i.e. all discrete values are numerical). If TRUE, all values will be mapped to discrete values in [0,1] using \code{sort(unique(x[,j]))}. If a list of vectors of length \code{ncol(x)}, values will be mapped to discrete values in [0,1] accoording to the provided list. See \code{\link{discrete_mapping}}.
#' @param ... Parameters for \code{\link{DGCFit}}.
#' @export
DGCPCA.default <- function(x, npc = 5, maps = TRUE, ...){
  res <- list()
  if (isFALSE(maps)) {
    mapped_data <- x
  } else {
    if (isTRUE(maps)) {
      maps <- lapply(x, function(x) sort(unique(x)))
      mapped_data <- discrete_mapping(x, maps)
    } else if (is.list(maps)) {
      mapped_data <- discrete_mapping(x, maps)
    } else {
      stop("maps must be logical or a list of vectors.")
    }
  }

  fit <- DGCFit(mapped_data, ...)
  return(DGCPCA.DGCFit(fit, npc = npc, ...))
}



