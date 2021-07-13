


drecover <- function(scores, rotmat,  Fhats, k_used){
  # browser()
  if (k_used < 0){
    stop("k_used need to be greater than 0")
  }
  if (k_used > 0){
    transformed_npc <- pnorm(scores[,1:k_used] %*% t(rotmat[, 1:k_used]))
  } else {
    transformed_npc <- array(0, dim = c(nrow(scores), nrow(rotmat)))
  }
  return(dtrans_inv(transformed_npc, Fhats))
}




#' Get the recovery accuracies
#'
#' Get the recovery accuracies using 0 to \code{dgcpca$npc} PCs.
#'
#' @param dgcpca An object of class \code{DGCPCA}.
#' @return A vector of recovery accuracies.
#' @export
get_recovery_accuray <- function(dgcpca){
  # browser()
  # npc <- gcpca$npc
  # Rtilda <- gcpca$fit$corr
  # mapped_data <- gcpca$fit$cdfs
  # first_npc_res <- firstnpc(npc, Rtilda, mapped_data)
  npc <- dgcpca$npc
  r_accuracy <- rep(0, npc + 1)
  for (k in 0:npc){
    recovered_data <- data_recover.DGCPCA(dgcpca, npc_used = k)
    r_accuracy[k+1] <- mean(recovered_data == dgcpca$fit$mapped_data)
  }
  return(r_accuracy)
}


#' Recover the mapped data (i.e. all discrete features are numerical) using first several PC scores
#'
#' @param x An \code{DGCPCA} object.
#' @param npc_used The number of PCs used for data recovery, no larger than \code{x$npc}.
#' @export
data_recover <- function(x, npc_used, ...){
  UseMethod("data_recover")
}

#' @describeIn data_recover S3 method for \code{DGCPCA}
#' @export
data_recover.DGCPCA <- function(x, npc_used, ...){
  if (missing(npc_used)) {
    npc_used <- x$npc
  }
  return(drecover(scores = x$scores, rotmat = x$rotmat, Fhats = x$fit$cdfs, k_used = npc_used))
}
