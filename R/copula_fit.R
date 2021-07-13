#' @import stats
#' @import foreach
NULL


#' Map every column of the original discrete data to new data with discrete values in [0,1]
#'
#' @param data The original discrete data
#' @param maps A list of the same length of ncol(data), each item is a vector containing the original image of the map. The vector should represent the ascending order of the discrete values in the original data, as the discrete values will be mapped to numbers in [0,1] in a ascending order.
#' @return A new dataset. All values are in [0,1].
#' @export
discrete_mapping <- function(data, maps){
  # maps is a list of the same length of ncol, each item contains the original image of the map
  # map every column to discrete values in [0,1]
  data_mapped <- array(dim = dim(data))
  for (j in 1:ncol(data)){
    #browser()
    mapped_j <- plyr::mapvalues(data[,j], from = maps[[j]], to = seq(0, 1, length.out = length(maps[[j]])))
    if (is.factor(mapped_j)){
      data_mapped[,j] <- as.numeric(levels(mapped_j))[mapped_j]
    }
    else {
      data_mapped[,j] <- as.numeric(mapped_j)
    }
  }
  return(data_mapped)
}

Fhat <- function(data){
  # returns a list of emperical cdfs. The list is of the same length of ncol(data)
  ecdfs <- list()

  for (j in 1:ncol(data)){
    ecdfs[[j]] <- ecdf(data[,j])
  }
  return(ecdfs)
}

equantile <- function(supp, probs){
  supp;
  probs;
  eqt <- function(x) stepfun(probs, supp, right = T)(x)
  return(eqt)
}

Fhat_inv <- function(Fhats){
  # returns a list of quantile functions. Inverse of Fhat.
  equantiles <- list()
  for (j in 1:length(Fhats)){
    supp <- knots(Fhats[[j]])
    probs <- Fhats[[j]](supp)
    equantiles[[j]] <- equantile(supp = supp, probs = probs[-length(probs)])
  }
  return(equantiles)
}

# library(EnvStats)

epmf <- function(supp, probs){
  supp;
  probs;
  epm <- function(x){
    y <- probs[match(x, supp)]
    y[is.na(y)] <- 0
    return(y)
  }
  return(epm)
}

pmf_hat <- function(Fhats){
  # returns a list of emperical pmf functions
  epmfs <- list()
  for (j in 1:length(Fhats)){
    #browser()
    supp <- knots(Fhats[[j]])
    probs <- diff(c(0, Fhats[[j]](supp)))
    epmfs[[j]] <- epmf(supp, probs)
  }
  return(epmfs)
}


#' Generalized distributional transformation of discrete data
#'
#' @param mapped_data A discrete dataset whose discrete values are all numeric or have been mapped to numbers.
#' @param Fhat A list of estimated marginal cdfs of the dataset. Each item in the list is a \code{stepfun} object and the support of the cdfs (\code{knots(Fhats[[j]])}) should be the superset of the unique discrete values of the corresponding feature (\code{unique(mapped_data[,j])}). If not specified, \code{ecdf} will be used to estimate the marginal cdfs.

gdtrans <- function(mapped_data, Fhats){
  # generalized distributional transformation
  transformed_data <- array(dim = dim(mapped_data))
  lower_bounds <- array(dim = dim(mapped_data))
  upper_bounds <- array(dim = dim(mapped_data))
  # Fhat_mapped_data <- Fhat(mapped_data)
  if (missing(Fhats)) {
    Fhats <- Fhat(mapped_data)
  }
  pmf_hats <- pmf_hat(Fhats)
  V <- matrix(runif(n = nrow(mapped_data)*ncol(mapped_data)), ncol = ncol(mapped_data))
  for(j in 1:ncol(mapped_data)){
    upper_bounds[,j] <- Fhats[[j]](mapped_data[,j])
    lower_bounds[,j] <- upper_bounds[,j] - pmf_hats[[j]](mapped_data[,j])
    transformed_data[,j] <- lower_bounds[,j] + V[,j] * (upper_bounds[,j] - lower_bounds[,j])
  }
  return(list(transformed_data = transformed_data, upper_bounds = upper_bounds, lower_bounds = lower_bounds))
}


dtrans_inv <- function(transformed_data, Fhats){
  Fhat_inv_mapped_data <- Fhat_inv(Fhats)
  data_back <- array(dim = dim(transformed_data))
  for (j in 1:ncol(transformed_data)){
    data_back[,j] <- Fhat_inv_mapped_data[[j]](transformed_data[,j])
  }
  return(data_back)
}
