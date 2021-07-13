fast_kendall <-function (x, y = NULL)
{
  cor = FALSE
  if (is.null(y)) {
    n <- nrow(x)
    if (!is.matrix(x) && !is.data.frame(x))
      stop("x must be either numeric vector, matrix or data.frame.")
    {
      p <- ncol(x)
      dn <- colnames(x)
      ret <- diag(p)
      dimnames(ret) <- list(dn, dn)
      for (i in 1:p) {
        if (i == p)
          return(ret)
        ord <- order(x[, i])
        cur.x <- x[ord, i]
        for (j in (i + 1):p) ret[i, j] <- ret[j, i] <- pcaPP:::.cor.fk.2d(cur.x,
                                                                          x[ord, j], cor)/(n*(n-1))
      }
      }
  }
  else {
    if (length(x) != length(y))
      stop("x and y must have same length.")
    n <- length(x)
    ord <- order(x)
    return(pcaPP:::.cor.fk.2d(x[ord], y[ord], cor)/(n*(n-1)))
  }
}
